"""各案 (現行法, 案 1, 案 2, 案 3) のパラメータを広く振る sweep. 3 段階を自動でつなぐ.

    A. モデルのパラメータ (omega, rho, lambda_rho_sp, sigma_sp) を, 各案の既定のグラフで全組み合わせ
    B. 各案・各画像について A の最良 (omega > 0) のモデルパラメータに固定し, グラフ構築のパラメータを全組み合わせ
    C. A, B を通じた最良のグラフ設定に固定し, omega を最良値のまわりで細かく振り直す

使い方 (Windows の PowerShell でもそのまま動く):
    uv run sweep_all.py --dry_run                 # 段階ごとの実行回数 (B, C は A の結果がないと仮の回数)
    uv run sweep_all.py                           # A -> B -> C を順に実行. 止めても再実行で続きから
    uv run sweep_all.py --levels A                # A だけ
    uv run sweep_all.py --proposals p2_window --images PaviaU
    uv run sweep_all.py --report                  # 結果の表 (result/spectral_graph_comparison/sweep_all/report.md にも保存)
    uv run sweep_all.py --profile quick --maxiter 50   # 動作確認用の小さいグリッド

グリッドは下の PROFILES を書き換えて調整する.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / "run_spectral_graph_comparison.py"

PROPOSAL_LABEL = {
    "baseline": "現行法",
    "p1_abs": "案1 絶対値 + 符号付き",
    "p2_window": "案2 窓 + 符号付き",
    "p3_coherence": "案3 コヒーレンス",
    "anchored": "提案法 (代表ノード + 辺)",
    "anchor_only": "代表ノードのみ",
    "residual_edges": "ずれに辺のみ",
}

# =============================================================================
# グリッド (ここを書き換える)
# =============================================================================
PROFILES = {
    "full": dict(
        # --- A: モデルのパラメータ (各案の既定グラフで全組み合わせ. omega = 0 はグラフによらないので 1 回だけ)
        model=dict(
            omega=[0, 10, 100, 1000, 10000],
            rho=[0.95, 0.98, 1.02],
            lambda_rho_sp=[0.8, 1.0],
            sigma_sp=["90", "med"],
        ),
        # --- B: グラフ構築のパラメータ (案ごと)
        graph=dict(
            baseline=dict(num_segments=[1, 4, 8], k_lap=[10, 30, 1000], order_filt=[1, 5]),
            p1_abs=dict(feature=["abs", "log"], k=[10, 30, 60], knn=["mutual", "union"],
                        tau_q=[0.3, 0.5, 0.7]),
            p2_window=dict(h=[1, 2, 4], k=[10, 30, 60], knn=["mutual", "union"], signed=[True, False]),
            p3_coherence=dict(f_max=[2, 4, 8], min_sep=[2, 4], k=[10, 30], prefilter=[0, 3]),
        ),
        # --- C: 最良の omega に掛ける倍率
        omega_factors=[0.3, 0.5, 2, 3],
    ),
    # rho を決めるための sweep (既存法の結果と並べて compare_methods.py で比べる). A だけ回せばよい
    "rho_select": dict(
        proposals=["baseline", "p2_window", "anchored", "anchor_only", "residual_edges"],
        model=dict(
            omega=[0, 10, 30, 100, 300, 1000],
            rho=[0.90, 0.95],
            lambda_rho_sp=[0.6, 0.7, 0.8, 0.9],
            sigma_sp=["90"],
        ),
        graph=dict(
            baseline=dict(num_segments=[1, 4, 8], k_lap=[10, 30, 1000], order_filt=[1, 5]),
            p2_window=dict(h=[1, 2], k=[30, 60, 100], knn=["mutual"]),
            anchored=dict(beta_q=[0.3, 1, 3], seg=["intensity", "spectral"], num_segments=[4, 8, 16],
                          rep=["median", "colmed"]),
            anchor_only=dict(seg=["intensity", "spectral"], num_segments=[4, 8, 16], rep=["median", "colmed"]),
            residual_edges=dict(seg=["intensity", "spectral"], num_segments=[4, 8, 16]),
        ),
        omega_factors=[0.5, 2],
    ),
    "quick": dict(
        model=dict(omega=[0, 100, 1000], rho=[0.98], lambda_rho_sp=[0.8], sigma_sp=["90"]),
        graph=dict(
            baseline=dict(num_segments=[1, 4]),
            p1_abs=dict(k=[10, 30]),
            p2_window=dict(h=[1, 2]),
            p3_coherence=dict(f_max=[2, 4]),
        ),
        omega_factors=[0.5, 2],
    ),
}


# =============================================================================
def fmt(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if v is None:
        return "none"
    return str(v)


def read_rows(csv_path: Path) -> list[dict]:
    if not csv_path.exists():
        return []
    return list(csv.DictReader(open(csv_path, newline="", encoding="utf-8")))


def best_row(rows, proposal, image, noise_name=None, require_positive_omega=True):
    cand = [r for r in rows if r["graph"] == proposal and r["image"] == image
            and (noise_name is None or r["noise"] == noise_name)
            and (not require_positive_omega or float(r["omega"]) > 0)]
    return max(cand, key=lambda r: float(r["mpsnr"])) if cand else None


def run(cmd: list[str], dry: bool):
    print("\n$ " + " ".join(cmd), flush=True)
    subprocess.run(cmd + (["--dry_run"] if dry else []), check=True)


def noise_names(rows, image):
    return sorted({r["noise"] for r in rows if r["image"] == image})


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profile", choices=list(PROFILES), default="full")
    p.add_argument("--levels", nargs="+", default=["A", "B", "C"], choices=["A", "B", "C"])
    p.add_argument("--proposals", nargs="+", default=None, choices=list(PROPOSAL_LABEL),
                   help="既定はプロファイルの proposals (無ければ現行法と案 1-3)")
    p.add_argument("--images", nargs="+", default=["JasperRidge", "PaviaU"])
    p.add_argument("--noise_idx", nargs="+", default=["9"])
    p.add_argument("--result_dir", default=None,
                   help="既定: result/spectral_graph_comparison/sweep_all (full) または sweep_<profile>")
    p.add_argument("--dry_run", action="store_true")
    p.add_argument("--report", action="store_true")
    a, passthrough = p.parse_known_args()  # --device, --maxiter, --compile などは run_spectral_graph_comparison.py へ

    cfg = PROFILES[a.profile]
    if a.proposals is None:
        a.proposals = cfg.get("proposals", ["baseline", "p1_abs", "p2_window", "p3_coherence"])
    if a.result_dir is None:
        name = "sweep_all" if a.profile == "full" else f"sweep_{a.profile}"
        a.result_dir = str(ROOT / "result" / "spectral_graph_comparison" / name)
    out = Path(a.result_dir)
    csv_path = out / "summary.csv"
    base = [sys.executable, str(RUNNER), "--result_dir", str(out), "--skip_existing", "--track_every", "100",
            "--noise_idx", *a.noise_idx, *passthrough]

    if a.report:
        report(read_rows(csv_path), a, out, passthrough)
        return

    m = cfg["model"]
    # ---------------- A ----------------
    if "A" in a.levels:
        print("\n########## A: モデルのパラメータ ##########")
        run(base + ["--images", *a.images, "--graphs", *a.proposals,
                    "--omega", *map(fmt, m["omega"]), "--rho", *map(fmt, m["rho"]),
                    "--lambda_rho_sp", *map(fmt, m["lambda_rho_sp"]), "--sigma_sp", *m["sigma_sp"]], a.dry_run)

    # ---------------- B ----------------
    if "B" in a.levels:
        print("\n########## B: グラフ構築のパラメータ ##########")
        rows = read_rows(csv_path)
        for image in a.images:
            for prop in a.proposals:
                for nn in (noise_names(rows, image) or [None]):
                    b = best_row(rows, prop, image, nn)
                    if b is None:
                        print(f"[B] {image} {prop}: A の結果がないので仮の値 (各グリッドの中央) で表示します"
                              if a.dry_run else f"[B] {image} {prop}: A の結果がないので飛ばします")
                        if not a.dry_run:
                            continue
                        b = dict(omega=m["omega"][len(m["omega"]) // 2], rho=m["rho"][0],
                                 lambda_rho_sp=m["lambda_rho_sp"][0], sigma_sp=m["sigma_sp"][0])
                    grid = cfg["graph"][prop]
                    gp = [f"{k}={','.join(fmt(v) for v in vals)}" for k, vals in grid.items()]
                    run(base + ["--images", image, "--graphs", prop,
                                "--omega", fmt(float(b["omega"])), "--rho", fmt(float(b["rho"])),
                                "--lambda_rho_sp", fmt(float(b["lambda_rho_sp"])), "--sigma_sp", b["sigma_sp"],
                                "--graph_param", *gp], a.dry_run)

    # ---------------- C ----------------
    if "C" in a.levels:
        print("\n########## C: 最良のグラフ設定で omega を振り直す ##########")
        rows = read_rows(csv_path)
        for image in a.images:
            for prop in a.proposals:
                for nn in (noise_names(rows, image) or [None]):
                    b = best_row(rows, prop, image, nn)
                    if b is None:
                        print(f"[C] {image} {prop}: 結果がないので飛ばします")
                        continue
                    om = float(b["omega"])
                    omegas = sorted({om * f for f in cfg["omega_factors"]})
                    kw = json.loads(b.get("graph_kw_json") or "{}")
                    gp = [f"{k}={fmt(v)}" for k, v in kw.items()]
                    run(base + ["--images", image, "--graphs", prop,
                                "--omega", *[f"{o:g}" for o in omegas], "--rho", fmt(float(b["rho"])),
                                "--lambda_rho_sp", fmt(float(b["lambda_rho_sp"])), "--sigma_sp", b["sigma_sp"]]
                        + (["--graph_param", *gp] if gp else []), a.dry_run)

    if not a.dry_run:
        report(read_rows(csv_path), a, out, passthrough)


def report(rows, a, out: Path, passthrough):
    if not rows:
        print("(結果がまだありません)")
        return
    maxiter = 20000
    if "--maxiter" in passthrough:
        maxiter = int(passthrough[passthrough.index("--maxiter") + 1])
    lines = ["# スペクトルグラフ sweep の結果", "",
             f"全 {len(rows)} 実行. 各案・各画像・各ノイズ条件で MPSNR が最良の設定 "
             f"(† は最大反復数に達して停止したもの).", ""]
    for image in sorted({r["image"] for r in rows}):
        for nn in noise_names(rows, image):
            sub = [r for r in rows if r["image"] == image and r["noise"] == nn]
            zero = [r for r in sub if float(r["omega"]) == 0]
            lines += [f"## {image} / {nn}", "",
                      "| 案 | MPSNR | MSSIM | SAM | omega | rho | λ_rho_sp | σ_sp | グラフ設定 | 反復 | 実行数 |",
                      "|---|---|---|---|---|---|---|---|---|---|---|"]
            if zero:
                z = max(zero, key=lambda r: float(r["mpsnr"]))
                lines.append(f"| (グラフ項なし ω=0) | {float(z['mpsnr']):.3f} | {float(z['mssim']):.4f} | "
                             f"{float(z['sam']):.3f} | 0 | {z['rho']} | {z['lambda_rho_sp']} | {z['sigma_sp']} | – | "
                             f"{z['iter']} | {len(zero)} |")
            for prop in PROPOSAL_LABEL:
                rs = [r for r in sub if r["graph"] == prop and float(r["omega"]) > 0]
                if not rs:
                    continue
                b = max(rs, key=lambda r: float(r["mpsnr"]))
                dag = "†" if int(b["iter"]) >= maxiter else ""
                lines.append(f"| {PROPOSAL_LABEL[prop]} | {float(b['mpsnr']):.3f} | {float(b['mssim']):.4f} | "
                             f"{float(b['sam']):.3f} | {float(b['omega']):g} | {b['rho']} | {b['lambda_rho_sp']} | "
                             f"{b['sigma_sp']} | {b['graph_kw']} | {b['iter']}{dag} | {len(rs)} |")
            # omega ごとの最良 (グラフ項が効いているかの確認)
            omegas = sorted({float(r["omega"]) for r in sub})
            lines += ["", "omega ごとの最良 MPSNR (他のパラメータは最良を選択):", "",
                      "| 案 | " + " | ".join(f"ω={o:g}" for o in omegas) + " |",
                      "|---|" + "---|" * len(omegas)]
            for prop in PROPOSAL_LABEL:
                if not any(r["graph"] == prop for r in sub):
                    continue
                vals = []
                for o in omegas:
                    rs = [r for r in sub if (r["graph"] == prop or o == 0) and float(r["omega"]) == o]
                    vals.append(f"{max(float(r['mpsnr']) for r in rs):.2f}" if rs else "–")
                lines.append(f"| {PROPOSAL_LABEL[prop]} | " + " | ".join(vals) + " |")
            lines.append("")
    text = "\n".join(lines)
    print(text)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.md").write_text(text, encoding="utf-8")
    print(f"\n-> {out / 'report.md'}")


if __name__ == "__main__":
    main()
