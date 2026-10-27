"""ノイズ条件 1-9 で, スペクトルグラフの候補を比べる (rho, lambda_rho_sp は固定, omega は条件ごとに振る).

    uv run sweep_noise_conditions.py --dry_run          # 実行回数
    uv run sweep_noise_conditions.py                    # 実行 (止めても再実行で続きから)
    uv run sweep_noise_conditions.py --report           # 集計 (result/spectral_graph_comparison/noise_conditions/report.md)
    uv run sweep_noise_conditions.py --noise_idx 1 2 3  # 一部の条件だけ

比べるもの (各条件・各画像で omega = 10, 30, 100, 300, 1000 を回し, 最良の omega の値で比較する):
  omega0      : スペクトルグラフ項なし (omega = 0)
  baseline    : 現行グラフ
  p2_mknn     : 案 2, mutual kNN (k = 60), h = 1, 符号は厳密版           ... 最良値が最も良かった設定
  p2_mknn_ns  : 案 2, mutual kNN (k = 60), h = 1, 符号なし               ... 符号の効果を見る比較
  p2_eps      : 案 2, eps 近傍 (eps_c = 0.5), h = 2, 符号は厳密版        ... omega への頑健性が最も高かった設定

ノイズ条件の番号は run_gasstv.NOISE_CONDITIONS (1: g0.1 / 2: g0.05 ps0.05 / 3: g0.1 ps0.05 / 4: g0.05 pt0.05 /
5: g0.1 pt0.05 / 6: g0.05 ps0.05 pt0.05 / 7: g0.1 ps0.05 pt0.05 / 8: g0.1 pd0.01 / 9: g0.1 ps0.05 pt0.05 pd0.01).
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / "run_spectral_graph_comparison.py"
OUT = ROOT / "result" / "spectral_graph_comparison" / "noise_conditions"
TIE_DB = 0.01  # 現行グラフとの差がこれ以内なら引き分け

# 名前 -> (グラフ, --graph_param)
VARIANTS = {
    "baseline": ("baseline", []),
    "p2_mknn": ("p2_window", ["knn=mutual", "k=60", "h=1", "sign_mode=exact"]),
    "p2_mknn_ns": ("p2_window", ["knn=mutual", "k=60", "h=1", "sign_mode=none"]),
    "p2_eps": ("p2_window", ["knn=eps", "eps_c=0.5", "k=0", "h=2", "sign_mode=exact"]),
}
LABEL = {
    "omega0": "スペクトルグラフ項なし",
    "baseline": "現行グラフ",
    "p2_mknn": "案2 mkNN 符号あり",
    "p2_mknn_ns": "案2 mkNN 符号なし",
    "p2_eps": "案2 eps 符号あり",
}


def variant_of(r: dict) -> str | None:
    """summary.csv の 1 行がどの比較対象か."""
    if float(r["omega"]) == 0:
        return "omega0"
    if r["graph"] == "baseline":
        return "baseline"
    kw = json.loads(r.get("graph_kw_json") or "{}")
    sm = kw.get("sign_mode")
    sm = "none" if sm in (None, "None") else sm
    if kw.get("knn") == "mutual" and kw.get("k") == 60 and kw.get("h") == 1:
        return "p2_mknn" if sm == "exact" else ("p2_mknn_ns" if sm == "none" else None)
    if kw.get("knn") == "eps" and abs(float(kw.get("eps_c", -1)) - 0.5) < 1e-9 and kw.get("h") == 2 and sm == "exact":
        return "p2_eps"
    return None


def report(out: Path):
    csv_path = out / "summary.csv"
    if not csv_path.exists():
        print("(結果がまだありません)")
        return
    rows = []
    for r in csv.DictReader(open(csv_path, newline="", encoding="utf-8")):
        v = variant_of(r)
        if v is None:
            continue
        r["variant"] = v
        for k in ("omega", "mpsnr", "mssim", "sam"):
            r[k] = float(r[k])
        rows.append(r)

    order = ["omega0", "baseline", "p2_mknn", "p2_mknn_ns", "p2_eps"]
    lines = ["# ノイズ条件ごとの比較", "",
             "各条件・各画像で, 比較対象ごとに omega を振った中で MPSNR が最良の値。",
             "Δ は現行グラフとの差 (MPSNR は大きいほど良い, SAM は小さいほど良い)。",
             "劣化 = 最良値 - omega=1000 の値 (omega を強めすぎたときの落ち込み. 小さいほど頑健)。", ""]
    summary = {v: dict(d_mpsnr=[], d_sam=[], drop=[], win=0, tie=0, loss=0) for v in order}
    for image in sorted({r["image"] for r in rows}):
        lines += [f"## {image}", "",
                  "| 条件 | " + " | ".join(LABEL[v] for v in order) + " |",
                  "|---|" + "---|" * len(order)]
        sam_lines = ["", "SAM / 劣化 (ω=1000):", "", "| 条件 | " + " | ".join(LABEL[v] for v in order) + " |",
                     "|---|" + "---|" * len(order)]
        for noise in sorted({r["noise"] for r in rows if r["image"] == image}):
            sub = [r for r in rows if r["image"] == image and r["noise"] == noise]
            best = {}
            drop = {}
            for v in order:
                rs = [r for r in sub if r["variant"] == v]
                if not rs:
                    continue
                best[v] = max(rs, key=lambda r: r["mpsnr"])
                at1000 = [r for r in rs if r["omega"] == 1000]
                drop[v] = best[v]["mpsnr"] - at1000[0]["mpsnr"] if at1000 else np.nan
            if "baseline" not in best:
                continue
            b = best["baseline"]
            top = max(best[v]["mpsnr"] for v in best)
            cells, scells = [], []
            for v in order:
                if v not in best:
                    cells.append("–")
                    scells.append("–")
                    continue
                r = best[v]
                mark = "**" if abs(r["mpsnr"] - top) < 1e-9 else ""
                if v == "baseline":
                    cells.append(f"{mark}{r['mpsnr']:.3f}{mark} (ω={r['omega']:g})")
                else:
                    cells.append(f"{mark}{r['mpsnr']:.3f}{mark} ({r['mpsnr'] - b['mpsnr']:+.3f}, ω={r['omega']:g})")
                d = drop.get(v, np.nan)
                scells.append(f"{r['sam']:.3f} / {d:.2f}" if v != "omega0" else f"{r['sam']:.3f}")
                if v not in ("omega0", "baseline"):
                    s = summary[v]
                    s["d_mpsnr"].append(r["mpsnr"] - b["mpsnr"])
                    s["d_sam"].append(r["sam"] - b["sam"])
                    dd = r["mpsnr"] - b["mpsnr"]
                    s["win" if dd > TIE_DB else ("loss" if dd < -TIE_DB else "tie")] += 1
                if v != "omega0" and np.isfinite(d):
                    summary[v]["drop"].append(d)
                if v == "omega0":
                    summary["omega0"]["d_mpsnr"].append(r["mpsnr"] - b["mpsnr"])
            lines.append(f"| {noise} | " + " | ".join(cells) + " |")
            sam_lines.append(f"| {noise} | " + " | ".join(scells) + " |")
        lines += sam_lines + [""]

    lines += ["## 全条件・全画像のまとめ", "",
              f"| 比較対象 | 現行グラフとの MPSNR 差 (平均 / 最小 / 最大) | 勝ち / 引き分け / 負け (±{TIE_DB} dB) | SAM 差 (平均) | 劣化 ω=1000 (平均) |",
              "|---|---|---|---|---|"]
    for v in order:
        s = summary[v]
        if v == "omega0":
            d = np.array(s["d_mpsnr"])
            if d.size:
                lines.append(f"| {LABEL[v]} | {d.mean():+.3f} / {d.min():+.3f} / {d.max():+.3f} | – | – | – |")
            continue
        dr = np.array(s["drop"])
        if v == "baseline":
            lines.append(f"| {LABEL[v]} | – | – | – | {dr.mean():.3f} |" if dr.size else f"| {LABEL[v]} | – | – | – | – |")
            continue
        d = np.array(s["d_mpsnr"])
        ds = np.array(s["d_sam"])
        if d.size:
            lines.append(f"| {LABEL[v]} | {d.mean():+.3f} / {d.min():+.3f} / {d.max():+.3f} | "
                         f"{s['win']} / {s['tie']} / {s['loss']} | {ds.mean():+.3f} | {dr.mean():.3f} |")
    text = "\n".join(lines)
    print(text)
    (out / "report.md").write_text(text, encoding="utf-8")
    print(f"\n-> {out / 'report.md'}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--images", nargs="+", default=["JasperRidge", "PaviaU"])
    p.add_argument("--noise_idx", nargs="+", default=[str(i) for i in range(1, 10)])
    p.add_argument("--variants", nargs="+", default=list(VARIANTS), choices=list(VARIANTS))
    p.add_argument("--rho", default="0.95")
    p.add_argument("--lambda_rho_sp", nargs="+", default=["0.8"])
    p.add_argument("--omega", nargs="+", default=["10", "30", "100", "300", "1000"])
    p.add_argument("--dry_run", action="store_true")
    p.add_argument("--report", action="store_true")
    a, passthrough = p.parse_known_args()

    if a.report:
        report(OUT)
        return

    common = [sys.executable, str(RUNNER), "--images", *a.images, "--noise_idx", *a.noise_idx,
              "--rho", a.rho, "--lambda_rho_sp", *a.lambda_rho_sp, "--result_dir", str(OUT),
              "--skip_existing", *passthrough]
    for v in a.variants:
        graph, gp = VARIANTS[v]
        omegas = (["0"] if v == "baseline" else []) + a.omega
        cmd = common + ["--graphs", graph, "--omega", *omegas] + (["--graph_param", *gp] if gp else [])
        cmd += ["--dry_run"] if a.dry_run else []
        print(f"\n##### {v} #####\n$ " + " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
    if not a.dry_run:
        report(OUT)


if __name__ == "__main__":
    main()
