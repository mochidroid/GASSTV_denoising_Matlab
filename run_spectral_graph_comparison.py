"""スペクトルグラフの作り方を変えて GASSTV 本体を回し, ノイズ除去性能を比べる (パラメータ sweep 対応).

空間グラフ・SSTV・ノイズモデル・epsilon/alpha/beta の式は run_gasstv.py (= exp_denoising_gstd_Base.m) と同じ.
複数の値を渡した引数はすべて全組み合わせで実行する.

例:
    # 何回実行されるかだけ確認
    uv run run_spectral_graph_comparison.py --graphs baseline window_h2 --omega 0 200 2000 --rho 0.95 0.98 --dry_run

    # モデル側のパラメータを振る
    uv run run_spectral_graph_comparison.py --graphs baseline --omega 0 200 2000 --rho 0.95 0.98 1.02 \\
        --lambda_rho_sp 0.8 1.0 --sigma_sp 90 med --noise_idx 1 9

    # グラフ構築側のパラメータを振る (そのキーを持つグラフにだけ効く)
    uv run run_spectral_graph_comparison.py --graphs window_h2 --omega 200 \\
        --graph_param h=1,2,4 k=10,30 signed=true,false

    # 集計 (既定は image, noise, graph ごとの最良. --group_by で切り口を変える)
    uv run run_spectral_graph_comparison.py --summary_only
    uv run run_spectral_graph_comparison.py --summary_only --group_by image graph omega rho

結果:
    result/spectral_graph_comparison/<image>/<noise>/<graph>/<params>.mat
    result/spectral_graph_comparison/summary.csv   (1 行 1 実行. 全パラメータを記録)
"""

from __future__ import annotations

import argparse
import csv
import json
import inspect
import itertools
import time
from pathlib import Path

import numpy as np
import torch

from run_gasstv import NOISE_CONDITIONS, get_observation
from src.data import DEFAULT_DATA_DIR, noise_radii, obsv_dirname
from src.io import save_mat
from src.metrics import calc_mpsnr, calc_mssim, calc_sam
from src.solver import SolverOptions, gasstv_condatvu, select_noise_model
from src.spectral_graphs import graph_abs_signed, graph_anchored, graph_coherence

ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# グラフのプリセット: 名前 -> (spec_graph, 既定のグラフ構築パラメータ)
# ---------------------------------------------------------------------------
GRAPH_PRESETS = {
    # 現行法 (run_gasstv.py の noisy プリセットと同じ)
    "baseline": ("baseline", dict(num_segments=4, sigma_l="med", order_filt=5, k_lap=1000)),
    "baseline_k10": ("baseline", dict(num_segments=4, sigma_l="med", order_filt=5, k_lap=10)),
    # 案 1: |delta| の類似度 + 符号付きラプラシアン
    "abs_h0": ("abs", dict(feature="log", h=0, k=30, knn="mutual")),
    # 案 2: |delta| の窓 (半幅 h) の類似度 + 符号付きラプラシアン
    "window_h1": ("abs", dict(feature="log", h=1, k=30, knn="mutual")),
    "window_h2": ("abs", dict(feature="log", h=2, k=30, knn="mutual")),
    "window_h2_unsigned": ("abs", dict(feature="log", h=2, k=30, knn="mutual", signed=False)),
    # 案 3: 低空間周波数コヒーレンス
    "coherence": ("coherence", dict(min_sep=2, k=30)),
    # sweep_all.py 用の名前 (中身は上と同じ. --graph_param で上書きして使う)
    "p1_abs": ("abs", dict(feature="log", h=0, k=30, knn="mutual")),
    "p2_window": ("abs", dict(feature="log", h=2, k=30, knn="mutual")),
    "p3_coherence": ("coherence", dict(min_sep=2, k=30)),
    # 提案法: 代表スペクトルにつないだスペクトルグラフ (src/spectral_graphs.py の graph_anchored)
    "anchored": ("anchored", dict(edges="delta", anchor="uniform", beta_q=1.0, num_segments=4, seg="intensity")),
    "anchor_only": ("anchored", dict(edges="none", anchor="uniform", num_segments=4, seg="intensity")),
    "residual_edges": ("anchored", dict(edges="delta", anchor="none", num_segments=4, seg="intensity")),
}

# --graph_param で振れるキー (spec_graph ごと)
GRAPH_KEYS = {
    "baseline": {"num_segments", "sigma_l", "order_filt", "k_lap"},
    "abs": set(inspect.signature(graph_abs_signed).parameters) - {"X", "name"},
    "coherence": set(inspect.signature(graph_coherence).parameters) - {"X", "name"},
    "anchored": set(inspect.signature(graph_anchored).parameters) - {"X", "name"},
}

SUMMARY_FIELDS = ["image", "noise", "graph", "guide", "omega", "lambda2_mode", "lambda_rho_sp", "sigma_sp",
                  "rho", "graph_kw", "graph_kw_json", "mpsnr", "mssim", "sam", "iter", "time_s",
                  "best_mpsnr_during", "best_iter", "path"]


# ---------------------------------------------------------------------------
def parse_scalar(s: str):
    t = s.strip()
    if t.lower() in ("true", "false"):
        return t.lower() == "true"
    if t.lower() in ("none", "null"):
        return None
    if t in ("med", "90"):  # sigma のモード指定は文字列のまま
        return t
    try:
        v = float(t)
        return int(v) if v.is_integer() and "." not in t and "e" not in t.lower() else v
    except ValueError:
        return t


def parse_graph_params(items: list[str]) -> dict[str, list]:
    """["h=1,2", "k=10,30"] -> {"h": [1, 2], "k": [10, 30]}."""
    out = {}
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"--graph_param は key=v1,v2,... の形で指定してください: {it}")
        k, v = it.split("=", 1)
        out[k.strip()] = [parse_scalar(x) for x in v.split(",")]
    return out


def graph_variants(graph: str, sweep: dict[str, list]) -> list[tuple[dict, dict]]:
    """プリセットに --graph_param の組み合わせを重ねる. 戻り値: [(全パラメータ, 上書きした分)]."""
    spec, base = GRAPH_PRESETS[graph]
    keys = [k for k in sweep if k in GRAPH_KEYS[spec]]
    if not keys:
        return [(dict(base), {})]
    out = []
    for combo in itertools.product(*[sweep[k] for k in keys]):
        ov = dict(zip(keys, combo))
        out.append(({**base, **ov}, ov))
    return out


def kw_text(d: dict) -> str:
    return "_".join(f"{k}{v}" for k, v in sorted(d.items()))


def savetext(omega, mode, lrs, ssp, rho, override, args) -> str:
    s = f"om{omega:g}_{mode[0]}_l{lrs:g}_sig{ssp}_r{rho:.2f}_stop1e-{args.stopcri_idx}"
    if override:
        s += "_" + kw_text(override)
    if args.maxiter != 20000:
        s += f"_it{args.maxiter}"
    return s


# ---------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_argument_group("実験条件")
    g.add_argument("--images", nargs="+", default=["JasperRidge", "PaviaU"])
    g.add_argument("--noise_idx", nargs="+", type=int, default=[9],
                   help=f"run_gasstv.NOISE_CONDITIONS の番号 (1..{len(NOISE_CONDITIONS)})")
    g.add_argument("--graphs", nargs="+", default=list(GRAPH_PRESETS), choices=list(GRAPH_PRESETS))
    g.add_argument("--guide", nargs="+", default=["noisy"], choices=["noisy", "oracle"],
                   help="グラフを作る画像. oracle は GT から (空間グラフ・スペクトルグラフとも)")

    g = p.add_argument_group("モデルのパラメータ (複数値で全組み合わせ)")
    g.add_argument("--omega", nargs="+", type=float, default=[20, 200, 2000],
                   help="GLR の重み. normalized のとき L/lam_max に掛かる. raw のときは lambda2 そのもの")
    g.add_argument("--lambda2_mode", choices=["normalized", "raw"], default="normalized")
    g.add_argument("--lambda_rho_sp", nargs="+", type=float, default=[0.8],
                   help="空間グラフ L1 ボールの半径の倍率 (lambda_sp = 倍率 x ||W_sp D_sp(参照画像)||_1)")
    g.add_argument("--sigma_sp", nargs="+", default=["90"], help='空間グラフの sigma: "90" | "med" | 数値')
    g.add_argument("--rho", nargs="+", type=float, default=[0.98],
                   help="epsilon (L2 ボール), alpha (スパース), beta (ストライプ) の半径の倍率")
    g.add_argument("--maxiter", type=int, default=20000)
    g.add_argument("--stopcri_idx", type=int, default=5, help="stopcri = 10^-idx")
    g.add_argument("--lambda_sp_ref", choices=["clean", "noisy", "guide"], default="clean")

    g = p.add_argument_group("スペクトルグラフ構築のパラメータ")
    g.add_argument("--graph_param", nargs="+", default=[], metavar="KEY=V1,V2",
                   help="例: h=1,2,4 k=10,30 signed=true,false tau_q=0.3,0.5 num_segments=1,4. "
                        "そのキーを持つグラフにだけ適用 (持たないグラフは既定値で 1 回だけ実行)")

    g = p.add_argument_group("実行環境・入出力")
    g.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    g.add_argument("--fast_l1", choices=["auto", "on", "off"], default="auto",
                   help="L1 ボール射影を Michelot 法にする (結果は同一. auto: CPU のときだけ)")
    g.add_argument("--compile", action="store_true", help="差分作用素などを torch.compile で融合 (結果は同一)")
    g.add_argument("--track_every", type=int, default=100, help="反復中の MPSNR を記録する間隔 (0 で記録しない)")
    g.add_argument("--seed", type=int, default=0)
    g.add_argument("--data_dir", default=DEFAULT_DATA_DIR)
    g.add_argument("--obsv", choices=["auto", "dataset", "python"], default="dataset")
    g.add_argument("--dataset_dir", default=str(ROOT / "dataset"))
    g.add_argument("--result_dir", default=str(ROOT / "result" / "spectral_graph_comparison"))
    g.add_argument("--skip_existing", action="store_true")
    g.add_argument("--dry_run", action="store_true", help="実行する組み合わせと回数だけ表示")
    g.add_argument("--quiet", action="store_true")

    g = p.add_argument_group("集計")
    g.add_argument("--summary_only", action="store_true")
    g.add_argument("--group_by", nargs="+", default=["image", "noise", "graph"],
                   help=f"集計の切り口 (summary.csv の列名: {', '.join(SUMMARY_FIELDS[:10])})")
    return p


# ---------------------------------------------------------------------------
def print_summary(csv_path: Path, group_by: list[str]):
    if not csv_path.exists():
        print("(no summary yet)")
        return
    rows = list(csv.DictReader(open(csv_path, newline="", encoding="utf-8")))
    best = {}
    for r in rows:
        k = tuple(r.get(c, "") for c in group_by)
        if k not in best or float(r["mpsnr"]) > float(best[k]["mpsnr"]):
            best[k] = r
    show = ["omega", "lambda_rho_sp", "sigma_sp", "rho", "graph_kw"]
    show = [c for c in show if c not in group_by]
    print("\n" + " | ".join(f"{c:>12s}" for c in group_by) + " |  MPSNR   MSSIM    SAM   iter | " + ", ".join(show))
    for k in sorted(best):
        r = best[k]
        rest = ", ".join(f"{c}={r.get(c, '')}" for c in show)
        print(" | ".join(f"{v:>12s}" for v in k)
              + f" | {float(r['mpsnr']):6.3f} {float(r['mssim']):6.4f} {float(r['sam']):6.3f} {int(r['iter']):6d} | {rest}")
    print(f"({len(rows)} runs)")


def main():
    args = build_parser().parse_args()
    out_root = Path(args.result_dir)
    csv_path = out_root / "summary.csv"
    if args.summary_only:
        print_summary(csv_path, args.group_by)
        return

    sweep = parse_graph_params(args.graph_param)
    unknown = [k for k in sweep if not any(k in GRAPH_KEYS[GRAPH_PRESETS[g][0]] for g in args.graphs)]
    if unknown:
        raise SystemExit(f"--graph_param のキー {unknown} はどのグラフにもありません. 使えるキー: "
                         + "; ".join(f"{s}: {sorted(v)}" for s, v in GRAPH_KEYS.items()))

    # 実行リスト
    jobs, seen_zero = [], set()
    for nidx, image, guide, graph in itertools.product(args.noise_idx, args.images, args.guide, args.graphs):
        for (kw, ov), omega, lrs, ssp, rho in itertools.product(
                graph_variants(graph, sweep), args.omega, args.lambda_rho_sp, args.sigma_sp, args.rho):
            if omega == 0:  # グラフ項なしはグラフによらず同じ結果なので 1 回だけ (最初のグラフで) 回す
                zkey = (nidx, image, guide, lrs, ssp, rho)
                if zkey in seen_zero:
                    continue
                seen_zero.add(zkey)
            jobs.append(dict(nidx=nidx, image=image, guide=guide, graph=graph, kw=kw, ov=ov,
                             omega=omega, lrs=lrs, ssp=ssp, rho=rho))
    print(f"total runs: {len(jobs)}")
    if args.dry_run:
        for j in jobs:
            print(f"  {j['image']:12s} n{j['nidx']} {j['guide']:6s} {j['graph']:20s} omega={j['omega']:g} "
                  f"lrs={j['lrs']:g} sigma_sp={j['ssp']} rho={j['rho']:g} {kw_text(j['ov'])}")
        return

    if args.fast_l1 == "on" or (args.fast_l1 == "auto" and args.device == "cpu"):
        from src import fast_ops
        fast_ops.enable()
    if args.compile:
        from src import fast_ops
        fast_ops.enable_compile()

    stopcri = 10.0 ** -args.stopcri_idx
    sig = lambda v: v if v in ("90", "med") else float(v)  # noqa: E731
    out_root.mkdir(parents=True, exist_ok=True)
    obs_cache = {}

    for n, j in enumerate(jobs, 1):
        key = (j["nidx"], j["image"])
        if key not in obs_cache:
            obs_cache.clear()
            cond = NOISE_CONDITIONS[j["nidx"] - 1]
            deg0 = dict(zip(["gaussian_sigma", "sparse_rate", "stripe_rate", "stripe_intensity", "deadline_rate"],
                            map(float, cond)))
            obs_cache[key] = get_observation(j["image"], deg0, args)
        HSI_clean, HSI_noisy, hsi, deg = obs_cache[key]
        noise_model = select_noise_model(deg)
        noise_name = obsv_dirname(deg)

        spec, _ = GRAPH_PRESETS[j["graph"]]
        st = savetext(j["omega"], args.lambda2_mode, j["lrs"], j["ssp"], j["rho"], j["ov"], args)
        gdir = j["graph"] + ("_OG" if j["guide"] == "oracle" else "")
        path = out_root / j["image"] / noise_name / gdir / f"{st}.mat"
        if args.skip_existing and path.exists():
            print(f"[{n}/{len(jobs)}] skip {path}")
            continue

        rho = j["rho"]
        params = dict(
            **noise_radii(rho, deg, hsi, HSI_noisy),  # deadline はスパースの項に含める
            maxiter=args.maxiter, stopcri=stopcri, rho_radius=rho,
            lambda2=j["omega"], lambda2_mode=args.lambda2_mode,
            lambda_rho_sp=j["lrs"], sigma_sp=sig(j["ssp"]), spec_graph=spec,
            sigma_l="med", num_segments=4, order_filt=5, k_lap=1000,
        )
        if spec == "baseline":
            params.update(j["kw"])
        else:
            params["spec_graph_kw"] = j["kw"]

        print(f"\n=== [{n}/{len(jobs)}] {j['image']} {noise_name} | {gdir} | omega={j['omega']:g} "
              f"lrs={j['lrs']:g} sigma_sp={j['ssp']} rho={rho:g} {kw_text(j['ov'])} ===", flush=True)
        opt = SolverOptions(oracle=(j["guide"] == "oracle"), lambda_sp_ref=args.lambda_sp_ref,
                            device=args.device, seed=args.seed, track_every=args.track_every,
                            disp_iters=tuple(range(1000, args.maxiter + 1, 1000)), verbose=not args.quiet)
        t0 = time.time()
        U, _, other = gasstv_condatvu(HSI_clean, HSI_noisy, params, noise_model, opt)
        elapsed = time.time() - t0
        mp, ms, sa = calc_mpsnr(U, HSI_clean), calc_mssim(U, HSI_clean), calc_sam(U, HSI_clean)
        traj = other["move_mpsnr"]
        ok = np.isfinite(traj).any()
        best_i = int(np.nanargmax(traj)) + 1 if ok else -1
        best_v = float(np.nanmax(traj)) if ok else np.nan
        print(f"-> MPSNR {mp:.3f}  MSSIM {ms:.4f}  SAM {sa:.3f}  iter {other['iter']}  time {elapsed:.0f}s  "
              f"(best during run {best_v:.3f} @ {best_i})", flush=True)

        save_mat(path, dict(
            HSI_restored=U, params={k: v for k, v in params.items() if k != "spec_graph_kw"},
            graph=j["graph"], graph_kw=str(j["kw"]), guide=j["guide"],
            val_mpsnr=mp, val_mssim=ms, val_sam=sa, iter=other["iter"],
            move_mpsnr=traj, converge_rate_U=other["converge_rate_U"], lam_max=other["lam_max"],
        ))
        new_csv = not csv_path.exists()
        with open(csv_path, "a", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
            if new_csv:
                w.writeheader()
            w.writerow(dict(image=j["image"], noise=noise_name, graph=j["graph"], guide=j["guide"],
                            omega=j["omega"], lambda2_mode=args.lambda2_mode, lambda_rho_sp=j["lrs"],
                            sigma_sp=j["ssp"], rho=rho, graph_kw=kw_text(j["kw"]),
                            graph_kw_json=json.dumps(j["kw"]),
                            mpsnr=mp, mssim=ms, sam=sa, iter=other["iter"], time_s=round(elapsed),
                            best_mpsnr_during=best_v, best_iter=best_i, path=str(path)))
    print_summary(csv_path, args.group_by)


if __name__ == "__main__":
    main()
