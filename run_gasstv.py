"""GASSTV (Condat-Vu 版) のデノイジング実験を回すスクリプト.

exp_denoising_gstd_Base.m の GASSTV_CondatVu / GASSTV_CondatVu_OG 部分に相当する.
下の「設定」を書き換えるか, コマンドライン引数で上書きして実行する.

例:
    uv run run_gasstv.py                                   # 既定設定 (非 Oracle, JasperRidge と PaviaU, ノイズ条件 1)
    uv run run_gasstv.py --graph oracle                    # Oracle Guide グラフ
    uv run run_gasstv.py --graph noisy oracle --images JasperRidge --noise_idx 1 7
    uv run run_gasstv.py --lambda_rho_sp 0.8 --lambda2 10 --sigma_sp 90 --maxiter 500
    uv run run_gasstv.py --noise 0.1 0.05 0.05 0.5 0       # ノイズ条件を直接指定 (g ps pt tint pd)
    uv run run_gasstv.py --dry_run                         # 実行内容の一覧だけ表示

結果:
    result/GASSTV_result/denoising_<image>/g<..>_ps<..>_pt<..>_pd<..>/<method>/<params>.mat      (全情報)
    result/GASSTV_comp_result/denoising_<image>/g<..>_ps<..>_pt<..>_pd<..>/<method>/<params>.mat (復元画像と指標のみ)
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import torch

from src.data import DEFAULT_DATA_DIR, generate_obsv, load_dataset, load_hsi, obsv_dirname
from src.io import save_mat
from src.metrics import calc_mpsnr, calc_mssim, calc_sam, psnr_per_band, ssim_per_band
from src.solver import SolverOptions, gasstv_condatvu, select_noise_model

ROOT = Path(__file__).resolve().parent
DATASET_DIR = ROOT / "dataset"  # Matlab/Generate_dataset_denoising.m で作成

# =============================================================================
# 設定 (MATLAB の exp_denoising_gstd_Base.m に合わせている)
# =============================================================================
NOISE_CONDITIONS = [
    # g     ps     pt     tint  pd
    (0.1,  0,     0,     0,    0),     # 1: g0.1
    (0.05, 0.05,  0,     0,    0),     # 2: g0.05 ps0.05
    (0.1,  0.05,  0,     0,    0),     # 3: g0.1 ps0.05
    (0.05, 0,     0.05,  0.5,  0),     # 4: g0.05 pt0.05
    (0.1,  0,     0.05,  0.5,  0),     # 5: g0.1 pt0.05
    (0.05, 0.05,  0.05,  0.5,  0),     # 6: g0.05 ps0.05 pt0.05
    (0.1,  0.05,  0.05,  0.5,  0),     # 7: g0.1 ps0.05 pt0.05
    (0.1,  0,     0,     0,    0.01),  # 8: g0.1 pd0.01
    (0.1,  0.05,  0.05,  0.5,  0.01),  # 9: g0.1 ps0.05 pt0.05 pd0.01
]

DEFAULTS = dict(
    images=["JasperRidge", "PaviaU"],
    noise_idx=[1],              # NOISE_CONDITIONS の番号 (1 始まり, MATLAB と同じ)
    graph=["noisy"],            # "noisy": 観測画像からグラフを作る / "oracle": GT 画像から作る (OG)
    rho=[0.98],                 # rho_radius (epsilon, alpha, beta の倍率)
    maxiter=20000,
    stopcri_idx=5,              # stopcri = 10^-stopcri_idx
)

# グラフの種類ごとのパラメータ (各値はリストで, 全組み合わせを実行する)
METHOD_PRESETS = {
    "noisy": dict(
        name="GASSTV_CondatVu",
        lambda_rho_sp=[0.6, 0.7, 0.8, 0.9, 1],
        lambda2=[0.1, 1, 10, 50],          # omega (GLR の重み)
        sigma_sp=["90", "med"],            # 空間グラフの sigma: "med" | "90" | 数値
        sigma_l=["med"],                   # スペクトルグラフの sigma: "med" | "90" | 数値
        num_segments=[4],
        order_filt=[5],                    # 代表スペクトルのメディアンフィルタ次数 (1 でフィルタなし)
        k_lap=[1000],                      # 相互 kNN の k (K-1 以上なら全結合)
    ),
    "oracle": dict(
        name="GASSTV_CondatVu_OG",
        lambda_rho_sp=[0.9, 1, 1.1],
        lambda2=[0.1, 1, 10],
        sigma_sp=["90", "med"],
        sigma_l=["med"],
        num_segments=[4],
        order_filt=[1],
        k_lap=[10],                        # MATLAB の OG 版は k_lap = 10 で固定
    ),
}
GRID_KEYS = ["lambda_rho_sp", "lambda2", "sigma_sp", "sigma_l", "num_segments", "order_filt", "k_lap"]


# =============================================================================
def parse_value(s: str):
    """"med" / "90" は文字列 (sigma のモード), それ以外は数値として解釈する."""
    if s.lower() in ("med", "90"):
        return s.lower()
    v = float(s)
    return int(v) if v.is_integer() and "." not in s and "e" not in s.lower() else v


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_argument_group("実験条件")
    g.add_argument("--images", nargs="+", default=DEFAULTS["images"],
                   help="JasperRidge, JasperRidge64, PaviaU, PaviaU120, PaviaU64, WashingtonDC, ...")
    g.add_argument("--noise_idx", nargs="+", type=int, default=DEFAULTS["noise_idx"],
                   help=f"ノイズ条件の番号 (1..{len(NOISE_CONDITIONS)})")
    g.add_argument("--noise", nargs=5, type=float, action="append", metavar=("G", "PS", "PT", "TINT", "PD"),
                   help="ノイズ条件を直接指定 (複数回指定可). 指定時は --noise_idx を無視")
    g.add_argument("--graph", nargs="+", choices=["noisy", "oracle"], default=DEFAULTS["graph"],
                   help="noisy: 観測からグラフ構築 / oracle: GT からグラフ構築 (Oracle Guide)")
    g.add_argument("--noise_model", choices=["auto", "g", "gs", "gt", "gst"], default="auto",
                   help="モデル (S, T の有無). auto はノイズ条件から MATLAB と同じ規則で選ぶ")

    g = p.add_argument_group("パラメータ (指定するとグラフの種類によらずこの値を使う. 複数値で全組み合わせ)")
    for k in GRID_KEYS:
        g.add_argument(f"--{k}", nargs="+", type=parse_value)
    g.add_argument("--rho", nargs="+", type=float, default=DEFAULTS["rho"], help="rho_radius")
    g.add_argument("--maxiter", type=int, default=DEFAULTS["maxiter"])
    g.add_argument("--stopcri_idx", type=int, default=DEFAULTS["stopcri_idx"], help="stopcri = 10^-idx")

    g = p.add_argument_group("アルゴリズムの細部")
    g.add_argument("--lambda_sp_ref", choices=["clean", "noisy", "guide"], default="clean",
                   help="空間グラフ制約の半径 lambda_sp を計算する画像. clean が MATLAB と同じ "
                        "(非 OG 版でも GT を使っている). guide はグラフと同じ画像")
    g.add_argument("--seed", type=int, default=0, help="ノイズ生成とセグメンテーションの乱数シード")

    g = p.add_argument_group("入出力・実行環境")
    g.add_argument("--data_dir", default=DEFAULT_DATA_DIR)
    g.add_argument("--obsv", choices=["auto", "dataset", "python"], default="auto",
                   help="観測の取得方法. dataset: dataset/<image>/<ノイズ条件>.mat を読む "
                        "(MATLAB の Matlab/Generate_dataset_denoising.m で作成, MATLAB 実験と同一のノイズ). "
                        "python: Python でノイズを生成. auto: dataset にあればそれを, なければ python")
    g.add_argument("--dataset_dir", default=str(DATASET_DIR))
    g.add_argument("--result_dir", default=str(ROOT / "result"))
    g.add_argument("--method_suffix", default="", help="手法名 (保存フォルダ名) の末尾に付ける文字列")
    g.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    g.add_argument("--track_every", type=int, default=1,
                   help="反復中の MPSNR/MSSIM を記録する間隔 (0 で記録しない. 大きくすると速い)")
    g.add_argument("--no_comp", action="store_true", help="GASSTV_comp_result への保存をしない")
    g.add_argument("--skip_existing", action="store_true", help="結果ファイルが既にあればスキップ")
    g.add_argument("--quiet", action="store_true", help="反復中の表示をしない")
    g.add_argument("--dry_run", action="store_true", help="実行する組み合わせを表示して終了")
    return p


def param_grid(lists: dict) -> list[dict]:
    """ParamsList2Comb 相当 (MATLAB の ndgrid と同じく先頭のパラメータが最も速く変わる順)."""
    keys = list(lists)
    combos = itertools.product(*[lists[k] for k in reversed(keys)])
    return [dict(zip(keys, reversed(c))) for c in combos]


def params_savetext(params: dict, graph: str, stopcri_idx: int, args) -> str:
    """MATLAB の get_params_savetext と同じファイル名."""
    sig = lambda v: v if isinstance(v, str) else f"{v:g}"  # noqa: E731
    s = (f"l{params['lambda_rho_sp']:.2g}_{params['lambda2']:.2g}"
         f"_sig{sig(params['sigma_sp'])}_{sig(params['sigma_l'])}"
         f"_ns{int(params['num_segments']):d}_fil{int(params['order_filt']):d}")
    if graph == "noisy" or params["k_lap"] != 10:
        s += f"_k{params['k_lap']:g}"
    s += f"_r{params['rho_radius']:.2f}_stop1e-{stopcri_idx:d}"
    # MATLAB と異なる設定のときは区別できるように付記する
    if args.lambda_sp_ref != "clean":
        s += f"_lsp{args.lambda_sp_ref}"
    if args.maxiter != DEFAULTS["maxiter"]:
        s += f"_it{args.maxiter}"
    return s


def get_observation(image: str, deg: dict, args) -> tuple[np.ndarray, np.ndarray, dict, dict]:
    path = Path(args.dataset_dir) / image / f"{obsv_dirname(deg)}.mat"
    use_dataset = args.obsv == "dataset" or (args.obsv == "auto" and path.exists())
    if use_dataset:
        if not path.exists():
            sys.exit(f"[error] {path} がありません. MATLAB で Matlab/Generate_dataset_denoising.m を実行してください.")
        HSI_clean, HSI_noisy, hsi, deg_ds = load_dataset(path)
        deg = {**deg_ds, **deg, "obsv_source": f"dataset:{image}/{path.name}"}
    else:
        HSI_clean, hsi = load_hsi(image, args.data_dir)
        HSI_noisy, deg = generate_obsv(HSI_clean, deg, seed=args.seed)
        deg["obsv_source"] = f"python:seed{args.seed}"
    return HSI_clean.astype(np.float32), HSI_noisy.astype(np.float32), hsi, deg


def main() -> None:
    args = build_parser().parse_args()
    stopcri = 10.0 ** -args.stopcri_idx
    disp_iters = tuple(sorted(set(range(1, 11)) | set(range(1000, args.maxiter + 1, 1000))))

    conditions = [tuple(c) for c in args.noise] if args.noise else [NOISE_CONDITIONS[i - 1] for i in args.noise_idx]

    # 手法ごとのパラメータの組み合わせ
    runs = []
    for graph in args.graph:
        preset = METHOD_PRESETS[graph]
        lists = {k: (getattr(args, k) if getattr(args, k) is not None else preset[k]) for k in GRID_KEYS}
        lists.update(maxiter=[args.maxiter], stopcri=[stopcri], rho_radius=args.rho)
        runs.append((graph, preset["name"] + args.method_suffix, param_grid(lists)))

    total = len(conditions) * len(args.images) * sum(len(r[2]) for r in runs)
    print(f"******* initium *******  (device: {args.device}, total runs: {total})")
    if args.dry_run:
        for graph, name, combos in runs:
            for prm in combos:
                print(f"  {name}: {params_savetext(prm, graph, args.stopcri_idx, args)}")
        print(f"  x noise conditions {conditions}\n  x images {args.images}")
        return

    result_dir = Path(args.result_dir)
    count = 0
    for cond in conditions:
        deg0 = dict(zip(["gaussian_sigma", "sparse_rate", "stripe_rate", "stripe_intensity", "deadline_rate"],
                        map(float, cond)))
        for image in args.images:
            HSI_clean, HSI_noisy, hsi, deg = get_observation(image, deg0, args)
            noise_model = select_noise_model(deg) if args.noise_model == "auto" else args.noise_model

            for graph, name, combos in runs:
                for idx_comb, prm in enumerate(combos, 1):
                    count += 1
                    params = dict(prm)
                    # epsilon, alpha, beta (exp_denoising_gstd_Base.m と同じ式)
                    ps, pd = deg["sparse_rate"], deg["deadline_rate"]
                    rate_except_sd = 1 - ps - 2 * pd + 2 * ps * pd
                    params["epsilon"] = params["rho_radius"] * deg["gaussian_sigma"] * np.sqrt(hsi["N"] * rate_except_sd)
                    params["alpha"] = params["rho_radius"] * (0.5 * hsi["N"] * (ps - 2 * ps * pd))
                    params["beta"] = (params["rho_radius"] * hsi["N"] * rate_except_sd
                                      * deg["stripe_rate"] * deg["stripe_intensity"] / 2)

                    savetext = params_savetext(params, graph, args.stopcri_idx, args)
                    sub = Path(f"denoising_{image}") / obsv_dirname(deg) / name / f"{savetext}.mat"
                    path_full = result_dir / "GASSTV_result" / sub
                    path_comp = result_dir / "GASSTV_comp_result" / sub

                    print("\n~~~ SETTINGS ~~~")
                    print(f"Method: {name} (graph: {graph}, model: {noise_model})")
                    print(f"Image: {image} Size: ({hsi['n1']}, {hsi['n2']}, {hsi['n3']}), obsv: {deg['obsv_source']}")
                    print("Noise: g {gaussian_sigma:g}, ps {sparse_rate:g}, pt {stripe_rate:g}, "
                          "tint {stripe_intensity:g}, pd {deadline_rate:g}".format(**deg))
                    print(f"Parameter settings: {savetext}")
                    print(f"Runs: ({count}/{total}), Params: ({idx_comb}/{len(combos)})")

                    if args.skip_existing and path_full.exists():
                        print("-> skip (exists)")
                        continue

                    opt = SolverOptions(
                        oracle=(graph == "oracle"), lambda_sp_ref=args.lambda_sp_ref,
                        device=args.device, seed=args.seed, track_every=args.track_every,
                        disp_iters=disp_iters, verbose=not args.quiet,
                    )
                    HSI_restored, removed_noise, other_result = gasstv_condatvu(
                        HSI_clean, HSI_noisy, params, noise_model, opt)

                    val_mpsnr = calc_mpsnr(HSI_restored, HSI_clean)
                    val_mssim = calc_mssim(HSI_restored, HSI_clean)
                    val_sam = calc_sam(HSI_restored, HSI_clean)
                    vals_psnr_per_band = psnr_per_band(HSI_restored, HSI_clean).numpy()
                    vals_ssim_per_band = ssim_per_band(HSI_restored, HSI_clean).numpy()
                    print("~~~ RESULTS ~~~")
                    print(f"MPSNR: {val_mpsnr:#.4g}\nMSSIM: {val_mssim:#.4g}\nSAM  : {val_sam:#.4g}")

                    save_mat(path_full, dict(
                        HSI_clean=HSI_clean, HSI_noisy=HSI_noisy, hsi=hsi, deg=deg, image=image,
                        HSI_restored=HSI_restored, removed_noise=removed_noise,
                        val_mpsnr=val_mpsnr, val_mssim=val_mssim, val_sam=val_sam,
                        vals_psnr_per_band=vals_psnr_per_band, vals_ssim_per_band=vals_ssim_per_band,
                        params=params, other_result=other_result,
                    ))
                    if not args.no_comp:
                        save_mat(path_comp, dict(
                            HSI_restored=HSI_restored, params=params,
                            val_mpsnr=val_mpsnr, val_mssim=val_mssim, val_sam=val_sam,
                        ))
                    print(f"Saved: {path_full}")

    print("******* finis *******")


if __name__ == "__main__":
    main()
