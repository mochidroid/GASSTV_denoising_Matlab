"""既存法 (MATLAB) と GASSTV 系 (Python) の結果を rho ごとに並べて比べる.

読むもの:
  既存法  : result/GASSTV_comp_result/denoising_<image>/<noise>/<method>/*.mat  (Matlab/run_existing_methods.m の出力)
  GASSTV系: result/spectral_graph_comparison/<sweep>/summary.csv               (run_spectral_graph_comparison.py / sweep_all.py の出力)

各 (画像, ノイズ条件, rho, 手法) について, その手法の他のパラメータを振った中で MPSNR が最良のものを表にする.
rho を使わない手法 (LRTDTV, TPTV, FastHyMix) は rho によらず 1 行で出す.

例:
    uv run compare_methods.py
    uv run compare_methods.py --sweeps sweep_rho_select sweep_all --rhos 0.9 0.95
    uv run compare_methods.py --noise g0.1_ps0.05_pt0.05_pd0.01
    uv run compare_methods.py --existing_dir result/old_results --out result/method_comparison_old.md  (古い既存法の結果)
結果は result/method_comparison.md にも保存する.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent

EXISTING = ["SSTV", "HSSTV_L1", "HSSTV_L12", "l0l1HTV", "LRTDTV", "TPTV", "FastHyMix"]
GRAPH_LABEL = {
    "baseline": "GASSTV (現行グラフ)",
    "p1_abs": "GASSTV 案1",
    "p2_window": "GASSTV 案2 (窓)",
    "p3_coherence": "GASSTV 案3",
    "anchored": "GASSTV 提案 (代表ノード + 辺)",
    "anchor_only": "GASSTV 代表ノードのみ",
    "residual_edges": "GASSTV ずれに辺のみ",
}


# ---------------------------------------------------------------------------
# MATLAB .mat の読み込み (v7.3 = HDF5 と, それ以前の形式の両方)
# ---------------------------------------------------------------------------
def _h5_value(obj):
    import h5py
    if isinstance(obj, h5py.Group):
        return {k: _h5_value(obj[k]) for k in obj.keys() if not k.startswith("#")}
    arr = obj[()]
    cls = obj.attrs.get("MATLAB_class", b"")
    cls = cls.decode() if isinstance(cls, bytes) else str(cls)
    if cls == "char":
        return "".join(chr(c) for c in np.asarray(arr).ravel())
    if cls == "string":  # MATLAB string は v7.3 では中身を読めない (参照になる)
        return None
    if isinstance(arr, np.ndarray):
        if arr.size == 1:
            return arr.ravel()[0].item()
        return arr
    return arr


def load_matlab_result(path: Path) -> dict:
    try:
        import h5py
        with h5py.File(path, "r") as f:
            return {k: _h5_value(f[k]) for k in f.keys() if not k.startswith("#")}
    except OSError:
        import scipy.io
        d = scipy.io.loadmat(path, simplify_cells=True)
        return {k: v for k, v in d.items() if not k.startswith("__")}


def collect_existing(result_dir: Path, methods, noises=None) -> list[dict]:
    rows = []
    base = result_dir / "GASSTV_comp_result"
    for img_dir in sorted(base.glob("denoising_*")):
        image = img_dir.name.replace("denoising_", "")
        for noise_dir in sorted(p for p in img_dir.iterdir() if p.is_dir()):
            if noises and noise_dir.name not in noises:
                continue
            for method in methods:
                mdir = noise_dir / method
                if not mdir.is_dir():
                    continue
                for f in sorted(mdir.glob("*.mat")):
                    try:
                        d = load_matlab_result(f)
                    except Exception as e:  # noqa: BLE001
                        print(f"[warn] 読めませんでした: {f} ({e})")
                        continue
                    prm = d.get("params") or {}
                    rho = prm.get("rho_radius") if isinstance(prm, dict) else None
                    rows.append(dict(image=image, noise=noise_dir.name, method=method,
                                     rho=None if rho is None else round(float(rho), 4),
                                     mpsnr=float(d["val_mpsnr"]), mssim=float(d["val_mssim"]),
                                     sam=float(d["val_sam"]), setting=f.stem))
    return rows


def collect_gasstv(sweep_dirs: list[Path], noises=None) -> list[dict]:
    rows = []
    for sd in sweep_dirs:
        csv_path = sd / "summary.csv"
        if not csv_path.exists():
            print(f"[warn] {csv_path} がありません")
            continue
        for r in csv.DictReader(open(csv_path, newline="", encoding="utf-8")):
            if noises and r["noise"] not in noises:
                continue
            if r.get("guide", "noisy") != "noisy":
                continue
            om = float(r["omega"])
            method = "GASSTV (スペクトルグラフ項なし)" if om == 0 else GRAPH_LABEL.get(r["graph"], r["graph"])
            rows.append(dict(image=r["image"], noise=r["noise"], method=method,
                             rho=round(float(r["rho"]), 4), mpsnr=float(r["mpsnr"]), mssim=float(r["mssim"]),
                             sam=float(r["sam"]),
                             setting=f"{sd.name}: omega={om:g} lrs={r['lambda_rho_sp']} sig={r['sigma_sp']} "
                                     f"{r.get('graph_kw', '')}"))
    return rows


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--result_dir", default=str(ROOT / "result"))
    p.add_argument("--existing_dir", default=None,
                   help="既存法の結果を読むフォルダ (GASSTV_comp_result の親. 既定: --result_dir). "
                        "古い結果は result/old_results")
    p.add_argument("--sweeps", nargs="+", default=["sweep_rho_select"],
                   help="result/spectral_graph_comparison/ の下のフォルダ名 (複数可)")
    p.add_argument("--existing", nargs="+", default=EXISTING)
    p.add_argument("--rhos", nargs="+", type=float, default=None, help="表に出す rho (既定: 見つかったもの全部)")
    p.add_argument("--noise", nargs="+", default=None, help="ノイズ条件のフォルダ名で絞る")
    p.add_argument("--metric", choices=["mpsnr", "mssim", "sam"], default="mpsnr", help="最良を選ぶ指標")
    p.add_argument("--out", default=str(ROOT / "result" / "method_comparison.md"))
    a = p.parse_args()

    res = Path(a.result_dir)
    rows = collect_existing(Path(a.existing_dir) if a.existing_dir else res, a.existing, a.noise)
    rows += collect_gasstv([res / "spectral_graph_comparison" / s for s in a.sweeps], a.noise)
    if not rows:
        print("結果が見つかりません")
        return
    better = (lambda x, y: x < y) if a.metric == "sam" else (lambda x, y: x > y)

    # (image, noise, method, rho) ごとに最良
    best = {}
    for r in rows:
        k = (r["image"], r["noise"], r["method"], r["rho"])
        if k not in best or better(r[a.metric], best[k][a.metric]):
            best[k] = r
    rhos = sorted({r["rho"] for r in rows if r["rho"] is not None})
    if a.rhos:
        rhos = [x for x in rhos if any(abs(x - y) < 1e-6 for y in a.rhos)]

    order = a.existing + ["GASSTV (スペクトルグラフ項なし)"] + list(GRAPH_LABEL.values())
    lines = ["# 既存法と GASSTV 系の比較 (rho ごと)", "",
             f"各手法の他のパラメータを振った中で {a.metric.upper()} が最良の値. 太字は各列の最良.", ""]
    for image in sorted({k[0] for k in best}):
        for noise in sorted({k[1] for k in best if k[0] == image}):
            methods = sorted({k[2] for k in best if k[0] == image and k[1] == noise},
                             key=lambda m: order.index(m) if m in order else len(order))
            lines += [f"## {image} / {noise}", ""]
            head = "| 手法 | " + " | ".join(f"ρ={x:g} MPSNR | MSSIM | SAM" for x in rhos) + " | ρ なし MPSNR | MSSIM | SAM |"
            lines += [head, "|---|" + "---|---|---|" * (len(rhos) + 1)]
            # 列ごとの最良
            col_best = {}
            for x in rhos + [None]:
                vals = [best[(image, noise, m, x)] for m in methods if (image, noise, m, x) in best]
                if vals:
                    col_best[x] = dict(mpsnr=max(v["mpsnr"] for v in vals), mssim=max(v["mssim"] for v in vals),
                                       sam=min(v["sam"] for v in vals))
            for m in methods:
                cells = []
                for x in rhos + [None]:
                    r = best.get((image, noise, m, x))
                    if r is None:
                        cells += ["–", "–", "–"]
                        continue
                    cb = col_best[x]
                    f = lambda key, fmt: (f"**{r[key]:{fmt}}**" if abs(r[key] - cb[key]) < 1e-9 else f"{r[key]:{fmt}}")  # noqa: E731
                    cells += [f("mpsnr", ".3f"), f("mssim", ".4f"), f("sam", ".3f")]
                lines.append(f"| {m} | " + " | ".join(cells) + " |")
            lines.append("")
            # rho ごとの最良の設定
            lines.append("<details><summary>最良の設定</summary>\n")
            for m in methods:
                for x in rhos + [None]:
                    r = best.get((image, noise, m, x))
                    if r is not None:
                        lines.append(f"- {m}, ρ={x if x is not None else '–'}: {r['setting']}")
            lines += ["\n</details>", ""]
    text = "\n".join(lines)
    print(text)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(text, encoding="utf-8")
    print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
