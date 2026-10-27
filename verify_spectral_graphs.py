"""スペクトルグラフの候補手法を, GT から作ったグラフと観測 (noisy) から作ったグラフで比較する.

使い方:
    uv run verify_spectral_graphs.py                                  # JasperRidge, PaviaU / 全ノイズ (条件 9)
    uv run verify_spectral_graphs.py --images PaviaU --noise g0.1_ps0.05_pt0.05_pd0.01

指標 (各手法について, GT から作ったグラフ G_gt と noisy から作ったグラフ G_noisy の両方で計算):
  [A] 一致度 (G_gt vs G_noisy)
      corr_W    : 符号付き隣接行列 (上三角) のピアソン相関
      jacc_top  : 重みの大きい上位 E = 5K 本のエッジ集合の Jaccard 係数
      dist_L    : ||L~_gt - L~_noisy||_F / ||L~_gt||_F   (L~ = L / lam_max)
  [B] 滑らかさ (GSP): 正規化レイリー商  R(z) = sum_p z_p^T L~ z_p / sum_p ||z_p||^2  in [0, 1]
      R_clean   : z_p = D_lambda u_p (GT のスペクトル差分). 小さいほど GT がそのグラフ上で滑らか
      R_gauss   : z_p = D_lambda n_p (ガウスノイズ成分). 大きいほどノイズを高周波として分離できる
      R_noise   : z_p = D_lambda (y_p - u_p) (混合ノイズ全体)
      sel       : R_gauss / R_clean (選択性)
  [C] 代理デノイジング: スペクトルグラフの事前分布だけで各画素を復元
      u_p = argmin ||u - y_p||^2 + mu (D u)^T L~ (D u)   を mu について最良にしたときの MPSNR
      psnr_g    : y = GT + ガウスノイズ成分 のみ (グラフの事前分布としての良さを見る)
      psnr_mix  : y = 観測 (混合ノイズ). L2 当てはめなのでスパース・ストライプは除けない. 参考値
  [D] 意図どおりか (GT の代表スペクトル差分 delta_gt の |.| 上位 10% を「急峻」ノードとする)
      w_ss / w_ff / w_sf : 全エッジ重みのうち 急峻-急峻 / 平坦-平坦 / 急峻-平坦 の割合
      sign_ok   : 急峻-急峻エッジのうち, エッジの符号が sign(delta_i delta_j) と合う重みの割合
      iso_steep : 次数がほぼ 0 の急峻ノードの割合 (つながず残す = 罰しない)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scipy.io

from src.spectral_graphs import SpectralGraph, build_graph

ROOT = Path(__file__).resolve().parent

# 手法名 -> (builder 名, パラメータ)
METHODS = {
    "baseline": ("baseline", dict(num_segments=4, sigma_l="med", k_lap=1000, order_filt=5)),   # 現行の既定値
    "baseline_k10": ("baseline", dict(num_segments=4, sigma_l="med", k_lap=10, order_filt=5)),
    "abs_h0": ("abs", dict(feature="log", h=0, k=30, knn="mutual")),                           # 案 1
    "window_h1": ("abs", dict(feature="log", h=1, k=30, knn="mutual")),                        # 案 2
    "window_h2": ("abs", dict(feature="log", h=2, k=30, knn="mutual")),                        # 案 2
    "window_h2_unsigned": ("abs", dict(feature="log", h=2, k=30, knn="mutual", signed=False)), # 符号の効果
    "coherence": ("coherence", dict(min_sep=2, k=30)),                                         # 案 3
}


# ---------------------------------------------------------------------------
def mpsnr(x, ref) -> float:
    mse = ((x - ref) ** 2).mean(axis=(0, 1))
    return float(np.mean(10 * np.log10(1.0 / mse)))


def diff_op(n3: int) -> np.ndarray:
    D = np.zeros((n3 - 1, n3))
    D[np.arange(n3 - 1), np.arange(n3 - 1)] = -1
    D[np.arange(n3 - 1), np.arange(1, n3)] = 1
    return D


def gt_delta(C: np.ndarray) -> np.ndarray:
    return np.diff(np.median(C.reshape(-1, C.shape[2]), axis=0))


# ---------------------------------------------------------------------------
# [A] 一致度
def consistency(g_gt: SpectralGraph, g_n: SpectralGraph) -> dict:
    A1, A2 = g_gt.signed_adjacency(), g_n.signed_adjacency()
    K = A1.shape[0]
    iu = np.triu_indices(K, 1)
    a1, a2 = A1[iu], A2[iu]
    corr = float(np.corrcoef(a1, a2)[0, 1]) if a1.std() > 0 and a2.std() > 0 else np.nan
    E = 5 * K
    t1 = set(np.argsort(-np.abs(a1))[:E][np.abs(a1)[np.argsort(-np.abs(a1))[:E]] > 0])
    t2 = set(np.argsort(-np.abs(a2))[:E][np.abs(a2)[np.argsort(-np.abs(a2))[:E]] > 0])
    jacc = len(t1 & t2) / max(1, len(t1 | t2))
    L1, L2 = g_gt.normalized_laplacian_mean(), g_n.normalized_laplacian_mean()
    dist = float(np.linalg.norm(L1 - L2) / np.linalg.norm(L1))
    return dict(corr_W=corr, jacc_top=jacc, dist_L=dist)


# [B] 滑らかさ
def rayleigh(g: SpectralGraph, Z: np.ndarray) -> float:
    """Z: (N, K) 画素ごとのスペクトル差分. 各画素はその画素のグラフ (labels) で評価."""
    lab = g.labels.reshape(-1)
    num = 0.0
    for s in range(g.S):
        m = lab == s
        if not m.any():
            continue
        Lt = g.L[s] / (g.lam_max[s] if g.lam_max[s] > 0 else 1.0)
        Zs = Z[m]
        num += float(np.sum((Zs @ Lt) * Zs))
    return num / float(np.sum(Z**2))


def smoothness(g: SpectralGraph, C, Y, Ng) -> dict:
    n3 = C.shape[2]
    Zc = np.diff(C, axis=2).reshape(-1, n3 - 1)
    Zg = np.diff(Ng, axis=2).reshape(-1, n3 - 1)
    Zn = np.diff(Y - C, axis=2).reshape(-1, n3 - 1)
    rc, rg, rn = rayleigh(g, Zc), rayleigh(g, Zg), rayleigh(g, Zn)
    return dict(R_clean=rc, R_gauss=rg, R_noise=rn, sel=rg / rc if rc > 0 else np.nan)


# [C] 代理デノイジング
MU_GRID = np.logspace(-2, 3, 26)


def proxy_denoise(g: SpectralGraph, Yobs: np.ndarray, C: np.ndarray) -> tuple[float, float]:
    n1, n2, n3 = C.shape
    D = diff_op(n3)
    lab = g.labels.reshape(-1)
    Yf = Yobs.reshape(-1, n3)
    eig = []
    for s in range(g.S):
        Lt = g.L[s] / (g.lam_max[s] if g.lam_max[s] > 0 else 1.0)
        lam, V = np.linalg.eigh(D.T @ Lt @ D)
        eig.append((np.maximum(lam, 0), V))
    best = (-np.inf, np.nan)
    for mu in MU_GRID:
        U = np.empty_like(Yf)
        for s, (lam, V) in enumerate(eig):
            m = lab == s
            if m.any():
                U[m] = ((Yf[m] @ V) / (1 + mu * lam)) @ V.T
        p = mpsnr(U.reshape(n1, n2, n3), C)
        if p > best[0]:
            best = (p, mu)
    return best


def tikhonov_reference(Yobs, C) -> tuple[float, float]:
    """参照: グラフなしの 1 次差分 Tikhonov (L = I)."""
    n3 = C.shape[2]
    K = n3 - 1
    g = SpectralGraph("tikhonov", L=np.eye(K)[None], W=np.zeros((1, K, K)), sign=np.ones((1, K, K)),
                      labels=np.zeros(C.shape[:2], dtype=np.int64), lam_max=np.ones(1))
    return proxy_denoise(g, Yobs, C)


# [D] 意図どおりか
def intent(g: SpectralGraph, delta_gt: np.ndarray, top: float = 0.1) -> dict:
    A = g.signed_adjacency()
    K = A.shape[0]
    steep = np.zeros(K, bool)
    steep[np.argsort(-np.abs(delta_gt))[: max(1, int(round(top * K)))]] = True
    W = np.abs(A)
    tot = W[np.triu_indices(K, 1)].sum()
    ss = W[np.ix_(steep, steep)][np.triu_indices(steep.sum(), 1)].sum()
    ff = W[np.ix_(~steep, ~steep)][np.triu_indices((~steep).sum(), 1)].sum()
    sf = W[np.ix_(steep, ~steep)].sum()
    Ass = A[np.ix_(steep, steep)]
    true_sign = np.sign(np.outer(delta_gt[steep], delta_gt[steep]))
    iu = np.triu_indices(steep.sum(), 1)
    wss = np.abs(Ass[iu])
    sign_ok = float((wss * (np.sign(Ass[iu]) == true_sign[iu])).sum() / wss.sum()) if wss.sum() > 0 else np.nan
    deg = W.sum(axis=1)
    iso = float((deg[steep] < 1e-3 * deg.max()).mean())
    return dict(w_ss=ss / tot, w_ff=ff / tot, w_sf=sf / tot, sign_ok=sign_ok, iso_steep=iso)


def sign_agreement(C: np.ndarray, delta_gt: np.ndarray, top: float = 0.1) -> dict:
    """急峻ノードで, 各画素の差分の符号が代表スペクトルの符号と一致する画素の割合."""
    n3 = C.shape[2]
    dC = np.diff(C, axis=2).reshape(-1, n3 - 1)
    steep = np.argsort(-np.abs(delta_gt))[: max(1, int(round(top * (n3 - 1))))]
    agree = (np.sign(dC[:, steep]) == np.sign(delta_gt[steep])).mean(axis=0)
    return dict(mean=float(agree.mean()), min=float(agree.min()), nodes=steep.tolist(), per_node=agree.tolist())


# ---------------------------------------------------------------------------
def plot_all(image, C, Y, graphs, delta_gt, out_dir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.spectral_graphs import representative_spectrum

    # 1) 代表スペクトル差分
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    K = delta_gt.size
    x = np.arange(1, K + 1)
    for ax, (rep, of) in zip(axes, [("median", 5), ("colmed", 5)]):
        dg = np.diff(representative_spectrum(C, rep, of))
        dn = np.diff(representative_spectrum(Y, rep, of))
        ax.plot(x, dg, color="#2a6fbb", lw=1.4, label="GT")
        ax.plot(x, dn, color="#d9622b", lw=1.0, label="noisy")
        ax.axhline(0, color="gray", lw=0.5)
        ax.set_title(f"{image}: delta of representative spectrum (rep={rep}, medfilt={of})   "
                     f"corr={np.corrcoef(dg, dn)[0, 1]:.3f}", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
    axes[-1].set_xlabel("difference node i (band i+1 - band i)")
    fig.tight_layout()
    fig.savefig(out_dir / "delta_representative.png", dpi=130)
    plt.close(fig)

    # 2) 符号付き隣接行列
    names = list(graphs)
    fig, axes = plt.subplots(2, len(names), figsize=(3.0 * len(names), 6.3))
    for j, nm in enumerate(names):
        for i, key in enumerate(["gt", "noisy"]):
            A = graphs[nm][key].signed_adjacency()
            v = np.abs(A).max() or 1.0
            ax = axes[i, j]
            ax.imshow(A, cmap="RdBu_r", vmin=-v, vmax=v, interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{nm}\n({key})", fontsize=9)
    fig.suptitle(f"{image}: signed adjacency (red +, blue -), node order = band index", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_dir / "adjacency.png", dpi=130)
    plt.close(fig)

    # 3) 急峻ノードの結合: delta_gt の上に, 各急峻ノードの最強エッジを線で描く
    steep = np.argsort(-np.abs(delta_gt))[: max(1, int(round(0.1 * K)))]
    fig, axes = plt.subplots(len(names), 2, figsize=(11, 1.9 * len(names)), sharex=True)
    for j, nm in enumerate(names):
        for i, key in enumerate(["gt", "noisy"]):
            ax = axes[j, i]
            A = graphs[nm][key].signed_adjacency()
            ax.plot(x, delta_gt, color="#444", lw=0.9)
            ax.scatter(x[steep], delta_gt[steep], s=10, color="#444", zorder=3)
            vmax = np.abs(A).max() or 1.0
            for a in steep:
                for b in np.argsort(-np.abs(A[a]))[:3]:
                    if np.abs(A[a, b]) < 1e-6 * vmax:
                        continue
                    col = "#c0392b" if A[a, b] > 0 else "#2471a3"
                    ax.plot([x[a], x[b]], [delta_gt[a], delta_gt[b]], color=col,
                            alpha=0.25 + 0.75 * abs(A[a, b]) / vmax, lw=0.8)
            ax.set_title(f"{nm} ({key}): top-3 edges of steep nodes", fontsize=8)
            ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "steep_edges.png", dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--images", nargs="+", default=["JasperRidge", "PaviaU"])
    p.add_argument("--noise", default="g0.1_ps0.05_pt0.05_pd0.01")
    p.add_argument("--methods", nargs="+", default=list(METHODS))
    p.add_argument("--dataset_dir", default=str(ROOT / "dataset"))
    p.add_argument("--out_dir", default=str(ROOT / "result" / "graph_verification"))
    p.add_argument("--no_plot", action="store_true")
    args = p.parse_args()

    all_rows = {}
    for image in args.images:
        d = scipy.io.loadmat(Path(args.dataset_dir) / image / f"{args.noise}.mat", simplify_cells=True)
        C = d["HSI_clean"].astype(np.float64)
        Y = d["HSI_noisy"].astype(np.float64)
        Ng = np.asarray(d["deg"]["gaussian_noise"], dtype=np.float64)
        Yg = C + Ng
        delta_gt = gt_delta(C)
        out_dir = Path(args.out_dir) / image / args.noise
        out_dir.mkdir(parents=True, exist_ok=True)

        rows = {}
        ref = dict(noisy_g=mpsnr(Yg, C), noisy_mix=mpsnr(Y, C),
                   tikhonov_g=tikhonov_reference(Yg, C), tikhonov_mix=tikhonov_reference(Y, C),
                   sign_agreement=sign_agreement(C, delta_gt))
        graphs = {}
        for nm in args.methods:
            builder, kw = METHODS[nm]
            g_gt = build_graph(builder, C, **kw)
            g_n = build_graph(builder, Y, **kw)
            graphs[nm] = dict(gt=g_gt, noisy=g_n)
            row = dict(consistency=consistency(g_gt, g_n))
            for key, g in (("gt", g_gt), ("noisy", g_n)):
                pg, mug = proxy_denoise(g, Yg, C)
                pm, mum = proxy_denoise(g, Y, C)
                row[key] = dict(**smoothness(g, C, Y, Ng), **intent(g, delta_gt),
                                psnr_g=pg, mu_g=mug, psnr_mix=pm, mu_mix=mum,
                                lam_max=float(g.lam_max.max()))
            rows[nm] = row
            print(f"[{image}] {nm:13s} corr_W={row['consistency']['corr_W']:.3f} "
                  f"jacc={row['consistency']['jacc_top']:.3f} | "
                  + " | ".join(f"{k}: Rc={row[k]['R_clean']:.3f} sel={row[k]['sel']:.2f} "
                               f"psnr_g={row[k]['psnr_g']:.2f} w_sf={row[k]['w_sf']:.2f} "
                               f"sign_ok={row[k]['sign_ok']:.2f}" for k in ("gt", "noisy")))
        print(f"[{image}] reference: noisy_g={ref['noisy_g']:.2f}, tikhonov_g={ref['tikhonov_g'][0]:.2f}, "
              f"noisy_mix={ref['noisy_mix']:.2f}, tikhonov_mix={ref['tikhonov_mix'][0]:.2f}, "
              f"sign agreement (steep) mean={ref['sign_agreement']['mean']:.3f}")
        all_rows[image] = dict(methods=rows, reference=ref)
        with open(out_dir / "metrics.json", "w") as f:
            json.dump(dict(methods=rows, reference=ref, method_params={k: METHODS[k] for k in args.methods}),
                      f, indent=1, default=float)
        if not args.no_plot:
            plot_all(image, C, Y, graphs, delta_gt, out_dir)
    return all_rows


if __name__ == "__main__":
    main()
