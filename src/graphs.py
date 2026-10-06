"""空間グラフ・スペクトルグラフの構築.

MATLAB 側の Create_SpatialGraphWeight.m / Create_SpectralGraphLaplacian.m に対応する.
グラフ構築は 1 回だけなので numpy (CPU, float64) で行う.
"""

from __future__ import annotations

import numpy as np

_EPS_SINGLE = float(np.finfo(np.float32).eps)


def prctile(x: np.ndarray, p: float) -> float:
    """MATLAB prctile と同じ補間 (numpy の 'hazen')."""
    return float(np.percentile(x, p, method="hazen"))


# ---------------------------------------------------------------------------
# 空間グラフ
# ---------------------------------------------------------------------------
def spatial_graph_weight(X: np.ndarray, sigma_sp) -> tuple[np.ndarray, float]:
    """Create_SpatialGraphWeight.m 相当.

    ガイド画像 (バンド平均) の 4 近傍 [下, 右, 右下, 右上] の差分からガウス重みを作る.
    戻り値 W は (n1, n2, 4). 全バンド共通なので MATLAB のように n3 方向へは複製しない.
    """
    g = X.astype(np.float32).mean(axis=2)
    n1, n2 = g.shape
    inf = np.float32(np.inf)

    grad = np.full((n1, n2, 4), inf, dtype=np.float32)
    grad[:-1, :, 0] = g[:-1, :] - g[1:, :]          # (i+1, j)
    grad[:, :-1, 1] = g[:, :-1] - g[:, 1:]          # (i, j+1)
    grad[:-1, :-1, 2] = g[:-1, :-1] - g[1:, 1:]     # (i+1, j+1) 右下
    grad[1:, :-1, 3] = g[1:, :-1] - g[:-1, 1:]      # (i-1, j+1) 右上

    finite = np.abs(grad[np.isfinite(grad)])
    if isinstance(sigma_sp, str):
        mode = sigma_sp.lower()
        if mode == "med":
            val = float(np.median(finite))
        elif mode == "90":
            g0 = prctile(finite, 90)  # 代表的な勾配の大きさ
            w0 = 0.2                  # g0 での重み
            val = max(g0 / np.sqrt(2 * np.log(1 / w0)), _EPS_SINGLE)
        else:
            raise ValueError('sigma_sp must be "med", "90", or a numeric value.')
    else:
        val = float(sigma_sp)

    W = np.exp(-(grad.astype(np.float64) ** 2) / (val**2) / 2).astype(np.float32)
    return W, val


# ---------------------------------------------------------------------------
# セグメンテーション (imsegkmeans の代替)
# ---------------------------------------------------------------------------
def kmeans_1d(x: np.ndarray, k: int, rng: np.random.Generator,
              num_attempts: int = 3, max_iter: int = 100) -> np.ndarray:
    """輝度値の k-means (k-means++ 初期化, 試行 num_attempts 回の最良).

    MATLAB の imsegkmeans(グレースケール画像, k) と同じ設定 (NumAttempts=3, MaxIterations=100)
    だが乱数系列が異なるため, ラベルが完全一致する保証はない.
    ラベルはクラスタ中心の輝度が小さい順に 0..k-1 に並べ替えて返す.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    n = x.size
    best_sse, best_labels, best_c = np.inf, None, None
    for _ in range(num_attempts):
        # k-means++
        c = [x[rng.integers(n)]]
        for _ in range(1, k):
            d2 = np.min((x[:, None] - np.asarray(c)[None, :]) ** 2, axis=1)
            p = d2 / d2.sum() if d2.sum() > 0 else np.full(n, 1 / n)
            c.append(x[rng.choice(n, p=p)])
        c = np.asarray(c)
        labels = None
        for _ in range(max_iter):
            new_labels = np.argmin((x[:, None] - c[None, :]) ** 2, axis=1)
            if labels is not None and np.array_equal(new_labels, labels):
                break
            labels = new_labels
            for j in range(k):
                m = labels == j
                if m.any():
                    c[j] = x[m].mean()
        sse = float(np.sum((x - c[labels]) ** 2))
        if sse < best_sse:
            best_sse, best_labels, best_c = sse, labels, c.copy()
    order = np.argsort(best_c)
    remap = np.empty(k, dtype=np.int64)
    remap[order] = np.arange(k)
    return remap[best_labels]


def medfilt1(x: np.ndarray, n: int) -> np.ndarray:
    """MATLAB medfilt1(x, n) (端はゼロ埋め) 相当."""
    n = int(n)
    if n <= 1:
        return x.copy()
    left = (n - 1) // 2 if n % 2 else n // 2
    right = n - 1 - left
    xp = np.concatenate([np.zeros(left), x, np.zeros(right)])
    win = np.lib.stride_tricks.sliding_window_view(xp, n)
    return np.median(win, axis=1)


# ---------------------------------------------------------------------------
# スペクトルグラフ (Segment-wise GLR)
# ---------------------------------------------------------------------------
def spectral_graph_laplacian(X: np.ndarray, num_segments: int, sigma_l, k_lap: int,
                             order_filt: int, rng: np.random.Generator,
                             labels: np.ndarray | None = None):
    """Create_SpectralGraphLaplacian.m 相当 (segment-wise spectral GLR graph).

    ノードは隣接バンド差分 (K = n3-1 個). セグメント s ごとに K x K のラプラシアン L_s を作る.

    Returns
    -------
    L : (S, K, K) float64
    lam_max : (S,) 各 L_s の最大固有値
    info : dict (labels は 0-based, (n1, n2))
    """
    n1, n2, n3 = X.shape
    K = n3 - 1
    S = int(num_segments)
    k_lap = max(0, min(int(k_lap), max(0, K - 1)))

    # 1) バンド平均画像をセグメンテーション
    guide = X.astype(np.float32).mean(axis=2)
    if labels is None:
        labels = kmeans_1d(guide, S, rng).reshape(n1, n2)

    # 2) 各セグメントの代表スペクトル (中央値 -> メディアンフィルタ)
    pix = X.reshape(-1, n3)
    lab = labels.reshape(-1)
    r = np.zeros((S, n3))
    counts = np.zeros(S, dtype=np.int64)
    for s in range(S):
        m = lab == s
        counts[s] = m.sum()
        if counts[s] > 0:
            r[s] = medfilt1(np.median(pix[m], axis=0).astype(np.float64), order_filt)

    # 3) 隣接バンド差分
    Rdiff = r[:, 1:] - r[:, :-1]  # (S, K)

    # 4) 全セグメント共通の距離 (sigma_l 決定用)
    GTG = Rdiff.T @ Rdiff
    nrm2 = np.diag(GTG)
    D2_global = nrm2[:, None] + nrm2[None, :] - 2 * GTG
    np.fill_diagonal(D2_global, 0)

    # 5) sigma_l
    d2 = D2_global[np.triu_indices(K, 1)]
    d2 = d2[d2 > 0]
    if d2.size == 0:
        d2 = np.array([1.0])
    if isinstance(sigma_l, str):
        mode = sigma_l.lower()
        if mode == "med":
            val_sigma_l = np.sqrt(np.median(d2) / S)
        elif mode == "90":
            val_sigma_l = np.sqrt(prctile(d2, 90))
        else:
            raise ValueError('sigma_l must be "med", "90", or a numeric value.')
    else:
        val_sigma_l = float(sigma_l)
    val_sigma_l = max(float(val_sigma_l), float(np.finfo(np.float64).eps))

    # 6) セグメントごとの GLR グラフ
    L = np.zeros((S, K, K))
    W_all = np.zeros((S, K, K))
    lam_max = np.zeros(S)
    for s in range(S):
        v = Rdiff[s]
        n2s = v**2
        D2_s = n2s[:, None] + n2s[None, :] - 2 * np.outer(v, v)
        np.fill_diagonal(D2_s, 0)
        W_s = np.exp(-D2_s / (2 * val_sigma_l**2))
        np.fill_diagonal(W_s, 0)

        # 相互 kNN による疎化 (対称化は max)
        if k_lap > 0 and K > 1:
            kk = min(k_lap, K - 1)
            Wk = np.zeros((K, K))
            idx = np.argsort(-W_s, axis=1, kind="stable")[:, :kk]
            rows = np.arange(K)[:, None]
            Wk[rows, idx] = W_s[rows, idx]
            W_s = np.maximum(Wk, Wk.T)

        L_s = np.diag(W_s.sum(axis=1)) - W_s
        L_s = (L_s + L_s.T) / 2
        ev = np.linalg.eigvalsh(L_s)
        ev = ev[np.isfinite(ev)]
        lam = float(ev.max()) if ev.size else 0.0
        lam_max[s] = max(lam, 0.0)
        L[s] = L_s
        W_all[s] = W_s

    info = dict(
        K=K, num_segments=S, order_filt=order_filt, segment_sizes=counts,
        labels=labels, representative=r, Rdiff=Rdiff, W_band_all=W_all,
        val_sigma_l=val_sigma_l, k_lap=k_lap,
    )
    return L, lam_max, info
