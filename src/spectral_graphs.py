"""スペクトルグラフ構築の候補手法 (検証用).

ノードはいずれも隣接バンド差分 (K = n3 - 1 個). すべての関数は共通の SpectralGraph を返す.

    baseline  : 現行法 (src.graphs.spectral_graph_laplacian). セグメントごと, 差分値の類似度, 符号なし
    abs       : 案 1. 代表スペクトル差分の絶対値 |delta| の類似度 + 符号付きラプラシアン (画像全体で 1 つ)
    window    : 案 2. |delta| の窓 (半幅 h) の類似度 + 符号付きラプラシアン. h = 0 で abs と一致
    coherence : 案 3. 差分画像どうしの低空間周波数コヒーレンス + 符号付きラプラシアン
    anchored  : 提案法. セグメントの代表スペクトルにつないだスペクトルグラフ (自己ループ付きの一般化ラプラシアン)
                R(u) = sum_p (D u_p - delta_s)^T (L_s + Q_s) (D u_p - delta_s)
                  L_s : 代表スペクトル差分 delta_s の値が近いノードどうしを結ぶ辺
                  Q_s : 各ノードを値 delta_s,i に固定した代表ノードへつなぐ自己ループ (セグメント内のばらつきが小さいほど強い)
                SpectralGraph.target = delta_s (ソルバーは D u_p - target に L をかける)

符号付きラプラシアン:  L = D - W o S,  D = diag(W 1),  S_ij = sigma_i sigma_j
    x^T L x = sum_{i<j} w_ij (x_i - S_ij x_j)^2  (半正定値. S が sigma sigma^T の形なら L = Sigma L_unsigned Sigma)
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import median_filter

from .graphs import spectral_graph_laplacian


# ---------------------------------------------------------------------------
# 共通の入れ物
# ---------------------------------------------------------------------------
@dataclass
class SpectralGraph:
    name: str
    L: np.ndarray                 # (S, K, K) ラプラシアン (符号付きの場合を含む)
    W: np.ndarray                 # (S, K, K) 重み (非負)
    sign: np.ndarray              # (S, K, K) エッジの符号 (+1 / -1). 符号なしグラフは全て +1
    labels: np.ndarray            # (n1, n2) 各画素のグラフ番号 (0..S-1). 画像全体で 1 つなら全て 0
    lam_max: np.ndarray           # (S,) 各 L の最大固有値
    info: dict = field(default_factory=dict)
    target: np.ndarray | None = None  # (S, K). あればソルバーは (D u_p - target[s]) に L_s をかける

    @property
    def S(self) -> int:
        return self.L.shape[0]

    @property
    def K(self) -> int:
        return self.L.shape[1]

    def signed_adjacency(self) -> np.ndarray:
        """画素数で重み付けして平均した符号付き隣接行列 (K, K). グラフ間の比較用."""
        counts = np.bincount(self.labels.ravel(), minlength=self.S).astype(float)
        a = counts / counts.sum()
        return np.einsum("s,sij->ij", a, self.W * self.sign)

    def normalized_laplacian_mean(self) -> np.ndarray:
        """L / lam_max を画素数で重み付けして平均したもの (K, K)."""
        counts = np.bincount(self.labels.ravel(), minlength=self.S).astype(float)
        a = counts / counts.sum()
        lam = np.where(self.lam_max > 0, self.lam_max, 1.0)
        return np.einsum("s,sij->ij", a, self.L / lam[:, None, None])


def _finalize(name, W, sign, labels, info) -> SpectralGraph:
    """W (S,K,K), sign (S,K,K) からラプラシアンと最大固有値を作る."""
    W = np.asarray(W, dtype=np.float64)
    sign = np.asarray(sign, dtype=np.float64)
    if W.ndim == 2:
        W, sign = W[None], sign[None]
    S, K, _ = W.shape
    L = np.zeros_like(W)
    lam = np.zeros(S)
    for s in range(S):
        Ws = (W[s] + W[s].T) / 2
        np.fill_diagonal(Ws, 0)
        Ls = np.diag(Ws.sum(axis=1)) - Ws * sign[s]
        Ls = (Ls + Ls.T) / 2
        L[s] = Ls
        lam[s] = max(float(np.linalg.eigvalsh(Ls).max()), 0.0)
        W[s] = Ws
    return SpectralGraph(name=name, L=L, W=W, sign=sign, labels=labels, lam_max=lam, info=info)


# ---------------------------------------------------------------------------
# 部品
# ---------------------------------------------------------------------------
def representative_spectrum(X: np.ndarray, rep: str = "colmed", order_filt: int = 5) -> np.ndarray:
    """画像全体の代表スペクトル (n3,).

    rep:
        "median" : 全画素の中央値 (現行法と同じ推定量)
        "colmed" : 列 (n2 方向の各列) ごとに中央値をとってから, 列方向に中央値.
                   ストライプ / デッドラインは列単位で乗るので, それらが外れ値の列として除かれる.
    order_filt: 波長方向のメディアンフィルタ次数 (端は nearest で延長. 現行の medfilt1 はゼロ埋め).
    """
    n3 = X.shape[2]
    X = np.asarray(X, dtype=np.float64)
    if rep == "median":
        r = np.median(X.reshape(-1, n3), axis=0)
    elif rep == "colmed":
        r = np.median(np.median(X, axis=0), axis=0)
    else:
        raise ValueError(f"unknown rep: {rep}")
    if order_filt > 1:
        r = median_filter(r, size=int(order_filt), mode="nearest")
    return r


def _sigma_from_d2(D2: np.ndarray, sigma) -> float:
    K = D2.shape[0]
    d2 = D2[np.triu_indices(K, 1)]
    d2 = d2[d2 > 0]
    if d2.size == 0:
        return 1.0
    if isinstance(sigma, str):
        if sigma == "med":
            return float(np.sqrt(np.median(d2)))
        if sigma == "90":
            return float(np.sqrt(np.percentile(d2, 90)))
        raise ValueError('sigma must be "med", "90" or numeric')
    return float(sigma)


def _sparsify(W: np.ndarray, k: int, knn: str = "mutual", w_min: float = 0.0,
              eps_c: float = 0.5) -> np.ndarray:
    """疎化.

    knn = "mutual" (kNN, min で対称化) | "union" (kNN, max. 現行法と同じ) | "none" (全結合)
        | "eps" (eps 近傍: 距離 d_ij < eps_c * sigma のペアだけ結ぶ. 次数はデータで決まり, 似た相手のいないノードは孤立する)
    W はガウス重み exp(-d^2 / 2 sigma^2) を想定している. このとき d < eps_c * sigma  <=>  w > exp(-eps_c^2 / 2).
    sigma を距離の中央値から決めている場合, eps_c = 1 で約半数のペアが残るので, 疎にするには eps_c < 1 にする.
    knn = "eps" で k > 0 のときは, 次数の上限として mutual kNN も併用する (k <= 0 なら上限なし).
    """
    K = W.shape[0]
    W = W.copy()
    np.fill_diagonal(W, 0)
    if knn == "eps":
        W[W <= np.exp(-eps_c**2 / 2)] = 0
        if 0 < k < K - 1:
            idx = np.argsort(-W, axis=1, kind="stable")[:, :k]
            Wk = np.zeros_like(W)
            rows = np.arange(K)[:, None]
            Wk[rows, idx] = W[rows, idx]
            W = np.minimum(Wk, Wk.T)
    elif knn != "none" and 0 < k < K - 1:
        idx = np.argsort(-W, axis=1, kind="stable")[:, :k]
        Wk = np.zeros_like(W)
        rows = np.arange(K)[:, None]
        Wk[rows, idx] = W[rows, idx]
        W = np.minimum(Wk, Wk.T) if knn == "mutual" else np.maximum(Wk, Wk.T)
    if w_min > 0:
        W[W < w_min] = 0
    return W


def _node_signs(delta: np.ndarray, tau: float) -> np.ndarray:
    """|delta| >= tau のノードに符号 +-1, それ以外 (符号が信頼できない) は 0."""
    s = np.sign(delta)
    s[np.abs(delta) < tau] = 0
    return s


def _node_signs_mode(delta: np.ndarray, mode: str, tau: float) -> np.ndarray:
    """ノードの符号.

    mode = "threshold" : |delta| >= tau のノードだけ +-1, それ以外 0 (0 を含む辺は S_ij = +1). 従来の既定
           "exact"     : sigma_i = sgn(delta_i), sgn(0) = +1. S = sigma sigma^T となり L = Sigma L_unsigned Sigma
                         (平衡な符号付きグラフ. 固有値は符号なしと同じ)
           "none"      : 符号なし (全て +1)
    """
    if mode is None:  # コマンドラインの sign_mode=none は None に変換されるため
        mode = "none"
    if mode == "threshold":
        return _node_signs(delta, tau)
    if mode == "exact":
        s = np.sign(delta)
        s[s == 0] = 1
        return s
    if mode == "none":
        return np.ones_like(delta)
    raise ValueError(f"unknown sign_mode: {mode}")


def _edge_signs_from_nodes(s: np.ndarray) -> np.ndarray:
    """S_ij = s_i s_j (どちらかが 0 なら +1)."""
    S = np.outer(s, s)
    S[S == 0] = 1
    return S


# ---------------------------------------------------------------------------
# 現行法
# ---------------------------------------------------------------------------
def graph_baseline(X: np.ndarray, num_segments: int = 4, sigma_l="med", k_lap: int = 1000,
                   order_filt: int = 5, seed: int = 0) -> SpectralGraph:
    """現行法 (src.graphs.spectral_graph_laplacian をそのまま使う)."""
    rng = np.random.default_rng(seed)
    L, lam, info = spectral_graph_laplacian(X, num_segments, sigma_l, k_lap, order_filt, rng)
    W = info["W_band_all"]
    g = SpectralGraph(name="baseline", L=L, W=W, sign=np.ones_like(W), labels=info["labels"],
                      lam_max=lam, info=dict(Rdiff=info["Rdiff"], val_sigma=info["val_sigma_l"]))
    return g


# ---------------------------------------------------------------------------
# 案 1 / 案 2: |delta| (の窓) の類似度 + 符号付きラプラシアン
# ---------------------------------------------------------------------------
def graph_abs_signed(X: np.ndarray, *, h: int = 0, rep: str = "colmed", order_filt: int = 5,
                     feature: str = "log", eps_q: float = 0.5, sigma="med",
                     k: int = 30, knn: str = "mutual", w_min: float = 0.0, eps_c: float = 0.5,
                     tau_q: float = 0.5, signed: bool = True, sign_mode: str = "threshold",
                     name: str | None = None) -> SpectralGraph:
    """|delta| の類似度でつなぐスペクトルグラフ (画像全体で 1 つ).

    feature:
        "abs" : f_i = |delta_i|
        "log" : f_i = log(|delta_i| + eps),  eps = |delta| の eps_q 分位点.
                変動の大きさを比で比べる (0.05 と 0.08 を「どちらも急峻」とみなす). eps 未満は 1 つの塊になる.
    h      : 窓の半幅. ノード i の特徴量を (f_{i-h}, ..., f_{i+h}) とする (端は鏡映). h = 0 で案 1.
    tau_q  : 符号を付けるノードの閾値 tau = |delta| の tau_q 分位点. それ未満のノードを含むエッジは S_ij = +1.
    signed : False なら通常の (符号なし) ラプラシアン (sign_mode = "none" と同じ).
    sign_mode : "threshold" (既定, tau_q 未満のノードの符号は未定) | "exact" (sgn(0) = +1 で全ノードに符号) | "none"
    knn, k, eps_c : 疎化 (_sparsify を参照). knn = "eps" で eps 近傍グラフ
    """
    if not signed or sign_mode is None:
        sign_mode = "none"
    r = representative_spectrum(X, rep, order_filt)
    delta = np.diff(r)
    K = delta.size
    a = np.abs(delta)
    if feature == "abs":
        f = a
        eps = 0.0
    elif feature == "log":
        eps = max(float(np.quantile(a, eps_q)), 1e-12)
        f = np.log(a + eps)
    else:
        raise ValueError(f"unknown feature: {feature}")

    if h > 0:
        fp = np.pad(f, h, mode="reflect")
        F = np.lib.stride_tricks.sliding_window_view(fp, 2 * h + 1)  # (K, 2h+1)
    else:
        F = f[:, None]
    G = F @ F.T
    n = np.diag(G)
    D2 = np.maximum(n[:, None] + n[None, :] - 2 * G, 0) / F.shape[1]
    np.fill_diagonal(D2, 0)

    sig = _sigma_from_d2(D2, sigma)
    W = np.exp(-D2 / (2 * sig**2))
    W = _sparsify(W, k, knn, w_min, eps_c)

    tau = float(np.quantile(a, tau_q))
    s_node = _node_signs_mode(delta, sign_mode, tau)
    S = _edge_signs_from_nodes(s_node)

    labels = np.zeros(X.shape[:2], dtype=np.int64)
    if name is None:
        name = ("window" if h > 0 else "abs") + ("" if sign_mode != "none" else "_unsigned")
    deg = W.sum(axis=1)
    info = dict(delta=delta, feature=feature, eps=eps, sigma=sig, tau=tau, h=h, rep=rep,
                order_filt=order_filt, k=k, knn=knn, eps_c=eps_c, sign_mode=sign_mode, node_sign=s_node,
                num_edges=int((np.triu(W, 1) > 0).sum()), num_isolated=int((deg == 0).sum()))
    return _finalize(name, W, S, labels, info)


# ---------------------------------------------------------------------------
# 案 3: 低空間周波数コヒーレンス
# ---------------------------------------------------------------------------
def graph_coherence(X: np.ndarray, *, patch: int = 20, step: int | None = None, f_max: int = 4,
                    exclude_col_axis: bool = True, prefilter: int = 0, min_sep: int = 2,
                    k: int = 30, knn: str = "mutual", w_min: float = 0.0,
                    signed: bool = True, name: str = "coherence") -> SpectralGraph:
    """差分画像 A_i = X(:,:,i+1) - X(:,:,i) どうしの低空間周波数コヒーレンスでつなぐ.

    patch, step : Welch 平均に使うパッチの大きさと間隔 (既定は 50% 重なり)
    f_max       : 使う空間周波数の上限 (パッチ内の DFT 番号で |k1|, |k2| <= f_max)
    exclude_col_axis : k1 = 0 (縦方向に一定な成分) を除く. ストライプとデッドラインは列単位 (縦に一定) なのでここに集中する.
    prefilter   : > 0 なら各バンドに空間メディアンフィルタ (prefilter x prefilter) をかけてから差分をとる
    min_sep     : |i - j| < min_sep のペアはつながない. 隣接差分 (|i-j| = 1) はバンドを共有するため
                  ノイズが相関 (-1/2) し, 見かけのコヒーレンスが生じる. min_sep = 1 で除外なし
    重み  w_ij = mean_{omega in Omega} C_ij(omega),  符号  S_ij = sign(Re sum_omega S_ij(omega))
    """
    X = np.asarray(X, dtype=np.float64)
    if prefilter > 0:
        X = median_filter(X, size=(prefilter, prefilter, 1), mode="nearest")
    A = np.diff(X, axis=2)  # (n1, n2, K)
    n1, n2, K = A.shape
    step = step or patch // 2

    win = np.outer(np.hanning(patch), np.hanning(patch))[:, :, None]
    k1 = np.fft.fftfreq(patch) * patch
    K1, K2 = np.meshgrid(k1, k1, indexing="ij")
    mask = (np.abs(K1) <= f_max) & (np.abs(K2) <= f_max)
    if exclude_col_axis:
        mask &= K1 != 0
    mask &= ~((K1 == 0) & (K2 == 0))

    coeffs = []
    for i0 in range(0, n1 - patch + 1, step):
        for j0 in range(0, n2 - patch + 1, step):
            P = A[i0:i0 + patch, j0:j0 + patch, :]
            P = (P - P.mean(axis=(0, 1), keepdims=True)) * win
            F = np.fft.fft2(P, axes=(0, 1))
            coeffs.append(F[mask])  # (M, K)
    Fq = np.stack(coeffs)  # (Q, M, K)
    Q = Fq.shape[0]
    Sxy = np.einsum("qmi,qmj->mij", Fq, Fq.conj()) / Q  # (M, K, K)
    P_ = np.real(np.einsum("mii->mi", Sxy))
    C = np.abs(Sxy) ** 2 / (P_[:, :, None] * P_[:, None, :] + 1e-30)
    W = C.mean(axis=0)
    ii, jj = np.indices((K, K))
    W[np.abs(ii - jj) < min_sep] = 0
    W = _sparsify(W, k, knn, w_min)

    if signed:
        S = np.sign(np.real(Sxy.sum(axis=0)))
        S[S == 0] = 1
    else:
        S = np.ones((K, K))
    labels = np.zeros((n1, n2), dtype=np.int64)
    info = dict(num_patches=Q, num_freq=int(mask.sum()), patch=patch, f_max=f_max, prefilter=prefilter,
                k=k, knn=knn, coherence_mean=W)
    return _finalize(name, W, S, labels, info)


# ---------------------------------------------------------------------------
# 名前で呼ぶための表
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 提案法: 代表スペクトルにつないだスペクトルグラフ
# ---------------------------------------------------------------------------
def segment_image(X: np.ndarray, num_segments: int, seg: str = "intensity", seed: int = 0,
                  n_pca: int = 5) -> np.ndarray:
    """セグメンテーション. 戻り値は (n1, n2) の 0..S'-1 (空のクラスタは詰める).

    seg:
        "intensity" : バンド平均画像の 1 次元 k-means (現行法と同じ)
        "spectral"  : 3x3 空間メディアンをかけたスペクトルの主成分 (n_pca 個) で k-means
    """
    from .graphs import kmeans_1d
    n1, n2, n3 = X.shape
    rng = np.random.default_rng(seed)
    if seg == "intensity":
        lab = kmeans_1d(np.asarray(X, np.float32).mean(axis=2), int(num_segments), rng).reshape(n1, n2)
    elif seg == "spectral":
        from scipy.cluster.vq import kmeans2
        Xs = median_filter(np.asarray(X, np.float64), size=(3, 3, 1), mode="nearest").reshape(-1, n3)
        Z = Xs - Xs.mean(axis=0)
        _, _, Vt = np.linalg.svd(Z, full_matrices=False)
        F = Z @ Vt[:n_pca].T
        _, lab = kmeans2(F, int(num_segments), minit="++", seed=seed)
        lab = lab.reshape(n1, n2)
    else:
        raise ValueError(f"unknown seg: {seg}")
    _, lab = np.unique(lab, return_inverse=True)
    return lab.reshape(n1, n2).astype(np.int64)


def _segment_representative(Xs_pix: np.ndarray, cols: np.ndarray, rep: str, order_filt: int) -> np.ndarray:
    """セグメント内の画素 (N_s, n3) から代表スペクトル (n3,).

    rep: "median" (全画素の中央値) | "colmed" (列ごとの中央値をとってから列方向に中央値. ストライプ・デッドラインに強い)
    """
    if rep == "colmed":
        ucols = np.unique(cols)
        meds = [np.median(Xs_pix[cols == c], axis=0) for c in ucols if np.sum(cols == c) >= 3]
        r = np.median(np.stack(meds), axis=0) if len(meds) >= 3 else np.median(Xs_pix, axis=0)
    elif rep == "median":
        r = np.median(Xs_pix, axis=0)
    else:
        raise ValueError(f"unknown rep: {rep}")
    if order_filt > 1:
        r = median_filter(r, size=int(order_filt), mode="nearest")
    return r


def graph_anchored(X: np.ndarray, *, num_segments: int = 4, seg: str = "intensity", rep: str = "median",
                   order_filt: int = 5, edges: str = "delta", sigma="med", k: int = 0, knn: str = "mutual",
                   eps_c: float = 0.5,
                   anchor: str = "uniform", beta_q: float = 1.0, tau_q: float = 1.0, seed: int = 0,
                   name: str | None = None) -> SpectralGraph:
    """提案法のスペクトルグラフ (セグメントごと).

    edges  : "delta" (代表スペクトル差分の値が近いノードを結ぶ. 現行法と同じ作り方) | "none" (自己ループのみ)
    sigma  : 辺のガウス重みの幅 ("med" | "90" | 数値). セグメントごとに決める
    k      : 辺の kNN 疎化 (0 なら全結合).  knn: "mutual" | "union" | "eps" (eps 近傍, eps_c. k > 0 なら次数の上限)
    anchor : "uniform" (q_i = 1. 既定)
             | "dispersion" (q_i = 1 / (sigma_i^2 + tau), sigma_i はセグメント内の差分の MAD 推定.
               観測から推定するとノイズ分散 (~0.02) が物質のばらつき (1e-6~1e-3) より桁違いに大きく,
               GT のばらつきとほぼ無相関になる (JasperRidge/PaviaU, 条件 9). oracle guide 用)
             | "none" (自己ループなし = 代表スペクトルからのずれに辺だけをかける)
    tau_q  : tau = tau_q * median_i(sigma_i^2)
    beta_q : 自己ループの強さ. 辺のラプラシアンを最大固有値 1 に, q を最大値 1 にそろえたうえで L + beta_q Q とする
    """
    X = np.asarray(X, dtype=np.float64)
    n1, n2, n3 = X.shape
    K = n3 - 1
    labels = segment_image(X, num_segments, seg, seed)
    S = int(labels.max()) + 1
    pix = X.reshape(-1, n3)
    lab = labels.reshape(-1)
    cols = np.tile(np.arange(n2), n1)  # 各画素の列番号 (C 順の reshape)

    L = np.zeros((S, K, K))
    W_all = np.zeros((S, K, K))
    lam = np.zeros(S)
    target = np.zeros((S, K))
    q_all = np.zeros((S, K))
    for s in range(S):
        m = lab == s
        P = pix[m]
        r = _segment_representative(P, cols[m], rep, order_filt)
        delta = np.diff(r)
        target[s] = delta

        if edges == "delta":
            D2 = (delta[:, None] - delta[None, :]) ** 2
            sig = _sigma_from_d2(D2, sigma)
            W = np.exp(-D2 / (2 * max(sig, 1e-12) ** 2))
            W = _sparsify(W, k, knn, 0.0, eps_c) if (k > 0 or knn == "eps") else W
            np.fill_diagonal(W, 0)
            Le = np.diag(W.sum(axis=1)) - W
            le = float(np.linalg.eigvalsh(Le).max())
            Le = Le / le if le > 0 else Le
        elif edges == "none":
            W = np.zeros((K, K))
            Le = np.zeros((K, K))
        else:
            raise ValueError(f"unknown edges: {edges}")

        if anchor == "dispersion":
            dP = np.diff(P, axis=1)
            mad = 1.4826 * np.median(np.abs(dP - np.median(dP, axis=0)), axis=0)
            v = mad**2
            q = 1.0 / (v + tau_q * np.median(v) + 1e-30)
            q = q / q.max()
        elif anchor == "uniform":
            q = np.ones(K)
        elif anchor == "none":
            q = np.zeros(K)
        else:
            raise ValueError(f"unknown anchor: {anchor}")
        q_all[s] = q

        A = Le + beta_q * np.diag(q)
        A = (A + A.T) / 2
        L[s] = A
        W_all[s] = W
        lam[s] = max(float(np.linalg.eigvalsh(A).max()), 0.0)

    if name is None:
        name = f"anchored_{edges}_{anchor}"
    info = dict(num_segments=S, seg=seg, rep=rep, order_filt=order_filt, edges=edges, anchor=anchor,
                beta_q=beta_q, tau_q=tau_q, k=k, q=q_all, segment_sizes=np.bincount(lab, minlength=S))
    return SpectralGraph(name=name, L=L, W=W_all, sign=np.ones_like(W_all), labels=labels,
                         lam_max=lam, info=info, target=target)


BUILDERS = {
    "baseline": graph_baseline,
    "anchored": graph_anchored,
    "abs": graph_abs_signed,
    "window": lambda X, **kw: graph_abs_signed(X, **{"h": 2, **kw}),
    "coherence": graph_coherence,
}


def build_graph(method: str, X: np.ndarray, **kwargs) -> SpectralGraph:
    return BUILDERS[method](X, **kwargs)
