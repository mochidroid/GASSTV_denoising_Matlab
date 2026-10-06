"""GASSTV (Condat-Vu 型 P-PDS) の torch 実装.

Matlab/methods/GASSTV_CondatVu/func_GASSTV_{g,gs,gt,gst}_for_denoising_CondatVu.m (非 OG 版) と
func_GASSTV_OG_{...}_for_denoising_CondatVu.m (Oracle Guide 版) を 1 つの関数にまとめたもの.

    min_{U,S,T}  ||D Dl(U)||_1  +  (omega/2) sum_p (Dl u_p)^T L_{seg(p)} (Dl u_p)
    s.t.  ||Wsp .* Dsp(U)||_1 <= lambda_sp,  ||S||_1 <= alpha,  ||T||_1 <= beta,  Dv(T) = 0,
          ||U + S + T - V||_2 <= epsilon,  U in [0, 1]

noise_model により S (スパースノイズ) と T (ストライプノイズ) の有無が変わる:
    "g": U のみ, "gs": U+S, "gt": U+T, "gst": U+S+T
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import torch

from . import ops
from .graphs import spatial_graph_weight, spectral_graph_laplacian
from .metrics import calc_mpsnr, calc_mssim

GAMMA2 = {"g": 1.0, "gs": 1 / 2, "gt": 1 / 2, "gst": 1 / 3}


def select_noise_model(deg: dict) -> str:
    """select_func_GASSTV_for_denoising_CondatVu.m と同じ分岐."""
    if deg["sparse_rate"] == 0 and deg["stripe_rate"] == 0:
        return "g"
    if deg["stripe_rate"] == 0:
        return "gs"
    if deg["sparse_rate"] == 0:
        return "gt"
    return "gst"


@dataclass
class SolverOptions:
    oracle: bool = False             # True: グラフを HSI_clean から作る (Oracle Guide)
    lambda_sp_ref: str = "clean"     # 空間グラフ制約の半径 lambda_sp を計算する画像 ("clean" | "noisy" | "guide")
    device: str = "cuda"
    seed: int = 0                    # セグメンテーション (k-means) の乱数シード
    track_every: int = 1             # 反復中の MPSNR/MSSIM を何反復ごとに記録するか (0: 記録しない)
    disp_iters: tuple = ()           # 進捗を表示する反復番号
    verbose: bool = True


def _norm(x: torch.Tensor) -> torch.Tensor:
    return torch.linalg.vector_norm(x)


def gasstv_condatvu(HSI_clean: np.ndarray, HSI_noisy: np.ndarray, params: dict,
                    noise_model: str, opt: SolverOptions):
    """GASSTV を実行する.

    Parameters
    ----------
    HSI_clean, HSI_noisy : (n1, n2, n3) float32
    params : lambda_rho_sp, lambda2, sigma_sp, sigma_l, num_segments, order_filt, k_lap,
             maxiter, stopcri, epsilon, alpha, beta
    noise_model : "g" | "gs" | "gt" | "gst"

    Returns
    -------
    HSI_restored (numpy), removed_noise (dict), output (dict)
    """
    log = print if opt.verbose else (lambda *a, **k: None)
    tag = "_OG" if opt.oracle else ""
    log(f"** Running GASSTV{tag}_{noise_model} (CondatVu, torch) **")

    dev = torch.device(opt.device)
    dt = torch.float32
    clean = torch.from_numpy(np.ascontiguousarray(HSI_clean, dtype=np.float32)).to(dev)
    noisy = torch.from_numpy(np.ascontiguousarray(HSI_noisy, dtype=np.float32)).to(dev)
    n1, n2, n3 = noisy.shape
    use_S = noise_model in ("gs", "gst")
    use_T = noise_model in ("gt", "gst")

    epsilon = float(params["epsilon"])
    alpha = float(params.get("alpha", 0.0))
    beta = float(params.get("beta", 0.0))
    lambda_rho_sp = float(params["lambda_rho_sp"])
    lambda2 = float(params["lambda2"])
    maxiter = int(params["maxiter"])
    stopcri = float(params["stopcri"])

    # ---------------- グラフ構築 ----------------
    guide = HSI_clean if opt.oracle else HSI_noisy
    guide = np.asarray(guide, dtype=np.float32)

    # (A) 空間グラフ
    log("~ Creating spatial graph ~")
    Wsp_np, val_sigma_sp = spatial_graph_weight(guide, params["sigma_sp"])
    log(f"sigma sp: {val_sigma_sp:5f}")
    Wsp = torch.from_numpy(Wsp_np).to(dev).unsqueeze(2)  # (n1, n2, 1, 4) -> バンド方向にブロードキャスト

    ref = {"clean": clean, "noisy": noisy,
           "guide": clean if opt.oracle else noisy}[opt.lambda_sp_ref]
    lambda_sp = float((Wsp * ops.D4(ref)).abs().sum(dtype=torch.float64)) * lambda_rho_sp

    # (B) スペクトルグラフ (セグメントごとのラプラシアン)
    log("~ Creating spectral graph ~")
    rng = np.random.default_rng(opt.seed)
    L_np, lam_max_vec, info_l = spectral_graph_laplacian(
        guide, int(params["num_segments"]), params["sigma_l"], int(params["k_lap"]),
        int(params["order_filt"]), rng)
    log(f"sigma_l: {info_l['val_sigma_l']:5f}")
    L = torch.from_numpy(L_np).to(dev, dt)  # (S, K, K)
    seg_flat = torch.from_numpy(info_l["labels"].reshape(-1)).to(dev)
    seg_idx = [torch.nonzero(seg_flat == s).squeeze(1) for s in range(L.shape[0])]
    seg_idx = [(s, idx) for s, idx in enumerate(seg_idx) if idx.numel() > 0]

    lam_max_vec = np.where(np.isfinite(lam_max_vec), lam_max_vec, 0)
    max_eig_L = float(lam_max_vec.max())

    def grad_glr(U: torch.Tensor) -> torch.Tensor:
        """compute_GLR_gradient: omega * Dl^T( L_{seg(p)} Dl(u_p) )."""
        V = ops.Dl_glr(U).reshape(n1 * n2, n3 - 1)
        W = torch.zeros_like(V)
        for s, idx in seg_idx:
            W[idx] = V[idx] @ L[s]
        return lambda2 * ops.Dlt_glr(W.reshape(n1, n2, n3 - 1))

    # ---------------- ステップサイズ ----------------
    # ||Dl^T L Dl|| <= 4 ||L||
    lipschitz_glr = 4 * lambda2 * max_eig_L
    opnorm_Dvh, opnorm_Dl, opnorm_Wsp, opnorm_Dsp, opnorm_Dv = 8, 4, 1, 16, 4
    gamma1_U = 1 / (lipschitz_glr / 2 + opnorm_Dvh * opnorm_Dl + opnorm_Wsp * opnorm_Dsp + 1)
    gamma1_S = 1.0
    gamma1_T = 1 / (opnorm_Dv + 1)
    gamma2 = GAMMA2[noise_model]

    # ---------------- 変数の初期化 ----------------
    z3 = lambda: torch.zeros((n1, n2, n3), dtype=dt, device=dev)  # noqa: E731
    U = z3()
    S = z3() if use_S else None
    T = z3() if use_T else None
    Y1 = torch.zeros((n1, n2, n3, 2), dtype=dt, device=dev)
    Y2 = torch.zeros((n1, n2, n3, 4), dtype=dt, device=dev)
    Y3 = z3() if use_T else None
    Y4 = z3()

    rec = {k: np.full(maxiter, np.nan, dtype=np.float32) for k in
           ("converge_rate_U", "converge_rate_S", "converge_rate_T", "converge_rate_N",
            "move_mpsnr", "move_mssim", "running_time", "l2ball")}

    # ---------------- メインループ (P-PDS) ----------------
    log("~~~ P-PDS STARTS ~~~")
    use_cuda = dev.type == "cuda"
    i = 0
    for i in range(1, maxiter + 1):
        t0 = time.perf_counter()

        # primal: U
        U_next = ops.proj_box(U - gamma1_U * (grad_glr(U) + ops.Dlt(ops.Dt(Y1))
                                              + ops.D4t(Wsp * Y2) + Y4), 0, 1)
        U_res = 2 * U_next - U
        sum_res = U_res

        # primal: S
        if use_S:
            S_next = ops.proj_l1ball(S - gamma1_S * Y4, alpha)
            S_res = 2 * S_next - S
            sum_res = sum_res + S_res

        # primal: T
        if use_T:
            T_next = ops.proj_l1ball(T - gamma1_T * (ops.Dvt(Y3) + Y4), beta)
            T_res = 2 * T_next - T
            sum_res = sum_res + T_res

        # dual: Y1 (SSTV)
        Y1_tmp = Y1 + gamma2 * ops.D(ops.Dl(U_res))
        Y1_next = Y1_tmp - gamma2 * ops.prox_l1(Y1_tmp / gamma2, 1 / gamma2)

        # dual: Y2 (空間グラフ L1 ボール制約)
        Y2_tmp = Y2 + gamma2 * Wsp * ops.D4(U_res)
        Y2_next = Y2_tmp - gamma2 * ops.proj_l1ball(Y2_tmp / gamma2, lambda_sp)

        # dual: Y3 (ストライプの縦方向平坦性 Dv(T) = 0)
        if use_T:
            Y3_next = Y3 + gamma2 * ops.Dv(T_res)

        # dual: Y4 (データ忠実度 L2 ボール)
        Y4_tmp = Y4 + gamma2 * sum_res
        Y4_next = Y4_tmp - gamma2 * ops.proj_l2ball(Y4_tmp / gamma2, noisy, epsilon)

        # 誤差
        N = noisy - U
        N_next = noisy - U_next
        if use_S:
            N, N_next = N - S, N_next - S_next
        if use_T:
            N, N_next = N - T, N_next - T_next

        cr_U = _norm(U_next - U) / _norm(U)
        stats = [cr_U, _norm(N_next - N) / _norm(N), _norm(N)]
        if use_S:
            stats.append(_norm(S_next - S) / _norm(S))
        if use_T:
            stats.append(_norm(T_next - T) / _norm(T))

        # 更新
        U, Y1, Y2, Y4 = U_next, Y1_next, Y2_next, Y4_next
        if use_S:
            S = S_next
        if use_T:
            T, Y3 = T_next, Y3_next

        if use_cuda:
            torch.cuda.synchronize()
        rec["running_time"][i - 1] = time.perf_counter() - t0

        vals = torch.stack(stats).float().cpu().numpy()
        rec["converge_rate_U"][i - 1] = vals[0]
        rec["converge_rate_N"][i - 1] = vals[1]
        rec["l2ball"][i - 1] = vals[2]
        j = 3
        if use_S:
            rec["converge_rate_S"][i - 1] = vals[j]
            j += 1
        if use_T:
            rec["converge_rate_T"][i - 1] = vals[j]

        show = i in opt.disp_iters
        stop = i >= 2 and vals[0] < stopcri
        if opt.track_every and (i % opt.track_every == 0 or show or stop or i == 1):
            rec["move_mpsnr"][i - 1] = calc_mpsnr(U, clean)
            rec["move_mssim"][i - 1] = calc_mssim(U, clean)

        if stop:
            break
        if show:
            log(f"Iter: {i}, Error: {vals[0]:0.6f}, MPSNR: {rec['move_mpsnr'][i-1]:#.4g}, "
                f"MSSIM: {rec['move_mssim'][i-1]:#.4g}, Time: {np.nansum(rec['running_time']):0.2f}.")

    mpsnr_last, mssim_last = calc_mpsnr(U, clean), calc_mssim(U, clean)
    log(f"Iter: {i}, Error: {rec['converge_rate_U'][i-1]:0.6f}, MPSNR: {mpsnr_last:#.4g}, "
        f"MSSIM: {mssim_last:#.4g}, Time: {np.nansum(rec['running_time']):0.2f}.")
    log("~~~ P-PDS ENDS ~~~")

    # ---------------- 出力の整理 (MATLAB の removed_noise / output と同じ構成) ----------------
    HSI_restored = U.cpu().numpy()
    noisy_np = np.asarray(HSI_noisy, dtype=np.float32)
    removed_noise = {"all_noise": noisy_np - HSI_restored}
    gaussian = noisy_np - HSI_restored
    if use_S:
        removed_noise["sparse_noise"] = S.cpu().numpy()
        gaussian = gaussian - removed_noise["sparse_noise"]
    if use_T:
        removed_noise["stripe_noise"] = T.cpu().numpy()
        gaussian = gaussian - removed_noise["stripe_noise"]
    removed_noise["gaussian_noise"] = gaussian

    keys = ["converge_rate_U", "converge_rate_N", "move_mpsnr", "move_mssim", "running_time", "l2ball"]
    if use_S:
        keys.append("converge_rate_S")
    if use_T:
        keys.append("converge_rate_T")
    output = {"iter": i, **{k: rec[k][:i] for k in keys}}
    output["Wsp"] = np.broadcast_to(Wsp_np[:, :, None, :], (n1, n2, n3, 4)).copy()
    output["L_delta"] = np.transpose(L_np, (1, 2, 0)).astype(np.float32)  # MATLAB と同じ [K x K x S]
    # 以下は Python 版で追加した情報
    output.update(
        noise_model=noise_model, oracle=opt.oracle, lambda_sp=lambda_sp,
        lambda_sp_ref=opt.lambda_sp_ref,
        val_sigma_sp=val_sigma_sp, val_sigma_l=info_l["val_sigma_l"],
        seg_labels=info_l["labels"] + 1,  # MATLAB と同じ 1-based
        lam_max=lam_max_vec, gamma1_U=gamma1_U, gamma1_S=gamma1_S, gamma1_T=gamma1_T, gamma2=gamma2,
    )
    return HSI_restored, removed_noise, output
