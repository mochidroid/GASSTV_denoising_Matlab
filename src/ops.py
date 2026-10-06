"""差分作用素と prox / 射影 (torch).

テンソルの軸は MATLAB と同じ (n1, n2, n3[, 方向]).
すべて Neumann 境界 (前進差分の最後の要素が 0).
"""

from __future__ import annotations

import torch


# ---------------------------------------------------------------------------
# 1 方向の前進差分とその随伴
# ---------------------------------------------------------------------------
def diff_fwd(z: torch.Tensor, dim: int) -> torch.Tensor:
    """z([2:end, end]) - z  (dim 方向)."""
    out = torch.zeros_like(z)
    n = z.shape[dim]
    out.narrow(dim, 0, n - 1).copy_(z.narrow(dim, 1, n - 1) - z.narrow(dim, 0, n - 1))
    return out


def diff_adj(z: torch.Tensor, dim: int) -> torch.Tensor:
    """diff_fwd の随伴: cat(-z(1), -z(2:end-1) + z(1:end-2), z(end-1))."""
    n = z.shape[dim]
    return torch.cat(
        [
            -z.narrow(dim, 0, 1),
            -z.narrow(dim, 1, n - 2) + z.narrow(dim, 0, n - 2),
            z.narrow(dim, n - 2, 1),
        ],
        dim=dim,
    )


# SSTV 用: 縦横 2 方向 (n1, n2, n3, 2)
def D(z):
    return torch.stack([diff_fwd(z, 0), diff_fwd(z, 1)], dim=-1)


def Dt(y):
    return diff_adj(y[..., 0], 0) + diff_adj(y[..., 1], 1)


def Dl(z):
    return diff_fwd(z, 2)


def Dlt(z):
    return diff_adj(z, 2)


def Dv(z):
    return diff_fwd(z, 0)


def Dvt(z):
    return diff_adj(z, 0)


# GLR 用のバンド差分 (パディングなし: n3 -> n3-1)
def Dl_glr(z):
    return z[:, :, 1:] - z[:, :, :-1]


def Dlt_glr(z):
    return torch.cat([-z[:, :, :1], -z[:, :, 1:] + z[:, :, :-1], z[:, :, -1:]], dim=2)


# 空間グラフ用: 4 方向 [下, 右, 右下, 右上] (D4_Neumann_GPU.m)
def D4(z):
    n1, n2 = z.shape[:2]
    out = torch.zeros((*z.shape, 4), dtype=z.dtype, device=z.device)
    out[..., 0] = diff_fwd(z, 0)
    out[..., 1] = diff_fwd(z, 1)
    out[: n1 - 1, : n2 - 1, :, 2] = z[1:, 1:] - z[:-1, :-1]
    out[1:, : n2 - 1, :, 3] = z[:-1, 1:] - z[1:, :-1]
    return out


def D4t(y):
    """D4t_Neumann_GPU.m."""
    out = diff_adj(y[..., 0], 0) + diff_adj(y[..., 1], 1)
    z_lt = y[..., 2].clone()
    z_lt[-1] = 0
    z_lt[:, -1] = 0
    out = out + torch.roll(z_lt, shifts=(1, 1), dims=(0, 1)) - z_lt
    z_lb = y[..., 3].clone()
    z_lb[0] = 0
    z_lb[:, -1] = 0
    out = out + torch.roll(z_lb, shifts=(-1, 1), dims=(0, 1)) - z_lb
    return out


# ---------------------------------------------------------------------------
# prox / 射影
# ---------------------------------------------------------------------------
def prox_l1(x: torch.Tensor, gamma) -> torch.Tensor:
    """ProxL1norm.m (ソフト閾値処理)."""
    return torch.sign(x) * torch.clamp(x.abs() - gamma, min=0)


def proj_box(x: torch.Tensor, a: float = 0.0, b: float = 1.0) -> torch.Tensor:
    return torch.clamp(x, a, b)


def proj_l1ball(x: torch.Tensor, alpha: float) -> torch.Tensor:
    """ProjFastL1Ball.m (半径 alpha の L1 ボールへの射影).

    ||x||_1 <= alpha のときは閾値が 0 になるので, ソートせずにそのまま返す (結果は同一).
    閾値計算の累積和は float64 で行う.
    """
    a = x.abs().reshape(-1)
    if a.sum(dtype=torch.float64) <= alpha:
        return x
    srt = torch.sort(a, descending=True).values.to(torch.float64)
    k = torch.arange(1, srt.numel() + 1, device=x.device, dtype=torch.float64)
    th = torch.clamp(((torch.cumsum(srt, 0) - alpha) / k).max(), min=0).to(x.dtype)
    return torch.sign(x) * torch.clamp(x.abs() - th, min=0)


def proj_l2ball(u: torch.Tensor, f: torch.Tensor, epsilon: float) -> torch.Tensor:
    """ProjL2ball.m (中心 f, 半径 epsilon の L2 ボールへの射影)."""
    t = u - f
    r = torch.linalg.vector_norm(t)
    if r > epsilon:
        return f + (epsilon / r) * t
    return u
