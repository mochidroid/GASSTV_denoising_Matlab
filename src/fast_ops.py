"""CPU 実行用の高速化 (結果は src.ops と同一).

proj_l1ball: ソートの代わりに Michelot (Newton) 法で閾値を求める. 閾値は単調増加し,
有限回で sort 版と同じ厳密解に到達する. 閾値以下になった要素は二度と戻らないので,
反復ごとに候補集合を縮める.
"""
from __future__ import annotations

import torch

from . import ops


def proj_l1ball_fast(x: torch.Tensor, alpha: float) -> torch.Tensor:
    ax = x.abs()
    a = ax.reshape(-1)
    total = a.sum(dtype=torch.float64)
    if total <= alpha:
        return x
    th = (total - alpha) / a.numel()
    cand = a
    for _ in range(1000):
        cand = cand[cand > th]
        th_new = (cand.sum(dtype=torch.float64) - alpha) / cand.numel()
        if th_new <= th:
            break
        th = th_new
    th = float(max(th, 0.0))
    return torch.sign(x) * torch.clamp(ax - th, min=0)


def enable() -> None:
    ops.proj_l1ball = proj_l1ball_fast


_COMPILED = ("D", "Dt", "Dl", "Dlt", "Dv", "Dvt", "Dl_glr", "Dlt_glr", "D4", "D4t", "prox_l1", "proj_box")


def enable_compile() -> None:
    """差分作用素などの要素演算を torch.compile で融合する (CPU で数倍速くなる. 結果は同一)."""
    import torch
    for name in _COMPILED:
        setattr(ops, name, torch.compile(getattr(ops, name), dynamic=False))
