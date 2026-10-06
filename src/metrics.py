"""評価指標 (func_metrics/*.m 相当). 入力は (n1, n2, n3) の torch.Tensor または numpy 配列."""

from __future__ import annotations

import math

import numpy as np
import torch
import torch.nn.functional as F


def _t(x, device=None) -> torch.Tensor:
    if isinstance(x, np.ndarray):
        x = torch.from_numpy(x)
    return x.to(device=device, dtype=torch.float64)


def psnr_per_band(restored, clean) -> torch.Tensor:
    """MATLAB psnr(A(:,:,l), ref(:,:,l)) (ピーク値 1)."""
    x, y = _t(restored), _t(clean, restored.device if torch.is_tensor(restored) else None)
    mse = ((x - y) ** 2).mean(dim=(0, 1))
    return 10 * torch.log10(1.0 / mse)


def calc_mpsnr(restored, clean) -> float:
    return float(psnr_per_band(restored, clean).mean())


_GAUSS_CACHE: dict = {}


def _gauss_kernel(device, radius: float = 1.5) -> torch.Tensor:
    key = (str(device), radius)
    if key not in _GAUSS_CACHE:
        r = math.ceil(3 * radius)
        ax = torch.arange(-r, r + 1, dtype=torch.float64, device=device)
        g = torch.exp(-(ax**2) / (2 * radius**2))
        _GAUSS_CACHE[key] = g / g.sum()
    return _GAUSS_CACHE[key]


def _gfilt(x: torch.Tensor, g: torch.Tensor) -> torch.Tensor:
    """x: (B, 1, H, W) に分離型ガウスフィルタ (境界は replicate)."""
    r = (g.numel() - 1) // 2
    x = F.pad(x, (r, r, r, r), mode="replicate")
    x = F.conv2d(x, g.view(1, 1, -1, 1))
    return F.conv2d(x, g.view(1, 1, 1, -1))


def ssim_per_band(restored, clean, dynamic_range: float = 1.0) -> torch.Tensor:
    """MATLAB ssim (Radius=1.5, replicate パディング, 既定の正則化定数) をバンドごとに計算."""
    x = _t(restored)
    y = _t(clean, x.device)
    g = _gauss_kernel(x.device)
    x = x.permute(2, 0, 1).unsqueeze(1)
    y = y.permute(2, 0, 1).unsqueeze(1)
    C1 = (0.01 * dynamic_range) ** 2
    C2 = (0.03 * dynamic_range) ** 2
    mux, muy = _gfilt(x, g), _gfilt(y, g)
    mux2, muy2, muxy = mux * mux, muy * muy, mux * muy
    sx2 = _gfilt(x * x, g) - mux2
    sy2 = _gfilt(y * y, g) - muy2
    sxy = _gfilt(x * y, g) - muxy
    num = (2 * muxy + C1) * (2 * sxy + C2)
    den = (mux2 + muy2 + C1) * (sx2 + sy2 + C2)
    return (num / den).mean(dim=(1, 2, 3))


def calc_mssim(restored, clean) -> float:
    return float(ssim_per_band(restored, clean).mean())


def calc_sam(restored, clean) -> float:
    """calc_SAM.m (度)."""
    tar = _t(restored)
    ref = _t(clean, tar.device)
    prod = (ref * tar).sum(dim=2)
    nrm = torch.sqrt((ref * ref).sum(dim=2) * (tar * tar).sum(dim=2))
    m = nrm != 0
    ang = torch.acos(torch.clamp(prod[m] / nrm[m], -1, 1))
    return float(ang.mean() * 180 / math.pi)
