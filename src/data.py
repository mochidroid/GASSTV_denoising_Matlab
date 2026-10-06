"""HSI の読み込みと観測 (ノイズ付加) の生成.

MATLAB 側の Load_HSI.m / normalize01.m / Generate_obsv_for_denoising.m に対応する.
配列はすべて MATLAB と同じ (n1, n2, n3) = (縦, 横, バンド) の並びで扱う.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import scipy.io

DEFAULT_DATA_DIR = "H:/マイドライブ/MATLAB_Share/HSIData"

# Load_HSI.m の各ケース: (ファイル, 変数名, 切り出し開始位置(1-based), 空間サイズ, バンド指定)
#   bands: None -> start_pos(3) から最後まで, list -> MATLAB の 1-based バンド範囲のリスト
_IMAGE_TABLE = {
    "PaviaU": ("PaviaU/PaviaU.mat", "paviaU", (170, 200, 5), (140, 140), None),
    "PaviaU120": ("PaviaU/PaviaU.mat", "paviaU", (170, 210, 5), (120, 120), None),
    "PaviaU64": ("PaviaU/PaviaU.mat", "paviaU", (211, 211, 5), (64, 64), None),
    "WashingtonDC": ("WashingtonDC/WashingtonDC_image.mat", "WashingtonDC", (657, 140, 1), (100, 100), None),
    "Beltsville": ("Beltsville.mat", "u_org", (145, 157, 1), (100, 100), None),
    "MoffettField128": ("MoffettField.mat", "I_REF", (21, 11, 1), (128, 128), None),
    "MoffettField64": ("MoffettField.mat", "I_REF", (11, 120, 1), (64, 64), None),
    "Salinas": ("Salinas_corrected.mat", "salinas_corrected", (221, 71, 7), (100, 100), None),
}

# JasperRidge64 で使うバンド (MATLAB: [1:102, 110:143, 147:end])
_JR64_BANDS = [(1, 102), (110, 143), (147, None)]


def load_mat(path: str | Path) -> dict[str, np.ndarray]:
    """.mat を読み込む (v5 は scipy, v7.3 は h5py). v7.3 は MATLAB と同じ軸順に戻す."""
    path = Path(path)
    try:
        return {k: v for k, v in scipy.io.loadmat(path).items() if not k.startswith("__")}
    except NotImplementedError:
        out = {}
        with h5py.File(path, "r") as f:
            for k, v in f.items():
                if isinstance(v, h5py.Dataset):
                    arr = v[()]
                    out[k] = arr.T if arr.ndim >= 2 else arr
        return out


def load_dataset(path: str | Path) -> tuple[np.ndarray, np.ndarray, dict, dict]:
    """dataset/<image>/<ノイズ条件>.mat (Matlab/Generate_dataset_denoising.m で作成) を読み込む.

    Returns: HSI_clean, HSI_noisy, hsi, deg (deg にはノイズ条件と各ノイズ成分が入る)
    """
    d = scipy.io.loadmat(path, simplify_cells=True)
    hsi = dict(d["hsi"])
    for k in ("n1", "n2", "n3", "N"):
        hsi[k] = int(hsi[k])
    deg = {k: (float(v) if np.isscalar(v) else v) for k, v in d["deg"].items()}
    return d["HSI_clean"], d["HSI_noisy"], hsi, deg


def normalize01(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    return (x - x.min()) / (x.max() - x.min())


def _select_bands(n3: int, spec) -> np.ndarray:
    idx = []
    for start, end in spec:
        end = n3 if end is None else end
        idx.extend(range(start - 1, end))
    return np.asarray(idx)


def load_hsi(image: str, data_dir: str | Path = DEFAULT_DATA_DIR) -> tuple[np.ndarray, dict]:
    """Load_HSI.m 相当. 戻り値の HSI_clean は [0, 1] 正規化済み (float64)."""
    data_dir = Path(data_dir)
    hsi: dict = {}

    if image in ("JasperRidge", "JasperRidge64"):
        Y = load_mat(data_dir / "JasperRidge/jasperRidge2_R198.mat")["Y"]
        # MATLAB: reshape(Y, [198, 100, 100]) -> permute([2, 3, 1])
        cube = np.reshape(Y, (198, 100, 100), order="F").transpose(1, 2, 0)
        if image == "JasperRidge":
            HSI_clean = normalize01(cube)
        else:
            start = (1, 37, 1)
            size = (64, 64, cube.shape[2] - start[2] + 1)
            hsi["start_pos"] = np.array(start)
            hsi["end_pos"] = np.array(start) + np.array(size) - 1
            sub = cube[start[0] - 1 : start[0] - 1 + 64, start[1] - 1 : start[1] - 1 + 64, :]
            HSI_clean = normalize01(sub[:, :, _select_bands(sub.shape[2], _JR64_BANDS)])
    elif image in _IMAGE_TABLE:
        fname, var, start, size2, _ = _IMAGE_TABLE[image]
        org = load_mat(data_dir / fname)[var]
        size = (*size2, org.shape[2] - start[2] + 1)
        hsi["start_pos"] = np.array(start)
        hsi["end_pos"] = np.array(start) + np.array(size) - 1
        sl = tuple(slice(s - 1, s - 1 + n) for s, n in zip(start, size))
        HSI_clean = normalize01(org[sl])
    else:
        raise ValueError(f"Unknown image: {image}")

    n1, n2, n3 = HSI_clean.shape
    hsi.update(sizeof=np.array([n1, n2, n3]), n1=n1, n2=n2, n3=n3, N=n1 * n2 * n3)
    return HSI_clean, hsi


def _imnoise_salt_pepper(shape, density: float, rng: np.random.Generator) -> np.ndarray:
    """imnoise(0.5*ones(shape), 'salt & pepper', d) 相当: 値は {0, 0.5, 1}."""
    x = rng.random(shape)
    out = np.full(shape, 0.5)
    out[x < density / 2] = 0.0
    out[(x >= density / 2) & (x < density)] = 1.0
    return out


def generate_obsv(HSI_clean: np.ndarray, deg: dict, seed: int = 0) -> tuple[np.ndarray, dict]:
    """Generate_obsv_for_denoising.m 相当 (ストライプ -> ガウス -> スパース -> デッドライン の順に付加).

    乱数生成器が MATLAB と異なるため, ノイズの値そのものは MATLAB と一致しない.
    MATLAB と完全に同じ観測を使いたい場合は Matlab/Generate_dataset_denoising.m で作った dataset/ を使う.
    """
    n1, n2, n3 = HSI_clean.shape
    rng = np.random.default_rng(seed)
    deg = dict(deg)

    # stripe noise
    if deg["stripe_intensity"] > 0:
        sparse_stripe = (
            2 * (_imnoise_salt_pepper((1, n2, n3), deg["stripe_rate"], rng) - 0.5)
            * rng.random((1, n2, n3))
            * np.ones((n1, n2, n3))
        )
        stripe_noise = deg["stripe_intensity"] * sparse_stripe / np.abs(sparse_stripe).max()
        HSI_noisy = HSI_clean + stripe_noise
        deg["stripe_noise"] = stripe_noise
    else:
        HSI_noisy = HSI_clean.copy()

    # Gaussian noise
    if deg["gaussian_sigma"] > 0:
        gaussian_noise = deg["gaussian_sigma"] * rng.standard_normal((n1, n2, n3))
        HSI_noisy = HSI_noisy + gaussian_noise
        deg["gaussian_noise"] = gaussian_noise

    # sparse (salt & pepper) noise
    if deg["sparse_rate"] > 0:
        HSI_tmp = HSI_noisy.copy()
        Sp = _imnoise_salt_pepper((n1, n2, n3), deg["sparse_rate"], rng)
        HSI_noisy[Sp == 0] = 0
        HSI_noisy[Sp == 1] = 1
        deg["sparse_noise"] = HSI_noisy - HSI_tmp

    # dead line noise
    if deg["deadline_rate"] > 0:
        HSI_tmp = HSI_noisy.copy()
        num = int(round(n2 * n3 * deg["deadline_rate"]))
        widths = rng.integers(1, 4, size=num)
        idc = rng.permutation(n2 * n3)[:num]
        mask2d = np.zeros((n2, n3))
        for w, lin in zip(widths, idc):
            col, band = np.unravel_index(lin, (n2, n3), order="F")
            mask2d[col : min(col + w, n2), band] = 1
        mask = np.broadcast_to(mask2d[None], (n1, n2, n3))
        HSI_noisy[mask == 1] = 0
        deg["deadline_mpose"] = HSI_noisy - HSI_tmp

    return HSI_noisy, deg


def obsv_dirname(deg: dict) -> str:
    """MATLAB の append("g", num2str(...), "_ps", ...) と同じフォルダ名."""
    return (
        f"g{num2str(deg['gaussian_sigma'])}_ps{num2str(deg['sparse_rate'])}"
        f"_pt{num2str(deg['stripe_rate'])}_pd{num2str(deg['deadline_rate'])}"
    )


def num2str(x) -> str:
    """MATLAB num2str (スカラ) の簡易版."""
    if isinstance(x, (int, np.integer)) or float(x).is_integer():
        return str(int(x))
    return f"{float(x):.4g}"
