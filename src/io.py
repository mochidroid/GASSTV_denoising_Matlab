"""結果の保存 (MATLAB で load できる .mat 形式)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.io


def _to_mat(v):
    """scipy.io.savemat で保存できる形に変換する (dict は MATLAB の struct になる)."""
    if v is None:
        return np.zeros((0, 0))
    if isinstance(v, dict):
        return {k: _to_mat(x) for k, x in v.items() if x is not None}
    if isinstance(v, (list, tuple)):
        if all(isinstance(x, str) for x in v):
            return np.array(v, dtype=object)
        return np.asarray(v)
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    return v


def _long_path(path: Path) -> Path:
    """Windows の 260 文字制限を超えるパスでも書き込めるように \\\\?\\ を付ける."""
    path = Path(path).resolve()
    if sys.platform == "win32" and len(str(path)) >= 250 and not str(path).startswith("\\\\?\\"):
        return Path("\\\\?\\" + str(path))
    return path


def save_mat(path: Path, variables: dict) -> None:
    path = _long_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    scipy.io.savemat(path, {k: _to_mat(v) for k, v in variables.items()},
                     long_field_names=True, do_compression=False, oned_as="row")
