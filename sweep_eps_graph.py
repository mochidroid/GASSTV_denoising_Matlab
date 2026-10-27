"""案 2 (窓) のスペクトルグラフで, 疎化を mutual kNN から eps 近傍グラフに替えた比較 (rho, lambda_rho_sp は固定).

    uv run sweep_eps_graph.py --dry_run     # 実行回数
    uv run sweep_eps_graph.py               # 実行 (止めても再実行で続きから)
    uv run sweep_eps_graph.py --summary     # 集計だけ

比べるもの (rho = 0.95, lambda_rho_sp = 0.8, omega = 10, 30, 100, 300, 1000):
  - 案 2 + eps 近傍  : eps_c = 0.5, 0.75, 1.0 (k = 0: 次数の上限なし) x h = 1, 2 x 符号 (exact / none)
  - 案 2 + mutual kNN: k = 30, 60 x h = 1, 2 x 符号 (exact / none)       ... 参照
  - 現行グラフ, omega = 0 (スペクトルグラフ項なし)                         ... 参照
eps_c は距離の閾値 (d_ij < eps_c * sigma, sigma は距離の中央値). eps_c = 1 で約半数のペアが残る.
グラフ単体の検証では eps_c <= 0.3 は疎すぎて (孤立ノードが多く) 代理デノイジングが大きく悪化したので範囲外にしている.

結果: result/spectral_graph_comparison/eps_graph/
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / "run_spectral_graph_comparison.py"
OUT = ROOT / "result" / "spectral_graph_comparison" / "eps_graph"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--images", nargs="+", default=["JasperRidge", "PaviaU"])
    p.add_argument("--noise_idx", nargs="+", default=["9"])
    p.add_argument("--rho", default="0.95")
    p.add_argument("--lambda_rho_sp", default="0.8")
    p.add_argument("--omega", nargs="+", default=["10", "30", "100", "300", "1000"])
    p.add_argument("--eps_c", default="0.5,0.75,1.0")
    p.add_argument("--dry_run", action="store_true")
    p.add_argument("--summary", action="store_true")
    a, passthrough = p.parse_known_args()

    common = [sys.executable, str(RUNNER), "--images", *a.images, "--noise_idx", *a.noise_idx,
              "--rho", a.rho, "--lambda_rho_sp", a.lambda_rho_sp, "--result_dir", str(OUT),
              "--skip_existing", *passthrough]
    runs = [
        # 案 2 + eps 近傍
        ["--graphs", "p2_window", "--omega", *a.omega,
         "--graph_param", "knn=eps", f"eps_c={a.eps_c}", "k=0", "h=1,2", "sign_mode=exact,none"],
        # 案 2 + mutual kNN (参照)
        ["--graphs", "p2_window", "--omega", *a.omega,
         "--graph_param", "knn=mutual", "k=30,60", "h=1,2", "sign_mode=exact,none"],
        # 現行グラフと omega = 0 (参照)
        ["--graphs", "baseline", "--omega", "0", *a.omega],
    ]
    if not a.summary:
        for r in runs:
            cmd = common + r + (["--dry_run"] if a.dry_run else [])
            print("\n$ " + " ".join(cmd), flush=True)
            subprocess.run(cmd, check=True)
        if a.dry_run:
            return
    for grp in (["image", "graph", "graph_kw"], ["image", "graph", "omega"]):
        subprocess.run([sys.executable, str(RUNNER), "--summary_only", "--result_dir", str(OUT),
                        "--group_by", *grp], check=False)


if __name__ == "__main__":
    main()
