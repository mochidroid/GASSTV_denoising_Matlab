"""スペクトルグラフ比較の sweep を段階ごとに実行する (run_spectral_graph_comparison.py を呼ぶだけ).

    uv run sweep_plan.py --stage 1 --dry_run      # 回数の確認
    uv run sweep_plan.py --stage 1                # 実行
    uv run sweep_plan.py --stage 2 --rho 1.0 --lambda_rho_sp 1.0
    uv run sweep_plan.py --stage 3 --rho 1.0 --lambda_rho_sp 1.0 --omega 200
    uv run sweep_plan.py --summary                # 段階ごとの集計

段階:
  1. スペクトルグラフの項が収束解に効く条件を探す (現行法のみ)
       omega = 0 (グラフ項なし) と omega > 0 の差が出る rho, lambda_rho_sp, ノイズ条件を見つける
  2. 全グラフの比較 (1 で見つけた rho, lambda_rho_sp で, omega を振る)
  3. 案 2 (窓) のグラフ構築パラメータを振る
  4. GT からグラフを作った場合 (oracle) との比較 (グラフの推定誤差がどれだけ効いているか)
同じ結果ファイルがあれば飛ばす (--skip_existing を常に付ける) ので, 途中で止めても再実行で続きから回る.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / "run_spectral_graph_comparison.py"

ALL_GRAPHS = ["baseline", "baseline_k10", "abs_h0", "window_h1", "window_h2", "window_h2_unsigned", "coherence"]


def stages(a) -> dict[int, list[list[str]]]:
    return {
        1: [[
            "--graphs", "baseline",
            "--omega", "0", "200", "2000",
            "--rho", *a.rho1,
            "--lambda_rho_sp", *a.lrs1,
            "--noise_idx", *a.noise1,
            "--result_dir", str(ROOT / "result" / "spectral_graph_comparison" / "stage1"),
        ]],
        2: [[
            "--graphs", *ALL_GRAPHS,
            "--omega", *a.omega2,
            "--rho", *a.rho, "--lambda_rho_sp", *a.lambda_rho_sp,
            "--noise_idx", *a.noise,
            "--result_dir", str(ROOT / "result" / "spectral_graph_comparison" / "stage2"),
        ]],
        3: [[
            "--graphs", "window_h2",
            "--omega", *a.omega,
            "--rho", *a.rho, "--lambda_rho_sp", *a.lambda_rho_sp,
            "--noise_idx", *a.noise,
            "--graph_param", "h=1,2,4", "k=10,30,60", "knn=mutual,union", "signed=true,false",
            "--result_dir", str(ROOT / "result" / "spectral_graph_comparison" / "stage3"),
        ]],
        4: [[
            "--graphs", "baseline", "window_h2",
            "--guide", "noisy", "oracle",
            "--omega", *a.omega,
            "--rho", *a.rho, "--lambda_rho_sp", *a.lambda_rho_sp,
            "--noise_idx", *a.noise,
            "--result_dir", str(ROOT / "result" / "spectral_graph_comparison" / "stage4"),
        ]],
    }


SUMMARY_GROUPS = {
    1: ["image", "noise", "rho", "lambda_rho_sp", "omega"],
    2: ["image", "noise", "graph"],
    3: ["image", "noise", "graph_kw"],
    4: ["image", "noise", "graph", "guide"],
}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", type=int, nargs="+", choices=[1, 2, 3, 4])
    p.add_argument("--summary", action="store_true", help="各段階の集計を表示")
    p.add_argument("--images", nargs="+", default=["JasperRidge", "PaviaU"])
    # 段階 1
    p.add_argument("--rho1", nargs="+", default=["0.95", "0.98", "1.02"])
    p.add_argument("--lrs1", nargs="+", default=["0.8", "1.0"])
    p.add_argument("--noise1", nargs="+", default=["1", "9"])
    # 段階 2-4 (段階 1 の結果を見て決める)
    p.add_argument("--rho", nargs="+", default=["0.98"])
    p.add_argument("--lambda_rho_sp", nargs="+", default=["0.8"])
    p.add_argument("--noise", nargs="+", default=["9"])
    p.add_argument("--omega2", nargs="+", default=["20", "200", "2000"])
    p.add_argument("--omega", nargs="+", default=["200"], help="段階 3, 4 の omega")
    a, passthrough = p.parse_known_args()  # --dry_run, --device, --maxiter などはそのまま渡す

    if a.summary:
        for s, grp in SUMMARY_GROUPS.items():
            d = ROOT / "result" / "spectral_graph_comparison" / f"stage{s}"
            if (d / "summary.csv").exists():
                print(f"\n##### stage {s} #####")
                subprocess.run([sys.executable, str(RUNNER), "--summary_only", "--result_dir", str(d),
                                "--group_by", *grp], check=False)
        return
    if not a.stage:
        p.error("--stage か --summary を指定してください")

    plan = stages(a)
    for s in a.stage:
        for extra in plan[s]:
            cmd = [sys.executable, str(RUNNER), "--images", *a.images, *extra, "--skip_existing", *passthrough]
            print(f"\n##### stage {s} #####\n" + " ".join(cmd), flush=True)
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
