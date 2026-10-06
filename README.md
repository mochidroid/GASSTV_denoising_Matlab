# GASSTV denoising

ハイパースペクトル画像の混合ノイズ除去手法 GASSTV (Condat-Vu 版) の PyTorch 実装と、
元の MATLAB 実装・比較手法 (MATLAB) をまとめたリポジトリ。

## フォルダ構成

```
.
├── run_gasstv.py          提案法 (GASSTV) の実験を回すスクリプト (Python)
├── src/                   提案法の PyTorch 実装
│   ├── solver.py          GASSTV 本体 (Condat-Vu 型 P-PDS)
│   ├── graphs.py          空間グラフ・スペクトルグラフの構築
│   ├── ops.py             差分作用素, prox, 射影
│   ├── metrics.py         MPSNR, MSSIM, SAM
│   ├── data.py            データ読み込み, ノイズ生成
│   └── io.py              .mat 保存
├── dataset/               GT と観測のペア (Matlab/Generate_dataset_denoising.m で作成, git 管理外)
├── result/                実験結果 (git 管理外)
└── Matlab/                MATLAB コード (元の提案法実装と比較手法)
    ├── run_existing_methods.m        比較手法の実験を回すスクリプト
    ├── Generate_dataset_denoising.m  dataset/ を作るスクリプト
    ├── methods/           各手法 (GASSTV_* が提案法, それ以外が比較手法)
    ├── sub_functions/, func_metrics/, Deep/, ...
    └── exp_*.m, output_*.m, ...      従来の実験・図出力スクリプト
```

## 1. データセットの作成 (MATLAB)

MATLAB で `Matlab/Generate_dataset_denoising.m` を実行すると、`Load_HSI` で読み込んだ GT にノイズを重畳したペアが保存される。
観測の作り方は従来の実験 (`exp_denoising_gstd_Base.m`) と同じなので、MATLAB の実験と同一のノイズになる。

```
dataset/<image>/g<..>_ps<..>_pt<..>_pd<..>.mat   (HSI_clean, HSI_noisy, hsi, deg, image)
```

## 2. 提案法 (Python)

```bash
uv sync                                    # 初回のみ
uv run run_gasstv.py                       # 既定設定 (非 Oracle, JasperRidge と PaviaU, ノイズ条件 1)
uv run run_gasstv.py --graph oracle        # Oracle Guide (GT からグラフを構築)
uv run run_gasstv.py --graph noisy oracle --images JasperRidge --noise_idx 1 7
uv run run_gasstv.py --help                # すべてのオプション
```

既定のパラメータは `run_gasstv.py` 冒頭の `DEFAULTS` / `METHOD_PRESETS` に書かれている
(MATLAB の `exp_denoising_gstd_Base.m` と同じ値)。コマンドライン引数で上書きでき、複数の値を渡すと全組み合わせを実行する。

観測は既定 (`--obsv auto`) で `dataset/` から読む。該当ファイルがなければ Python でノイズを生成するが、
乱数生成器が違うためノイズの値は MATLAB と一致しない。

## 3. 比較手法 (MATLAB)

MATLAB で `Matlab/run_existing_methods.m` を実行する (SSTV, HSSTV_L1/L12, l0l1HTV, LRTDTV, TPTV, FastHyMix)。
パラメータは `exp_denoising_gstd_Base.m` と同じ。観測は `dataset/` から読み、結果は Python 版と同じ `result/` に保存する。

従来の `exp_*.m` などのスクリプトも、MATLAB の作業フォルダを `Matlab/` にすればこれまでどおり動く。

既存法のコード (`Matlab/methods/` 以下) には手を加えず、`run_existing_methods.m` 側で次の対応をしている。

- SSTV / HSSTV / l0l1HTV には g, gt 用の関数が無いため、gs, gst 用で代用 (スパースノイズ無しでは alpha = 0 で S = 0 に固定されるので等価)。
- FastHyMix は single 入力だと NaN になるため double で渡す。
  ただし JasperRidge のスパース+ストライプ条件では FastHyMix 本体 (.p) が内部エラーを出し、観測がそのまま返る。

## 結果の保存先

提案法・比較手法とも同じ階層 (従来の MATLAB 実験と同じ構成):

```
result/GASSTV_result/denoising_<image>/g<..>_ps<..>_pt<..>_pd<..>/<method>/<params>.mat       全情報
result/GASSTV_comp_result/denoising_<image>/g<..>_ps<..>_pt<..>_pd<..>/<method>/<params>.mat  復元画像と指標のみ
```

`<method>` は `GASSTV_CondatVu` (提案法, 非 Oracle), `GASSTV_CondatVu_OG` (提案法, Oracle), `SSTV`, `HSSTV_L1`, ...

## Python 版と MATLAB 版 (提案法) の対応

| Python | MATLAB |
|---|---|
| `src/solver.py` | `Matlab/methods/GASSTV_CondatVu/func_GASSTV_{,OG_}{g,gs,gt,gst}_for_denoising_CondatVu.m` |
| `src/graphs.py` | `Create_SpatialGraphWeight.m`, `Create_SpectralGraphLaplacian.m` |
| `src/ops.py` | 差分作用素 (`D4_Neumann_GPU.m` など), `ProjFastL1Ball.m`, `ProjL2ball.m`, `ProxL1norm.m` |
| `src/metrics.py` | `calc_MPSNR.m`, `calc_MSSIM.m`, `calc_SAM.m`, `calc_PSNR_SSIM_per_band.m` |
| `src/data.py` | `Load_HSI.m`, `Generate_obsv_for_denoising.m` |

### 元コードからの修正

- Y2 (空間グラフの L1 ボール制約) の双対更新を Moreau 分解どおり
  `Y2_tmp - gamma2 * Proj(Y2_tmp / gamma2)` に修正している。
  元コードは `Proj(Y2_tmp)` と `/gamma2` が抜けており, gamma2 ≠ 1 (gs/gt/gst) のとき
  制約の内側で ‖W_sp D_sp U‖² の項が余分に加わった問題を解いてしまっていた。
  gamma2 = 1 の g モデルでは元コードと同じ。

### MATLAB 版との一致

上記修正前の状態で、同じ観測とセグメンテーションを与えると MATLAB 版と復元結果が
数値誤差の範囲 (最大差 ~5e-7) で一致することを確認済み。
ほかの差は `imsegkmeans` の代わりに使っている k-means で、ラベルが完全一致するとは限らない
(JasperRidge で一致率 99.8%)。

### オプション

- `--lambda_sp_ref`: 空間グラフ制約の半径 λ_sp を計算する画像。既定の `clean` は元コードと同じで、
  非 Oracle 版でも GT (`HSI_clean`) を使う。
