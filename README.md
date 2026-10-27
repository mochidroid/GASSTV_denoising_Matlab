# GASSTV denoising

ハイパースペクトル画像の混合ノイズ除去手法 GASSTV (Condat-Vu 版) の PyTorch 実装と、
元の MATLAB 実装・比較手法 (MATLAB) をまとめたリポジトリ。

## フォルダ構成

```
.
├── run_gasstv.py          提案法 (GASSTV) の実験を回すスクリプト (Python)
├── verify_spectral_graphs.py         スペクトルグラフ単体の検証 (GT から作ったグラフと noisy から作ったグラフの比較)
├── run_spectral_graph_comparison.py  スペクトルグラフの作り方を変えて GASSTV 本体を回す (パラメータ sweep 対応)
├── sweep_all.py           各案のパラメータを広く振る sweep (A: モデル -> B: グラフ構築 -> C: omega の振り直し)
├── compare_methods.py     既存法 (MATLAB) と GASSTV 系 (Python) の結果を rho ごとに 1 つの表にする
├── sweep_eps_graph.py     案2 の疎化を mutual kNN から eps 近傍グラフに替えた比較
├── sweep_noise_conditions.py  ノイズ条件 1-9 でのスペクトルグラフの比較
├── src/                   提案法の PyTorch 実装
│   ├── solver.py          GASSTV 本体 (Condat-Vu 型 P-PDS)
│   ├── graphs.py          空間グラフ・スペクトルグラフ (現行法) の構築
│   ├── spectral_graphs.py スペクトルグラフの候補手法 (現行法 / 案1 絶対値 / 案2 窓 / 案3 コヒーレンス / 提案法 anchored)
│   ├── fast_ops.py        CPU 実行の高速化 (L1 ボール射影の Michelot 法, torch.compile. 結果は同一)
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

## 0. 環境構築 (Python)

Python 側はパッケージ管理に [uv](https://docs.astral.sh/uv/) を使う。
Python 本体 (3.12, `.python-version` で指定) と依存パッケージ (CUDA 12.8 版 PyTorch など, `pyproject.toml`) は uv がまとめて用意するので,
事前に Python を入れておく必要はない。

### 0.1 uv のインストール

Windows (PowerShell) では次のどちらか:

```powershell
# 公式インストーラ
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

# または winget
winget install --id=astral-sh.uv -e
```

インストール後, **ターミナルを開き直してから** 確認する:

```powershell
uv --version
```

`uv` が見つからないと言われる場合は, ターミナル (VS Code なら VS Code ごと) を再起動する。
それでもだめなら `%USERPROFILE%\.local\bin` を環境変数 PATH に追加する。

macOS / Linux の場合は `curl -LsSf https://astral.sh/uv/install.sh | sh`。

### 0.2 GPU (NVIDIA ドライバ) の確認

```powershell
nvidia-smi
```

右上の `CUDA Version` が 12.8 以上なら, `pyproject.toml` で指定している CUDA 12.8 版 PyTorch が動く。
これより古い場合は NVIDIA ドライバを更新する (CUDA Toolkit 自体のインストールは不要。PyTorch に同梱されている)。

### 0.3 依存パッケージのインストール

リポジトリのフォルダで:

```powershell
cd C:\Users\shint\Documents\Github\GASSTV_denoising_Matlab
uv sync
```

初回は Python 3.12 と PyTorch (CUDA 版, 数 GB) をダウンロードするので時間がかかる。`.venv/` が作られ, 以降の `uv run ...` はこの環境で動く
(`activate` は不要)。

### 0.4 動作確認

```powershell
uv run python -c "import torch; print(torch.__version__, torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else '-')"
```

`True` と GPU 名 (例: `NVIDIA GeForce RTX 3090`) が出れば OK。`False` のときは 0.2 のドライバを確認する。

### 0.5 (任意) 図を描くパッケージ

`verify_spectral_graphs.py` の図の出力には matplotlib が要る:

```powershell
uv add matplotlib
```

(`pyproject.toml` と `uv.lock` に追記される。図が不要なら `--no_plot` を付ければ入れなくてよい。)

### よくあるつまずき

| 症状 | 対処 |
|---|---|
| `uv : 用語 'uv' は ... 認識されません` | ターミナルを開き直す。だめなら PATH に `%USERPROFILE%\.local\bin` を追加 |
| インストーラ実行時に「スクリプトの実行が無効」 | 上のコマンドの `-ExecutionPolicy ByPass` を付けて実行しているか確認 (または winget を使う) |
| `torch.cuda.is_available()` が `False` | `nvidia-smi` の CUDA Version を確認しドライバを更新。`uv sync` をやり直す |
| `uv sync` が途中で失敗する (ネットワーク) | もう一度 `uv sync` (ダウンロード済みの分はキャッシュされる) |
| `.venv` が壊れた / 別 PC からコピーした | `.venv` フォルダを削除して `uv sync` |

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

MATLAB で `Matlab/run_existing_methods.m` を実行する。既定は rho を使う 4 手法 (SSTV, HSSTV_L1, HSSTV_L12, l0l1HTV),
rho ∈ {0.90, 0.95}, ノイズ条件 9, 両画像, 結果があれば飛ばす。リポジトリ直下で:

```powershell
matlab -batch "run('Matlab/run_existing_methods.m')"
# 設定の上書き (指定しない項目は既定値)
matlab -batch "cfg.noise = [1 9]; cfg.rhos = {0.90}; cfg.methods = {'SSTV','HSSTV_L1'}; run('Matlab/run_existing_methods.m')"
```

`cfg` の項目: `noise` (ノイズ条件の番号), `images` (画像の番号), `rhos` (cell), `methods` (SSTV, HSSTV_L1, HSSTV_L12, l0l1HTV, LRTDTV, TPTV, FastHyMix),
`skip_existing`。MATLAB の GUI から回すときは, ファイル冒頭の既定値を書き換えて実行してもよい。
パラメータは `exp_denoising_gstd_Base.m` と同じ。観測は `dataset/` から読み、結果は Python 版と同じ `result/` に保存する。

従来の `exp_*.m` などのスクリプトも、MATLAB の作業フォルダを `Matlab/` にすればこれまでどおり動く。

既存法のコード (`Matlab/methods/` 以下) には手を加えず、`run_existing_methods.m` 側で次の対応をしている。

- SSTV / HSSTV / l0l1HTV には g, gt 用の関数が無いため、gs, gst 用で代用 (スパースノイズも deadline ノイズも無いときは alpha = 0 で S = 0 に固定されるので等価)。
- epsilon, alpha, beta は `Matlab/sub_functions/calc_noise_radii.m` (Python は `src/data.py` の `noise_radii`) で計算し、deadline ノイズはスパースノイズの項に含める。
- FastHyMix は single 入力だと NaN になるため double で渡す。
  ただし JasperRidge のスパース+ストライプ条件では FastHyMix 本体 (.p) が内部エラーを出し、観測がそのまま返る。

## 4. スペクトルグラフの作り方の比較

スペクトルグラフ (隣接バンド差分をノードとするグラフ) の作り方を変えて比較する。候補 (`src/spectral_graphs.py`):

| 名前 | 内容 |
|---|---|
| `baseline` | 現行法. セグメントごとの代表スペクトル差分の値の類似度 (符号なしラプラシアン) |
| `p1_abs` (= `abs_h0`) | 案1. 代表スペクトル差分の絶対値の類似度 + 符号付きラプラシアン (画像全体で 1 つ) |
| `p2_window` (= `window_h2`) | 案2. 差分の絶対値の窓 (半幅 h) の類似度 + 符号付きラプラシアン |
| `p3_coherence` (= `coherence`) | 案3. 差分画像どうしの低空間周波数コヒーレンス |

### 4.1 グラフ単体の検証

```powershell
uv run verify_spectral_graphs.py            # 結果: result/graph_verification/
```

GT から作ったグラフと noisy から作ったグラフについて, 一致度, GT のスペクトルの滑らかさ (レイリー商),
スペクトルグラフだけでの代理デノイジング, 急峻な変動どうしがつながっているか, を出す。

### 4.2 GASSTV 本体での比較 (sweep)

```powershell
uv run sweep_all.py --dry_run    # 実行回数の確認
uv run sweep_all.py              # A -> B -> C を順に実行 (止めても再実行で続きから)
uv run sweep_all.py --report     # 結果の表 (result/spectral_graph_comparison/sweep_all/report.md)
```

- A: 各案の既定グラフで, モデルのパラメータ (omega, rho, lambda_rho_sp, sigma_sp) を全組み合わせ
- B: 各案・各画像で A の最良 (omega > 0) に固定し, グラフ構築のパラメータを全組み合わせ
- C: 最良のグラフ設定で omega を最良値のまわりで振り直す

グリッドは `sweep_all.py` 冒頭の `PROFILES` で変える。`--levels A`, `--proposals p2_window`, `--images PaviaU` で絞れる。
動作確認は `uv run sweep_all.py --profile quick --maxiter 50`。

個別に回すときは `run_spectral_graph_comparison.py` を直接使う (複数の値はすべて全組み合わせで実行):

```powershell
uv run run_spectral_graph_comparison.py --graphs baseline p2_window --omega 0 100 1000 --rho 0.95 0.98 --lambda_rho_sp 0.8 1.0
uv run run_spectral_graph_comparison.py --graphs p2_window --omega 1000 --graph_param h=1,2,4 k=10,30 signed=true,false
uv run run_spectral_graph_comparison.py --summary_only --group_by image graph omega
```

主なパラメータ:

| 引数 | 意味 |
|---|---|
| `--omega` | スペクトルグラフ項の重み. 既定 (`--lambda2_mode normalized`) では L を最大固有値で割ってから掛ける. 現行法の raw の lambda2 とは `lambda2 ≈ omega / lam_max` (lam_max は JasperRidge で約 180, PaviaU で約 80) |
| `--rho` | epsilon (L2 ボール), alpha (スパース), beta (ストライプ) の半径の倍率 |
| `--lambda_rho_sp` | 空間グラフ L1 ボールの半径の倍率 |
| `--sigma_sp` | 空間グラフの重みの幅 (`90` / `med` / 数値) |
| `--guide` | グラフを作る画像 (`noisy` / `oracle`) |
| `--graph_param KEY=V1,V2` | グラフ構築のパラメータ. 現行法: `num_segments`, `sigma_l`, `k_lap`, `order_filt` / 案1・2: `h`, `feature`, `eps_q`, `sigma`, `k`, `knn`, `w_min`, `tau_q`, `signed`, `rep`, `order_filt` / 案3: `patch`, `step`, `f_max`, `exclude_col_axis`, `prefilter`, `min_sep`, `k`, `knn`, `w_min`, `signed` |

停止条件は `run_gasstv.py` と同じで, U の相対変化 `||U_{k+1} - U_k|| / ||U_k||` が `10^-stopcri_idx` (既定 1e-5) を下回るか,
`--maxiter` (既定 20000) に達したとき。omega が大きいと U のステップ幅が約 1/(2 omega) まで小さくなるので,
解が動いている途中でも相対変化が小さくなって止まることがある点に注意 (report の † は最大反復数で止まった実行)。

CPU で回す場合は L1 ボール射影が自動で高速版になる。さらに `--compile` を付けると約 1.4 倍速い (結果はどちらも同一)。

### 4.3 提案法: 代表スペクトルにつないだスペクトルグラフ

セグメント $s$ の代表スペクトルの差分を $\delta_s$, 各画素の差分のずれを $e_p = D_\lambda u_p - \delta_s$ として

$$R(u) = \sum_p e_p^\top (L_s + \beta_q Q_s)\, e_p$$

- $Q_s$ (自己ループ): 各ノードを値 $\delta_{s,i}$ に固定した代表ノードへつなぐ. 「同じセグメントのスペクトルは代表スペクトルに等倍で近づく」
- $L_s$ (辺): $\delta_s$ の値が近いノードを結ぶ (現行法と同じ作り方). 定数倍のずれ $e = (c-1)\delta_s$ はこのグラフ上で滑らかなので罰が小さい
- 明るさのオフセットは差分で消えるので罰しない

`run_spectral_graph_comparison.py` のグラフ名: `anchored` (辺 + 自己ループ), `anchor_only` (自己ループのみ), `residual_edges` (ずれに辺のみ)。
主な `--graph_param`: `beta_q` (自己ループの強さ), `seg` (`intensity` = 平均輝度の k-means / `spectral` = スペクトル特徴量の k-means),
`num_segments`, `rep` (`median` / `colmed`), `edges` (`delta` / `none`), `anchor` (`uniform` / `dispersion` / `none`), `k`。

### 4.4 rho の決め方 (既存法との比較)

rho は既存法も含めて全手法で共通の値 (0.90 か 0.95) に固定する。両方で回して比べてから決める:

```powershell
# 1. 既存法 (MATLAB): SSTV, HSSTV_L1, HSSTV_L12, l0l1HTV を rho = 0.90, 0.95 で
matlab -batch "run('Matlab/run_existing_methods.m')"

# 2. GASSTV 系 (Python): 現行グラフ, 案2, 提案法 (と比較実験) を rho = 0.90, 0.95, lambda_rho_sp = 0.6-0.9, omega = 0-1000 で
uv run sweep_all.py --profile rho_select --levels A --dry_run   # 416 回
uv run sweep_all.py --profile rho_select --levels A

# 3. 表にする (result/method_comparison.md)
uv run compare_methods.py
```

1 と 2 は別のターミナルで同時に回してよい (GPU を取り合うので遅くはなる)。
rho を決めたら, 提案法のグラフ構築パラメータを詰める:

```powershell
uv run sweep_all.py --profile rho_select --levels B C --proposals anchored anchor_only residual_edges
```

(B は A の結果のうち各案で最良のモデルパラメータを使う。rho を 1 つに固定したいときは `sweep_all.py` の `PROFILES["rho_select"]` の `rho` を 1 つにしてから A を回し直す。)

### 4.5 eps 近傍グラフ (案2 の疎化の変更)

kNN は全ノードに k 本の辺を強制するので, 似た変動を持つ相手がいない (急峻な) ノードも結ばれてしまう。
eps 近傍グラフは距離 $d_{ij} < \varepsilon_c \sigma$ のペアだけを結び ($\sigma$ は距離の中央値), 次数はデータで決まる。

```powershell
uv run sweep_eps_graph.py --dry_run   # 212 回 (rho = 0.95, lambda_rho_sp = 0.8 で固定)
uv run sweep_eps_graph.py
uv run sweep_eps_graph.py --summary
```

個別に使うときは `--graph_param knn=eps eps_c=0.75 k=0` (k > 0 にすると次数の上限として mutual kNN も併用)。
符号は `sign_mode` で選ぶ: `threshold` (従来. |delta| が中央値未満のノードの符号は未定) / `exact` (sgn(0)=+1 で全ノードに符号.
L = Sigma L_unsigned Sigma となる平衡な符号付きグラフ) / `none` (符号なし)。

### 4.6 ノイズ条件ごとの検証

rho = 0.95, lambda_rho_sp = 0.8 で固定し, 各条件・各画像で omega = 10-1000 を振って最良値を比べる
(スペクトルグラフ項なし / 現行グラフ / 案2 mkNN 符号あり / 案2 mkNN 符号なし / 案2 eps 符号あり)。

```powershell
uv run sweep_noise_conditions.py --dry_run   # 378 回
uv run sweep_noise_conditions.py
uv run sweep_noise_conditions.py --report    # result/spectral_graph_comparison/noise_conditions/report.md

# 既存法も同じ条件で (MATLAB) -> 1 つの表に
matlab -batch "cfg.noise = 1:9; cfg.rhos = {0.95}; run('Matlab/run_existing_methods.m')"
uv run compare_methods.py --sweeps noise_conditions --rhos 0.95
```

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
