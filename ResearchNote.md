# Research Note: GASSTV のスペクトルグラフ設計 (2026-10)

最終更新: 2026-10-27。1 か月後に読み返して「何を・なぜ・どこまで」やったかが分かるように書いた記録。
詳細な数値は各 report / CSV にあり、ここには要点だけを書く。

---

## 0. 3 行でいうと

- GASSTV (SSTV + 空間グラフ制約 + スペクトルグラフ正則化) の **スペクトルグラフの作り方** を変えて比べた。
  - 現行グラフ (`baseline`): セグメントごとに、代表スペクトル差分の **値** が近いバンドをつなぐ。
  - 案2 (`p2_window`): 代表スペクトル差分の **絶対値の窓** が近いバンドを、**符号付き** でつなぐ。
- 採用候補は **案2 mutual kNN (k=60, h=1, 符号 exact) = `p2_mknn`**。
  - 現行グラフに負けない (最良 ω どうしで 12 勝 2 分 0 敗)。
  - ω を大きくしたときの劣化が小さい。
  - ただし差が 0.05 dB を超えたのは 14 組中 4 組だけ (σ=0.1 でスパースノイズを含む条件)。
- 未完: ノイズ条件 1 と 8 の GASSTV 系、条件 8 の既存法 (result/ 側)、現行の式での条件 9 の回し直し、複数のノイズ実現での確認 (9 節)。

---

## 1. 背景と目的

- 対象: ハイパースペクトル画像 (HSI) の混合ノイズ除去。混合ノイズはガウス、スパース (salt & pepper)、ストライプ、デッドラインの 4 種類。
- 手法 GASSTV は、元は MATLAB 実装 (`Matlab/methods/GASSTV_CondatVu/`)。それを PyTorch に移植して (`src/solver.py`)、スペクトルグラフの設計を検証している。
- 動機: 現行のスペクトルグラフには次の 2 つの問題があるのではないか。
  1. 差分の **値** の類似度でつなぐので、向きが逆の急峻な変動 (吸収帯の立ち下がりと立ち上がりなど) どうしは結ばれない。
  2. 平坦なバンドと急峻なバンドが、ガウス重みの裾でつながる。
- これを、
  - (i) 変動の **大きさ** (|δ|) でつなぎ、
  - (ii) **符号** で向きの違いを表し、
  - (iii) 前後のバンドを含む **窓** で比べる、
  
  ように作り替えたのが案2。

---

## 2. データとノイズ条件

| 画像 | サイズ (n1 × n2 × n3) | N |
|---|---|---|
| JasperRidge | 100 × 100 × 198 | 1,980,000 |
| PaviaU | 140 × 140 × 99 | 1,940,400 |

- 観測は `dataset/<image>/<ノイズ条件>.mat`。`Matlab/Generate_dataset_denoising.m` で作った GT と観測のペアで、MATLAB の `rng("default")` を使う。
  - Python 側も MATLAB 側も同じ観測を読む。**ノイズの実現は 1 通り**。
- ノイズ条件 (`run_gasstv.NOISE_CONDITIONS`、MATLAB 側も同じ番号):

| 番号 | フォルダ名 | 内容 |
|---|---|---|
| 1 | g0.1_ps0_pt0_pd0 | ガウス σ=0.1 |
| 2 | g0.05_ps0.05_pt0_pd0 | ガウス 0.05 + スパース 5% |
| 3 | g0.1_ps0.05_pt0_pd0 | ガウス 0.1 + スパース 5% |
| 4 | g0.05_ps0_pt0.05_pd0 | ガウス 0.05 + ストライプ 5% (強度 0.5) |
| 5 | g0.1_ps0_pt0.05_pd0 | ガウス 0.1 + ストライプ |
| 6 | g0.05_ps0.05_pt0.05_pd0 | ガウス 0.05 + スパース + ストライプ |
| 7 | g0.1_ps0.05_pt0.05_pd0 | ガウス 0.1 + スパース + ストライプ |
| 8 | g0.1_ps0_pt0_pd0.01 | ガウス 0.1 + デッドライン 1% |
| 9 | g0.1_ps0.05_pt0.05_pd0.01 | 全部 |

- デッドライン: (列, バンド) の組をランダムに選び、その列から幅 1〜3 列 (一様、平均 2 列) を、全行にわたって 0 にする。

---

## 3. 最適化問題 (GASSTV)

### 3.1 定式化

観測を $V \in \mathbb{R}^{n_1 \times n_2 \times n_3}$、復元画像を $U$、スパースノイズを $S$、ストライプノイズを $T$ とする。

$$
\begin{aligned}
\min_{U,S,T}\;& \|\mathbf{D}\,\mathbf{D}_\lambda U\|_1 \;+\; \frac{\omega}{2}\sum_{p} (\mathbf{D}_\lambda u_p)^\top \tilde L_{\mathrm{seg}(p)} (\mathbf{D}_\lambda u_p) \\
\text{s.t.}\;& \|W_{sp}\odot \mathbf{D}_4 U\|_1 \le \lambda_{sp}, \quad \|S\|_1 \le \alpha, \quad \|T\|_1 \le \beta, \quad \mathbf{D}_v T = 0,\\
& \|U + S + T - V\|_2 \le \varepsilon, \quad U \in [0,1]^{n_1 n_2 n_3}
\end{aligned}
$$

| 記号 | 意味 | コード (`src/ops.py`) |
|---|---|---|
| $\mathbf{D}_\lambda$ | 波長方向の前進差分 (バンド $i+1$ − バンド $i$, Neumann 境界) | `Dl`, `Dl_glr` |
| $\mathbf{D}$ | 空間の縦・横の前進差分。$\|\mathbf{D}\mathbf{D}_\lambda U\|_1$ が SSTV | `D` |
| $u_p \in \mathbb{R}^{n_3}$ | 画素 $p$ のスペクトル。$\mathbf{D}_\lambda u_p \in \mathbb{R}^{K}$, $K = n_3 - 1$ | |
| $\tilde L_s$ | セグメント $s$ のスペクトルグラフのラプラシアン ($K \times K$)。4 節 | |
| $\mathbf{D}_4$, $W_{sp}$ | 空間 4 方向 (下・右・右下・右上) の差分と、その重み | `D4` |
| $\mathbf{D}_v$ | 縦方向の差分。ストライプは列ごとに縦に一定、という制約 | `Dv` |

- 第 1 項 (SSTV): スペクトル差分の空間的な滑らかさ。
- 第 2 項 (GLR, graph Laplacian regularization): 各画素のスペクトル差分が、スペクトルグラフ上で滑らかであること。**今回変えたのはここだけ。**
- 空間グラフ制約 (L1 ボール): ガイド画像のエッジでは重みを小さくして、エッジを保つ。

### 3.2 ノイズモデル (どの変数を使うか)

`select_noise_model` (`src/solver.py`) で決める。

| モデル | 使う変数 | 条件 |
|---|---|---|
| g | U | スパースもデッドラインもストライプも無い (条件 1) |
| gs | U, S | スパースまたはデッドラインがあり、ストライプが無い (2, 3, 8) |
| gt | U, T | ストライプのみ (4, 5) |
| gst | U, S, T | 両方 (6, 7, 9) |

デッドラインはスパースの項 S で扱う (2026-10-20 に変更、3.4 節)。既存法 (SSTV / HSSTV / l0l1HTV) は g, gt 用の関数が無いので、gs, gst 用で代用している。

### 3.3 パラメータ (半径と重み)

- **ω (omega)**: GLR の重み。グラフ間で意味を揃えるため、$\tilde L = L / \lambda_{\max}(L)$ と正規化してから掛ける (`lambda2_mode=normalized`)。
  - 現行の raw の lambda2 との関係は $\text{lambda2} \approx \omega/\lambda_{\max}$。$\lambda_{\max}$ は JasperRidge で約 180、PaviaU で約 80。
- **λ_sp** = `lambda_rho_sp` × $\|W_{sp}\odot \mathbf{D}_4 U_{\mathrm{GT}}\|_1$。
  - **既定では GT から計算する** (`--lambda_sp_ref clean`、元の MATLAB コードと同じ)。非 Oracle 版でも GT の情報を 1 つ使っている点に注意。
- **W_sp**: バンド平均画像 $g$ の 4 方向差分 $\nabla g$ に対し、$w = \exp(-(\nabla g)^2 / 2\sigma_{sp}^2)$。
  - `sigma_sp="90"` のとき: $|\nabla g|$ の 90 パーセンタイル $g_0$ で重みが 0.2 になるように、$\sigma_{sp} = g_0/\sqrt{2\ln 5}$ とする。
- **ε, α, β** (`rho_radius` = ρ 倍)。式は `calc_noise_radii.m` / `src/data.py::noise_radii`。
  - デッドラインの被覆率: $c_{dl} = 1 - e^{-2 p_d}$
  - $\mu$: 観測の (0, 1) 内の画素値の平均
  - $p_s$: スパース率、$p_t$: ストライプ率、$t$: ストライプ強度、$\sigma$: ガウスの標準偏差

$$
\varepsilon = \rho\,\sigma\sqrt{N(1-p_s)(1-c_{dl})},\qquad
\alpha = \rho\,N\,(0.5\,p_s + \mu\,c_{dl}),\qquad
\beta = \rho\,N(1-p_s)(1-c_{dl})\,\frac{p_t\, t}{2}
$$

### 3.4 半径の式の経緯 (重要)

| 時期 | alpha の式 | 使った実験 |
|---|---|---|
| 2025-09 (commit `4a66b09`) | 上の式 (デッドラインをスパースの項に含める) | 過去の MATLAB 実験 (`H:\マイドライブ\...`) の条件 8, 9 |
| 2025-12 〜 2026-10-20 | $\alpha = \rho \cdot 0.5 N (p_s - 2 p_s p_d)$。デッドラインを含めない | Python の sweep_all, rho_select, eps_graph、既存法の result/ の条件 9 |
| 2026-10-20 〜 現在 | 上の式に戻した (MATLAB と Python の両方) | これから回す実験 |

- デッドラインが無い条件 (1-7) では、2 つの式は一致する。
- **条件 9 の既存の結果は、すべて旧式で計算されている。**
- 条件 8 は旧式だと alpha = 0 になり、デッドラインを除けない。

### 3.5 アルゴリズム

Condat–Vũ 型の P-PDS (`src/solver.py`)。

- 主変数 U, S, T の更新: 勾配 (GLR) と、双対変数からの寄与を足してから射影する (U は箱 [0, 1]、S と T は L1 ボール)。
- 双対変数 Y1 (SSTV の L1)、Y2 (空間グラフの L1 ボール)、Y3 (Dv T = 0)、Y4 (L2 ボール) の更新: Moreau 分解で行う。
- ステップ幅:
  - $\gamma_{1,U} = 1/(\tfrac{1}{2}\cdot 4\omega\lambda_{\max}(\tilde L) + 8\cdot 4 + 16 + 1)$
  - $\gamma_{1,S} = 1$
  - $\gamma_{1,T} = 1/5$
  - $\gamma_2$ は g で 1、gs と gt で 1/2、gst で 1/3
- 停止: $\|U_{k+1}-U_k\|/\|U_k\| < 10^{-5}$、または 20000 反復。
  - ω が大きいと U のステップ幅が約 $1/(2\omega)$ まで小さくなり、解が動いている途中で止まることがある。
- 元の MATLAB コードからの修正: Y2 の双対更新を $Y_2^{tmp} - \gamma_2\,\mathrm{Proj}(Y_2^{tmp}/\gamma_2)$ にした。元のコードでは $\gamma_2 \neq 1$ のとき別の問題を解いていた。
  - 修正前の状態では MATLAB 版と最大差 5e-7 で一致することを確認済み。

---

## 4. スペクトルグラフの設計

### 4.1 共通の枠組み

- ノード: 隣接バンド差分 $i = 1,\dots,K$ ($K = n_3 - 1$)。グラフは $K \times K$ で、全画素が共有する (現行グラフはセグメントごとに 1 つ)。
- 正則化 $\frac{\omega}{2}\sum_p z_p^\top \tilde L z_p$, $z_p = \mathbf{D}_\lambda u_p$ は、つながったバンドどうしで差分が近くなるよう促す。
- 符号付きラプラシアン (`src/spectral_graphs.py`):

$$
L = \mathrm{diag}(W\mathbf 1) - W \odot \Sigma, \qquad
x^\top L x = \sum_{i<j} w_{ij}\,(x_i - \Sigma_{ij}\,x_j)^2
$$

- $\Sigma_{ij} = \sigma_i \sigma_j$ (ノードの符号 $\sigma_i = \mathrm{sgn}(\delta_i)$) のとき、$L = \mathrm{diag}(\sigma)\,L_{\text{unsigned}}\,\mathrm{diag}(\sigma)$ となる (平衡な符号付きグラフ)。
  - 固有値は符号なしと同じ。
  - 向きが逆の急峻な変動 ($\delta_i > 0, \delta_j < 0$) に対しては、$x_i \approx -x_j$、つまり「同じ大きさで逆向き」であることを促す。

### 4.2 現行グラフ `baseline` (`src/graphs.py::spectral_graph_laplacian`)

1. バンド平均画像を 1 次元 k-means で S=4 個のセグメントに分ける。
2. 各セグメントの代表スペクトル $r_s$ = 画素の中央値に、波長方向のメディアンフィルタ (次数 5、ゼロ埋め) をかけたもの。$\delta_s = \mathrm{diff}(r_s)$。
3. 重み $w_{ij} = \exp(-(\delta_{s,i}-\delta_{s,j})^2 / 2\sigma_l^2)$。
   - $\sigma_l = \sqrt{\mathrm{median}(d^2)/S}$ (`sigma_l="med"`)。$d^2$ は全セグメントの差分ベクトルから作るバンド間距離。
4. `k_lap=1000` ≥ K−1 なので全結合。符号なし。

問題意識: 値の類似度なので、急峻な変動どうしでも向きが逆なら遠い。向きが同じでも大きさが違えば遠い。平坦なバンドどうしは全部近い。

### 4.3 案1 `p1_abs` と案2 `p2_window` (`graph_abs_signed`)

画像全体で 1 つのグラフを作る (セグメントなし)。

1. 代表スペクトル $r$ (`rep="colmed"`): 列ごとに中央値をとってから、列方向に中央値をとる。ストライプとデッドラインは列単位で乗るので、外れ値の列として除かれる。その後、メディアンフィルタ (次数 5、端は nearest)。$\delta = \mathrm{diff}(r)$。
2. 特徴量 (`feature="log"`): $f_i = \log(|\delta_i| + \epsilon)$。$\epsilon$ は $|\delta|$ の中央値。
   - 変動の大きさを比で比べる。ε 未満の平坦なバンドは 1 つの塊になる。
3. 窓 (案2): ノード $i$ の特徴ベクトルを $F_i = (f_{i-h},\dots,f_{i+h})$ とする (端は鏡映)。h = 0 なら案1。
   - 距離 $d_{ij}^2 = \|F_i - F_j\|^2/(2h+1)$。
   - $\sigma$ = 距離の中央値。$w_{ij} = \exp(-d_{ij}^2/2\sigma^2)$。
4. 疎化 (`_sparsify`):
   - `knn=mutual`: 互いに k 近傍のときだけ残す (min で対称化)。採用は k = 60。
   - `knn=eps`: $d_{ij} < \epsilon_c\sigma$ のペアだけ残す。次数はデータで決まり、似た相手の無いノードは孤立する。
5. 符号 (`sign_mode`):
   - `exact`: $\sigma_i = \mathrm{sgn}(\delta_i)$ ($\mathrm{sgn}(0) = +1$)。
   - `threshold`: $|\delta|$ の中央値未満のノードの符号は未定 (+1 扱い)。初期の既定。
   - `none`: 符号なし。

### 4.4 その他の候補 (途中で外したもの)

- 案3 `p3_coherence`: 差分画像どうしの、低空間周波数のコヒーレンスでつなぐ。
  - パッチ 20 で Welch 平均。縦方向に一定な成分 (ストライプ) は除く。
  - 隣接差分はノイズが相関するので、min_sep = 2 で除く。
- 提案法 `anchored`: $R(u) = \sum_p (\mathbf{D}_\lambda u_p - \delta_s)^\top (L_s + \beta_q Q_s)(\mathbf{D}_\lambda u_p - \delta_s)$。
  - 自己ループ $Q_s$ で各ノードを代表値に固定する。
  - 変種として `anchor_only` (自己ループのみ) と `residual_edges` (ずれに辺のみ)。

### 4.5 グラフ単体の検証 (`verify_spectral_graphs.py`, 条件 9)

GASSTV を回す前に、グラフだけで良し悪しを見た。

- Rc: GT のスペクトル差分の正規化レイリー商。小さいほど GT がそのグラフ上で滑らか。
- sel: ガウスノイズ成分のレイリー商 ÷ Rc。大きいほどノイズを分離できる。
- psnr_g: グラフ事前分布だけで行う代理デノイジングの MPSNR。
- sgn: 急峻ノードどうしの辺のうち、符号が合っている割合。

| グラフ (noisy から作成) | JasperRidge Rc / sel / psnr_g / sgn | PaviaU Rc / sel / psnr_g / sgn |
|---|---|---|
| baseline (S=4, 全結合) | 0.129 / 6.43 / 32.79 / 0.65 | 0.467 / 1.77 / 33.57 / 1.00 |
| 案1 (h=0, mkNN k=30) | 0.125 / 5.21 / 30.28 / 1.00 | 0.321 / 2.26 / 33.37 / 1.00 |
| 案2 (h=1, mkNN k=30) | 0.077 / 8.28 / 33.10 / 0.92 | 0.222 / 2.86 / 34.02 / 1.00 |
| 案2 (h=2, mkNN k=30) | 0.061 / 9.10 / 32.89 / 1.00 | 0.196 / 2.89 / 34.03 / 1.00 |
| 案2 (h=2, 符号なし) | 0.075 / 7.49 / 32.16 / 0.91 | 0.208 / 2.87 / 33.89 / 1.00 |
| 案3 (コヒーレンス, k=30) | 0.443 / 1.03 / 30.68 / 0.47 | 0.559 / 0.98 / 31.98 / 0.87 |

参照: 単位行列 (Tikhonov) による psnr_g は、JasperRidge 30.49、PaviaU 32.15。

→ 案2 (窓あり) は、GT の滑らかさ・選択性・代理デノイジングのどれでも現行グラフ以上。案3 は単体でも悪い。

---

## 5. コードとコマンド (どこに何があるか)

| ファイル | 役割 |
|---|---|
| `src/solver.py` | GASSTV 本体 (P-PDS) |
| `src/graphs.py` | 空間グラフ、現行のスペクトルグラフ |
| `src/spectral_graphs.py` | 案1-3、anchored、疎化、符号 |
| `src/data.py` | データ読み込み、`noise_radii` (半径の式) |
| `run_gasstv.py` | GASSTV を 1 回ずつ回す (MATLAB の exp 相当) |
| `run_spectral_graph_comparison.py` | グラフを変えて GASSTV を回す本体。`summary.csv` に 1 実行 1 行 |
| `sweep_all.py` / `sweep_eps_graph.py` / `sweep_noise_conditions.py` | 各段階の sweep (6 節) |
| `verify_spectral_graphs.py` | グラフ単体の検証 |
| `compare_methods.py` | 既存法と GASSTV 系を 1 つの表にする (`--existing_dir` で古い結果も読める) |
| `Matlab/run_existing_methods.m` | 既存法 (SSTV, HSSTV_L1/L12, l0l1HTV, ...) を回す |
| `Matlab/sub_functions/calc_noise_radii.m` | 半径の式 (MATLAB 側) |
| `Matlab/import_existing_results.m` | 過去の MATLAB 実験の結果を `result/old_results/` に取り込む |
| `analysis/noise_conditions/` | ノイズ条件ごとの分析 (スクリプト、report.md、CSV、図) |

結果の置き場所:

- `result/spectral_graph_comparison/<sweep>/summary.csv` (GASSTV 系)
- `result/GASSTV_{result,comp_result}/denoising_<image>/<noise>/<method>/*.mat` (既存法)
- `result/old_results/...` (過去の既存法の結果)

---

## 6. 実験の経緯 (時系列)

| 段階 | 内容 | 条件 | 結論 |
|---|---|---|---|
| ① sweep_all | 現行 / 案1 / 案2 / 案3 で ω, ρ, λ_ρsp, σ_sp と、グラフのパラメータを広く振る (668 実行) | 9 | 案3 は外す。案1 と案2 はほぼ同等で、案2 の方が ω に強い |
| ② sweep_rho_select | ρ ∈ {0.90, 0.95}、提案法 anchored 系を追加 (416 実行) | 9 | **ρ = 0.95 に固定**。anchored 系は現行グラフを超えないので外す |
| ③ 既存法 | MATLAB で ρ ∈ {0.90, 0.95} | 9 | どの手法も ρ = 0.95 の方が良い |
| ④ eps_graph | 案2 の疎化を mutual kNN と eps 近傍で比較、符号 exact / none (212 実行) | 9 | 最良値はどれもほぼ同じ (案2 の 20 設定で幅 J 0.067 / P 0.039 dB)。eps (ε_c=0.5, h=2) が ω に最も頑健 |
| ⑤ noise_conditions | 5 比較対象 × ω ∈ {10, 30, 100, 300, 1000} (252 実行) | 2-7 | 7 節 |
| ⑥ 既存法 | ρ = 0.95、4 手法 | 1-7, 9 | 7 節 |
| ⑦ 過去の結果の取り込み (10-20) | 過去の MATLAB の結果を `result/old_results/` に取り込む (152 件) | 1-9 | 条件 8 の既存法はここにしかない |
| ⑧ 半径の式の変更 (10-20) | デッドラインをスパースの項に含める式に戻す (3.4 節) | 8, 9 に影響 | |
| ⑨ 分析 (10-20) | `analysis/noise_conditions/` | 2-7, 9 | 7 節 |

固定したパラメータ (⑤ 以降): ρ = 0.95, λ_ρsp = 0.8, σ_sp = "90", `lambda2_mode=normalized`, maxiter 20000, stopcri 1e-5。

---

## 7. 実験結果のまとめ

### 7.1 条件 9 での候補選び (①②④、旧式の alpha)

最良の MPSNR [dB] (各案のパラメータを振った中で最良):

| 案 | JasperRidge | PaviaU | 備考 |
|---|---|---|---|
| グラフ項なし (ω=0) | 35.166 | 35.383 | |
| 現行グラフ | 35.396 | 35.692 | ① (PaviaU の最良は k_lap=10, ω=200) |
| 案1 | 35.380 | 35.753 | ① |
| 案2 | 35.396 | 35.745 | ① |
| 案3 | 35.219 | 35.467 | ① 外した |
| anchored (代表ノード + 辺) | 35.384 | 35.472 | ② 外した |
| residual_edges | 35.404 | 35.666 | ② |

- ρ: 0.90 にすると、どの手法も 1.6〜4.2 dB 下がる。
  - 例: 現行グラフは J 33.274 / P 33.976 (ρ = 0.95 では 35.392 / 35.586)。
  - 既存法 (HSSTV_L12) は J 31.136 / P 31.932 (ρ = 0.95 では 35.061 / 35.540)。
- 疎化と符号 (④、最良 / ω=1000 での劣化):

| 設定 | JasperRidge | PaviaU |
|---|---|---|
| 現行グラフ | 35.390 / 1.136 | 35.586 / 0.434 |
| 案2 mkNN k=60 h=1 符号あり (`p2_mknn`) | 35.404 / 1.024 | 35.729 / 0.264 |
| 同 符号なし (`p2_mknn_ns`) | 35.388 / 1.529 | 35.715 / 0.301 |
| 案2 eps ε_c=0.5 h=2 符号あり (`p2_eps`) | 35.343 / 0.495 | 35.704 / 0.162 |

→ ⑤ では `omega0`, `baseline`, `p2_mknn`, `p2_mknn_ns`, `p2_eps` の 5 つを比べることにした。

### 7.2 ノイズ条件ごとの比較 (⑤⑥⑨、詳細は [analysis/noise_conditions/report.md](analysis/noise_conditions/report.md))

最良の MPSNR [dB]。GASSTV 系は ω を振った中の最良、既存法は設定を振った中の最良。太字は列の最良。

JasperRidge

| 手法 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| omega0 | – | 39.515 | 34.742 | 39.092 | 34.771 | **39.290** | 34.877 | – | 35.166 |
| baseline | – | 39.482 | 35.119 | 39.456 | 35.441 | 39.248 | 35.360 | – | 35.390 |
| p2_mknn | – | 39.487 | **35.329** | **39.487** | **35.477** | 39.259 | 35.500 | – | **35.404** |
| p2_eps | – | 39.501 | 35.250 | 39.383 | 35.280 | 39.289 | **35.803** | – | 35.343 |
| 最良の既存法 | **34.673** | **39.537** | 34.786 | 39.369 | 34.936 | 39.240 | 35.062 | 35.197† | 35.061 |

PaviaU

| 手法 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| omega0 | – | 39.340 | 33.271 | 39.408 | 34.996 | 39.354 | 33.597 | – | 35.383 |
| baseline | – | 39.693 | 33.440 | 39.854 | 35.798 | 39.421 | 34.464 | – | 35.586 |
| p2_mknn | – | **39.725** | 33.461 | 39.855 | **35.816** | 39.469 | **34.663** | – | **35.729** |
| p2_eps | – | 39.718 | 33.440 | 39.704 | 35.489 | **39.472** | 34.506 | – | 35.704 |
| 最良の既存法 | **34.876** | 39.342 | **34.014** | **39.871** | 35.208 | 39.457 | 34.196 | 35.375† | 35.540 |

- 最良の既存法は、ほとんどが HSSTV_L12。JasperRidge の 2, 3, 7, 8 だけ l0l1HTV。
- † は過去の実験から取り込んだ結果 (`result/old_results`)。条件 8 の HSSTV は 3 設定中 2 設定だけ。
- 条件 9 は GASSTV 系も既存法も旧式の alpha。
- p2_mknn_ns の行は省いた ([report.md](analysis/noise_conditions/report.md) にある)。

### 7.3 主な数値

- **現行グラフとの差 (最良 ω どうし、14 組 = 7 条件 × 2 画像)**:

| 比較対象 | ΔMPSNR 平均 / 最小 / 最大 | 勝 / 分 / 負 (±0.01) |
|---|---|---|
| omega0 (グラフ項なし) | −0.354 / −0.867 / +0.042 | 2 / 0 / 12 |
| p2_mknn | +0.065 / +0.001 / +0.210 | 12 / 2 / 0 |
| p2_mknn_ns | −0.009 / −0.308 / +0.145 | 6 / 4 / 4 |
| p2_eps | +0.009 / −0.309 / +0.444 | 8 / 1 / 5 |

- **同じ ω どうし (p2_mknn − baseline)**: ω = 10 で平均 +0.011、ω = 100 で +0.057、ω = 1000 で +0.179。ω ≥ 100 では 14 組すべてで p2_mknn が上。
- **ω を全条件で 1 つに固定したときの平均 MPSNR** (14 組):

| 比較対象 | 最も良い共通 ω | 平均 MPSNR |
|---|---|---|
| p2_mknn | 100 | 36.904 |
| p2_eps | 100 | 36.879 |
| baseline | 30 | 36.848 |
| p2_mknn_ns | 30 | 36.835 |
| omega0 | — | 36.629 |

- **符号の効果 (p2_mknn − p2_mknn_ns)**:
  - JasperRidge では ω が大きいほど大きい (ω = 1000 で +0.40〜+1.14)。PaviaU では ω = 1000 でも ≤ +0.055。
  - JasperRidge の σ=0.1 の条件 (3, 7) では最良 ω が大きいので、最良 ω どうしでも +0.33〜+0.45 残る。
- **ノイズの種類**:
  - グラフ項そのものの効果 (baseline − omega0) が大きいのは、ストライプのみ (4, 5) と σ=0.1 のスパース+ストライプ (7) (+0.36〜+0.87)。
  - σ=0.05 でスパースノイズを含む JasperRidge (2, 6) では ≤ 0。
  - 案2 の効果 (p2_mknn − baseline) が 0.05 dB を超えるのは、σ=0.1 でスパースを含む条件 (J3, J7, P7, P9) だけ。
- **既存法との比較**: σ=0.1 の条件 (3, 5, 7, 9) では、GASSTV 系 (p2_mknn) が最良の既存法を +0.19〜+0.61 dB 上回る。例外は PaviaU 条件 3 で、HSSTV_L12 に −0.552 dB 負ける。

### 7.4 結論 (2026-10-27 時点)

- **採用: `p2_mknn`** (案2, `feature=log`, `h=1`, `knn=mutual`, `k=60`, `sign_mode=exact`, `rep=colmed`)。
  - 負けが無い。
  - 同じ ω なら常に現行グラフ以上。
  - ω を固定した比較でも平均が最良。
  - ω を大きくしたときの劣化が小さい。
- **ただし、有意な改善 (> 0.05 dB) は σ=0.1 でスパースノイズを含む条件に限られる。** 他は現行グラフと同等。
- p2_eps は、条件によって大きく良くなる組と悪くなる組があり、ストライプのみの条件では明確に悪い。主案にはしない。
- 符号は残す。符号なしにすると JasperRidge で ω を大きくしたときに大きく劣化する。

---

## 8. 注意点・既知の問題

- **ノイズの実現は 1 通りだけ。** 0.05 dB 未満の差は有意とみなさない。
- **条件 9 の結果 (GASSTV 系・既存法とも) は旧式の alpha**。現行コードで回し直すと値が変わる。
- result/ の既存法の条件 9 (旧式) は、2026-10-20 の時点で「そのまま残す」と決めた。回し直すには `cfg.skip_existing = false` が要り、上書きになる。
- λ_sp を GT から計算している (`--lambda_sp_ref clean`)。元のコードと同じだが、非 Oracle 版でも GT の情報を 1 つ使っている。
- l0l1HTV の `sr0.9999` は全条件で最大反復 (20000) に達している。最良の既存法に選ばれている条件がある。
- 最良 ω が格子の端にある組: J2, J6, J9 (ω=10) と P7 (ω=1000)。
- 過去の実験の l0l1HTV の条件 3 は `stopcri = 1e-4` で、現行 (1e-5) と違う。ファイル名に stopcri が入らないので取り込んでいない。
- `compare_methods.py` の出力をファイルにリダイレクトすると cp932 で落ちる。`PYTHONIOENCODING=utf-8` を付ける。

---

## 9. 次にやること (コマンド付き、未実行)

1. 欠けている条件を埋める (現行の式で回る)。

```
uv run sweep_noise_conditions.py --noise_idx 1 8
matlab -batch "cfg.noise = 8; cfg.rhos = {0.95}; run('Matlab/run_existing_methods.m')"
```

2. 条件 9 を現行の式で回し直す。GASSTV 系は noise_conditions/ に保存される。既存法は上書きになるので、先に決める。

```
uv run sweep_noise_conditions.py --noise_idx 9
matlab -batch "cfg.noise = 9; cfg.rhos = {0.95}; cfg.skip_existing = false; run('Matlab/run_existing_methods.m')"
```

3. ω の格子を広げる。

```
uv run sweep_noise_conditions.py --noise_idx 2 6 9 --images JasperRidge --omega 1 3
uv run sweep_noise_conditions.py --noise_idx 7 --images PaviaU --omega 3000
```

4. 複数のノイズ実現で確認する (Python で生成。既存法の MATLAB 側は今はシードを変えられない)。

```
uv run sweep_noise_conditions.py --noise_idx 3 7 9 --obsv python --seed 1 --result_dir result/spectral_graph_comparison/noise_conditions_seed1
```

5. 集計し直す。

```
uv run analysis/noise_conditions/analyze_noise_conditions.py
```

- 論文向けには、ω を 1 つに固定した比較 (7.3 節) を主にする。条件ごとに ω を選ぶ比較は有利すぎる。

---

## 付録: 用語

| 用語 | 意味 |
|---|---|
| δ (delta) | 代表スペクトルの隣接バンド差分 |
| 急峻なノード | \|δ\| が大きいバンド差分 (吸収帯の縁など) |
| mkNN | mutual kNN。互いに k 近傍のときだけ辺を張る |
| eps 近傍 | 距離が ε_c σ 未満のペアだけ辺を張る |
| sign_mode=exact | ノードの符号 sgn(δ_i) で辺に符号を付ける (平衡な符号付きグラフ) |
| omega0 | ω = 0 (スペクトルグラフ項なし) |
| J3, P7 など | JasperRidge 条件 3, PaviaU 条件 7 |
| 旧式 / 現行の式 | 半径 alpha にデッドラインを含めない式 / 含める式 (3.4 節) |
