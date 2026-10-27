# 指示書: ノイズ条件ごとのスペクトルグラフ比較の分析

このリポジトリ (GASSTV: ハイパースペクトル画像の混合ノイズ除去) で, **既に出ている実験結果を読み込んで比較・分析し, レポートを作る**。
新しい実験は回さない (足りない結果は, 回すためのコマンドを示すところまで)。

---

## 1. 背景 (分析に必要な前提)

GASSTV は SSTV + 空間グラフ制約 + **スペクトルグラフ正則化** (GLR) を P-PDS で解く手法。今回はスペクトルグラフの作り方だけを変えて比べている。

- ノードは隣接バンド差分 (K = バンド数 − 1)。正則化は $\frac{\omega}{2}\sum_p (D_\lambda u_p)^\top \tilde L (D_\lambda u_p)$, $\tilde L = L/\lambda_{\max}(L)$ (`lambda2_mode=normalized`)。
- 比較対象 (5 つ):

| 名前 | 内容 | summary.csv での見分け方 |
|---|---|---|
| `omega0` | スペクトルグラフ項なし | `omega == 0` |
| `baseline` | 現行グラフ (セグメントごと, 代表スペクトル差分の値の類似度, 符号なし) | `graph == "baseline"` かつ `omega > 0` |
| `p2_mknn` | 案2 (代表スペクトル差分の絶対値の窓の類似度), mutual kNN k=60, h=1, 符号付き (厳密版) | `graph == "p2_window"`, `graph_kw_json` が knn=mutual, k=60, h=1, sign_mode=exact |
| `p2_mknn_ns` | 上と同じで符号なし | 同上で sign_mode が null (= none) |
| `p2_eps` | 案2, eps 近傍 eps_c=0.5, h=2, 符号付き (厳密版) | knn=eps, eps_c=0.5, h=2, sign_mode=exact |

  見分け方は `sweep_noise_conditions.py` の `variant_of()` をそのまま使ってよい。
- 固定パラメータ: rho = 0.95, lambda_rho_sp = 0.8, sigma_sp = "90"。omega は {10, 30, 100, 300, 1000} を振っている。
- これまでの結果 (条件 9 = g0.1 ps0.05 pt0.05 pd0.01 のみ) で分かっていること (検証対象):
  1. 最良の omega どうしでは, 案2 の各設定の差は 0.06 dB 以内。現行グラフとの差は JasperRidge で ±0, PaviaU で +0.15 dB, SAM も PaviaU で改善。
  2. omega を大きくしたとき (1000) の劣化は, 現行グラフより案2 (特に eps) が小さい = omega への頑健性。
  3. 符号の効果は最良 omega 付近ではほぼ 0。omega が大きいと JasperRidge で +0.3〜0.6 dB (逆符号の急峻な変動がある画像), PaviaU ではほぼ 0。
- ノイズ条件の番号とフォルダ名 (`run_gasstv.NOISE_CONDITIONS`):

| 番号 | フォルダ名 | 内容 |
|---|---|---|
| 1 | g0.1_ps0_pt0_pd0 | ガウス 0.1 |
| 2 | g0.05_ps0.05_pt0_pd0 | ガウス 0.05 + スパース |
| 3 | g0.1_ps0.05_pt0_pd0 | ガウス 0.1 + スパース |
| 4 | g0.05_ps0_pt0.05_pd0 | ガウス 0.05 + ストライプ |
| 5 | g0.1_ps0_pt0.05_pd0 | ガウス 0.1 + ストライプ |
| 6 | g0.05_ps0.05_pt0.05_pd0 | ガウス 0.05 + スパース + ストライプ |
| 7 | g0.1_ps0.05_pt0.05_pd0 | ガウス 0.1 + スパース + ストライプ |
| 8 | g0.1_ps0_pt0_pd0.01 | ガウス 0.1 + デッドライン |
| 9 | g0.1_ps0.05_pt0.05_pd0.01 | 全部 |

## 2. 読むデータ

| 場所 | 中身 | 読み方 |
|---|---|---|
| `result/spectral_graph_comparison/noise_conditions/summary.csv` | 上の 5 比較対象の本体実行 (1 行 1 実行) | csv。列 `image, noise, graph, omega, rho, lambda_rho_sp, graph_kw_json, mpsnr, mssim, sam, iter, best_mpsnr_during, ...` |
| `result/spectral_graph_comparison/eps_graph/summary.csv` | 条件 9 の同じ設定を含む (rho=0.95, lrs=0.8) | 同上。**条件 9 はここから補う** (5 比較対象に該当する行だけ) |
| `result/GASSTV_comp_result/denoising_<image>/<noise>/<method>/*.mat` | 既存法 (SSTV, HSSTV_L1, HSSTV_L12, l0l1HTV) の結果 | MATLAB v7.3 (HDF5)。`compare_methods.py` の `load_matlab_result()` / `collect_existing()` を使う。`params.rho_radius` で rho を確認し, **0.95 だけ使う** |
| `result/old_results/import_log.csv` | 既存法の一部は過去の実験から取り込んだもの (`status`: exists / no_old など) | 取り込んだ結果を使う場合はレポートに由来を明記する |

**2026-10-20 時点で確認した欠け (自分でも確認すること)**:
- `noise_conditions/summary.csv` は条件 2-7 のみ (各画像・各条件 21 行)。条件 1, 8, 9 がない (9 は eps_graph から補える)。
- 既存法: 条件 8 が両画像ともない。PaviaU の条件 7 は SSTV と HSSTV_L1 (1 設定) しかない。

## 3. やること

### 3.0 完全性の確認 (最初に)
- (画像, 条件, 比較対象 / 既存法) ごとの実行数と, 期待される数 (GASSTV 系: omega0 は 1, 他は 5 / 既存法: SSTV 1, HSSTV_L1 3, HSSTV_L12 3, l0l1HTV 4) を表にし, 欠けを列挙する。
- 欠けを埋めるコマンドを示す (実行はしない):
  - GASSTV 系: `uv run sweep_noise_conditions.py --noise_idx 1 8` (など)
  - 既存法: `matlab -batch "cfg.noise = [7 8]; cfg.rhos = {0.95}; run('Matlab/run_existing_methods.m')"` (必要な条件・画像に絞る)
- 以降の分析は, 揃っている条件だけで行い, 欠けている条件は表で「–」とする。
- `iter == 20000` (最大反復で停止) の実行があれば列挙する。

### 3.1 主表: 条件 × 画像 × 手法
- 各 (画像, 条件) で, 各比較対象は omega を振った中の最良 MPSNR, 既存法は各手法のパラメータを振った中の最良 MPSNR (rho=0.95 のみ)。
- MPSNR / MSSIM / SAM を並べ, 各列の最良を太字。最良を選んだ omega (既存法は設定名) も併記。

### 3.2 現行グラフとの対応比較
- 各条件で `Δ = (比較対象の最良) − (baseline の最良)` を MPSNR と SAM について。`omega0` も同様 (グラフ項そのものの効果)。
- 全条件・両画像での平均 / 最小 / 最大, 勝ち / 引き分け / 負け (|Δ| ≤ 0.01 dB を引き分け)。
- 同じ omega どうしの差 (omega ごとの Δ) も出す (最良 omega の選び方に左右されない比較)。

### 3.3 omega への頑健性
- 図: 横軸 omega (対数), 縦軸 MPSNR。画像 × 条件ごとに 5 比較対象の曲線 (omega0 は水平線)。2 画像 × 条件数のグリッド 1 枚 (`fig_omega_curves.png`)。
- 指標: 劣化 = 最良 − (omega=1000 の値), 良好域の幅 = 最良から 0.1 dB 以内に入る omega の個数。条件ごとと平均。

### 3.4 符号の効果
- `p2_mknn − p2_mknn_ns` を omega ごと・条件ごとに (MPSNR, SAM)。
- 「符号の効果は omega が大きいときに, 逆符号の急峻な変動を持つ JasperRidge で表れる」という仮説が全条件で成り立つか。

### 3.5 ノイズの種類・強さとの関係
- 条件をノイズの種類 (ガウスのみ / +スパース / +ストライプ / +両方 / +デッドライン) と強さ (σ = 0.05 / 0.1) で整理し, グラフ項の効果 (omega0 との差) と案2 の効果 (baseline との差) がどう変わるか。

### 3.6 omega を固定したときの比較 (論文向け)
- 各比較対象について, 全条件・両画像で共通の 1 つの omega を使った場合 (5 つの omega それぞれ) の平均 MPSNR を出し, 最も良い共通 omega と, その omega での条件ごとの表を作る (条件ごとに omega を選ぶ比較は有利すぎるため)。

### 3.7 結論
- 採用すべき設定 (p2_mknn / p2_eps / baseline のどれか) と理由。
- これまでの知見 (1 節の 1-3) のうち, 全条件で支持されたもの・されなかったもの。
- 既存法に対する位置づけ (条件ごとに最良の既存法との差)。
- 注意点: ノイズの実現は 1 通り。0.05 dB 未満の差は有意とみなさない。
- 次にやるべき実験 (欠けの補完, 複数シードなど) をコマンド付きで。

## 4. 出力

- `analysis/noise_conditions/analyze_noise_conditions.py`: 上の集計・図をすべて再現するスクリプト (`uv run analysis/noise_conditions/analyze_noise_conditions.py` で動くこと)。
- `analysis/noise_conditions/report.md`: レポート (日本語)。表は Markdown。図はファイル名で参照。
- `analysis/noise_conditions/*.csv`: 主表・差分表など, レポートの表の元データ。
- `analysis/noise_conditions/*.png`: 図 (`fig_omega_curves.png` ほか必要なもの)。

## 5. 守ること

- **既存の結果ファイル・スクリプトは変更しない, 削除しない** (読むだけ)。新しいファイルは `analysis/noise_conditions/` の下に作る。
- 実験 (`run_*.py`, `sweep_*.py` の実行, MATLAB) は回さない。必要ならコマンドを示して止まる。
- 数値の桁: MPSNR 小数 3 桁, MSSIM 4 桁, SAM 3 桁。
- レポートは事実 (数値) と解釈を分けて書く。比喩・誇張・前置きは不要。簡潔に。
- 図を描くのに matplotlib が無ければ `uv add matplotlib`。
- 分からない点 (比較対象の見分けがつかない行, rho が読めない .mat など) は推測で埋めず, レポートに「未確認」として列挙する。
