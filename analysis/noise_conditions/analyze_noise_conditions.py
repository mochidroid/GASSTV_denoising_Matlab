"""ノイズ条件ごとのスペクトルグラフ比較の集計 (TASK_noise_condition_analysis.md).

既存の結果を読むだけで, 実験は回さない. 出力はすべてこのファイルと同じフォルダに書く:
  completeness.csv, iter_max.csv          : 完全性の確認 (3.0)
  main_table.csv                          : 主表 (3.1)
  delta_best.csv, delta_best_summary.csv  : 最良 omega どうしの現行グラフとの差 (3.2)
  delta_per_omega.csv, delta_per_omega_summary.csv : 同じ omega どうしの差 (3.2)
  robustness.csv, robustness_summary.csv  : omega への頑健性 (3.3)
  sign_effect.csv                         : 符号の効果 (3.4)
  noise_type_effect.csv                   : ノイズの種類・強さとの関係 (3.5)
  fixed_omega_mean.csv, fixed_omega_table.csv : 共通 omega での比較 (3.6)
  vs_existing.csv                         : 最良の既存法との差 (3.7)
  tables.md                               : 上の表を Markdown にしたもの (report.md の元)
  fig_omega_curves.png, fig_fixed_omega.png, fig_delta_per_omega.png

    uv run analysis/noise_conditions/analyze_noise_conditions.py
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))

from compare_methods import collect_existing  # noqa: E402
from src.data import obsv_dirname  # noqa: E402
from sweep_noise_conditions import variant_of  # noqa: E402

# run_gasstv.NOISE_CONDITIONS と同じ (run_gasstv は torch を読み込むので値だけ写す. 下で一致を確認する)
NOISE_CONDITIONS = [
    (0.1, 0, 0, 0, 0), (0.05, 0.05, 0, 0, 0), (0.1, 0.05, 0, 0, 0),
    (0.05, 0, 0.05, 0.5, 0), (0.1, 0, 0.05, 0.5, 0), (0.05, 0.05, 0.05, 0.5, 0),
    (0.1, 0.05, 0.05, 0.5, 0), (0.1, 0, 0, 0, 0.01), (0.1, 0.05, 0.05, 0.5, 0.01),
]
NOISE_TYPE = {1: "ガウスのみ", 2: "+スパース", 3: "+スパース", 4: "+ストライプ", 5: "+ストライプ",
              6: "+両方", 7: "+両方", 8: "+デッドライン", 9: "全部"}

IMAGES = ["JasperRidge", "PaviaU"]
VARIANTS = ["omega0", "baseline", "p2_mknn", "p2_mknn_ns", "p2_eps"]
OMEGAS = [10.0, 30.0, 100.0, 300.0, 1000.0]
EXISTING = ["SSTV", "HSSTV_L1", "HSSTV_L12", "l0l1HTV"]
EXPECTED_SETTINGS = {
    "SSTV": ["r0.95_stop1e-5"],
    "HSSTV_L1": [f"o{o}_r0.95_stop1e-5" for o in ("0.01", "0.03", "0.05")],
    "HSSTV_L12": [f"o{o}_r0.95_stop1e-5" for o in ("0.01", "0.03", "0.05")],
    "l0l1HTV": [f"sr{s}_th{t}_r0.95_maxiter20000" for s in ("0.999", "0.9999") for t in ("0.02", "0.03")],
}
RHO, LRS, SIGMA_SP = 0.95, 0.8, "90"
TIE_DB = 0.01      # |Δ| がこれ以内なら引き分け
GOOD_DB = 0.1      # 良好域: 最良から 0.1 dB 以内
MAXITER = 20000

SWEEP_DIR = ROOT / "result" / "spectral_graph_comparison"
RESULT = ROOT / "result"
OLD = ROOT / "result" / "old_results"

# 図の色: 基準パレットの categorical slot 1-4 を固定順で使う. omega0 は灰色の水平線 (参照線)
COLOR = {"baseline": "#2a78d6", "p2_mknn": "#eb6834", "p2_mknn_ns": "#1baf7a", "p2_eps": "#eda100", "omega0": "#8a8a85"}
MARKER = {"baseline": "o", "p2_mknn": "s", "p2_mknn_ns": "D", "p2_eps": "^"}
LABEL = {"omega0": "omega0 (グラフ項なし)", "baseline": "baseline (現行グラフ)", "p2_mknn": "p2_mknn (案2 mkNN 符号あり)",
         "p2_mknn_ns": "p2_mknn_ns (案2 mkNN 符号なし)", "p2_eps": "p2_eps (案2 eps 符号あり)"}

NOISE_NAME = {i + 1: obsv_dirname(dict(zip(["gaussian_sigma", "sparse_rate", "stripe_rate", "stripe_intensity",
                                            "deadline_rate"], map(float, c))))
              for i, c in enumerate(NOISE_CONDITIONS)}
COND_OF = {v: k for k, v in NOISE_NAME.items()}
unconfirmed: list[str] = []


# ---------------------------------------------------------------------------
# 読み込み
# ---------------------------------------------------------------------------
def check_noise_conditions():
    """run_gasstv.py の NOISE_CONDITIONS とここの写しが一致するかをソースから確認."""
    import ast
    src = (ROOT / "run_gasstv.py").read_text(encoding="utf-8")
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "NOISE_CONDITIONS":
            ref = [tuple(float(x) for x in t) for t in ast.literal_eval(node.value)]
            assert ref == [tuple(float(x) for x in c) for c in NOISE_CONDITIONS], "NOISE_CONDITIONS が違う"
            return
    unconfirmed.append("run_gasstv.py に NOISE_CONDITIONS が見つからない")


def load_gasstv() -> list[dict]:
    """noise_conditions/summary.csv と, 条件 9 は eps_graph/summary.csv から 5 比較対象の行を読む."""
    rows, seen = [], {}
    for sweep, only_cond in (("noise_conditions", None), ("eps_graph", 9)):
        path = SWEEP_DIR / sweep / "summary.csv"
        for r in csv.DictReader(open(path, newline="", encoding="utf-8")):
            cond = COND_OF.get(r["noise"])
            if cond is None:
                unconfirmed.append(f"{sweep}: 条件番号の分からない noise {r['noise']}")
                continue
            if only_cond is not None and cond != only_cond:
                continue
            v = variant_of(r)
            if v is None:
                if sweep == "noise_conditions":
                    unconfirmed.append(f"{sweep}: 比較対象の分からない行 {r['image']} {r['noise']} {r['graph']} "
                                       f"omega={r['omega']} {r['graph_kw_json']}")
                continue
            if (abs(float(r["rho"]) - RHO) > 1e-9 or abs(float(r["lambda_rho_sp"]) - LRS) > 1e-9
                    or r["sigma_sp"] != SIGMA_SP or r.get("guide", "noisy") != "noisy"
                    or r.get("lambda2_mode", "normalized") != "normalized"):
                continue
            key = (r["image"], cond, v, float(r["omega"]))
            if key in seen:
                unconfirmed.append(f"重複 {key}: {seen[key]} と {sweep} (先に読んだ方を使用)")
                continue
            seen[key] = sweep
            rows.append(dict(image=r["image"], cond=cond, noise=r["noise"], variant=v, omega=float(r["omega"]),
                             mpsnr=float(r["mpsnr"]), mssim=float(r["mssim"]), sam=float(r["sam"]),
                             iter=int(float(r["iter"])), source=sweep))
    return rows


def _iteration(full_path: Path):
    """GASSTV_result の .mat から other_result.iteration を読む (無ければ None)."""
    try:
        with h5py.File(full_path, "r") as f:
            return int(np.asarray(f["other_result"]["iteration"][()]).ravel()[0])
    except Exception:  # noqa: BLE001
        return None


def load_existing() -> list[dict]:
    """既存法 (rho=0.95). result/ を優先し, 欠けている設定だけ result/old_results から補う."""
    out, have = [], set()
    for base, source in ((RESULT, "result"), (OLD, "old_results")):
        for r in collect_existing(base, EXISTING):
            cond = COND_OF.get(r["noise"])
            if cond is None:
                continue
            if r["rho"] is None:
                unconfirmed.append(f"rho が読めない: {source}/{r['image']}/{r['noise']}/{r['method']}/{r['setting']}")
                continue
            if abs(r["rho"] - RHO) > 1e-6:
                continue
            key = (r["image"], cond, r["method"], r["setting"])
            if key in have:
                continue
            have.add(key)
            full = base / "GASSTV_result" / f"denoising_{r['image']}" / r["noise"] / r["method"] / f"{r['setting']}.mat"
            out.append(dict(image=r["image"], cond=cond, noise=r["noise"], method=r["method"], setting=r["setting"],
                            mpsnr=r["mpsnr"], mssim=r["mssim"], sam=r["sam"], iter=_iteration(full), source=source))
    return out


# ---------------------------------------------------------------------------
# 表の書式
# ---------------------------------------------------------------------------
def f(x, nd=3, sign=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def md_table(head, rows):
    return ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(map(str, r)) + " |" for r in rows]


def write_csv(name, rows, fields):
    with open(HERE / name, "w", newline="", encoding="utf-8") as fp:
        w = csv.DictWriter(fp, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def bold_best(vals: dict, key, higher=True, nd=3):
    """vals: name -> number or None. 最良を太字にした文字列の dict を返す."""
    ok = {k: v for k, v in vals.items() if v is not None}
    if not ok:
        return {k: "–" for k in vals}
    best = max(ok.values()) if higher else min(ok.values())
    return {k: ("–" if v is None else (f"**{v:.{nd}f}**" if abs(v - best) < 1e-12 else f"{v:.{nd}f}")) for k, v in vals.items()}


# ---------------------------------------------------------------------------
def main():
    check_noise_conditions()
    gas = load_gasstv()
    ext = load_existing()
    md = ["# 集計表 (analyze_noise_conditions.py の出力)", ""]

    # 索引
    G = defaultdict(dict)   # (image, cond, variant) -> omega -> row
    for r in gas:
        G[(r["image"], r["cond"], r["variant"])][r["omega"]] = r
    E = defaultdict(dict)   # (image, cond, method) -> setting -> row
    for r in ext:
        E[(r["image"], r["cond"], r["method"])][r["setting"]] = r

    # ---------------- 3.0 完全性 ----------------
    comp_rows = []
    for image in IMAGES:
        for cond in range(1, 10):
            for v in VARIANTS:
                got = sorted(G.get((image, cond, v), {}))
                exp = [0.0] if v == "omega0" else OMEGAS
                miss = [o for o in exp if o not in got]
                srcs = sorted({G[(image, cond, v)][o]["source"] for o in got}) if got else []
                comp_rows.append(dict(image=image, cond=cond, noise=NOISE_NAME[cond], target=v, kind="GASSTV",
                                      n=len(got), expected=len(exp), missing=" ".join(f"ω={o:g}" for o in miss),
                                      source="+".join(srcs)))
            for m in EXISTING:
                got = E.get((image, cond, m), {})
                miss = [s for s in EXPECTED_SETTINGS[m] if s not in got]
                srcs = sorted({x["source"] for x in got.values()})
                n_old = sum(1 for x in got.values() if x["source"] == "old_results")
                comp_rows.append(dict(image=image, cond=cond, noise=NOISE_NAME[cond], target=m, kind="既存法",
                                      n=len(got), expected=len(EXPECTED_SETTINGS[m]), missing=" ".join(miss),
                                      source="+".join(srcs), n_old=n_old))
    write_csv("completeness.csv", comp_rows,
              ["image", "cond", "noise", "kind", "target", "n", "expected", "missing", "source", "n_old"])

    md += ["## 3.0 完全性", "", "実行数 / 期待数。由来 old は result/old_results から補ったもの。", ""]
    for image in IMAGES:
        head = ["条件"] + VARIANTS + EXISTING
        body = []
        for cond in range(1, 10):
            cells = [f"{cond}: {NOISE_NAME[cond]}"]
            for t in VARIANTS + EXISTING:
                c = next(r for r in comp_rows if r["image"] == image and r["cond"] == cond and r["target"] == t)
                s = f"{c['n']}/{c['expected']}"
                if c.get("n_old"):
                    s += f" (old {c['n_old']})"
                cells.append(s if c["n"] == c["expected"] else f"**{s}**")
            body.append(cells)
        md += [f"### {image}", ""] + md_table(head, body) + [""]
    miss_lines = [f"- {r['image']} 条件{r['cond']} {r['target']}: {r['n']}/{r['expected']}"
                  + (f" (欠け: {r['missing']})" if r["n"] else "") for r in comp_rows if r["n"] < r["expected"]]
    md += ["### 欠けの一覧", ""] + miss_lines + [""]

    it_rows = [dict(kind="GASSTV", image=r["image"], cond=r["cond"], target=r["variant"], setting=f"omega={r['omega']:g}",
                    iter=r["iter"]) for r in gas if r["iter"] >= MAXITER]
    it_rows += [dict(kind="既存法", image=r["image"], cond=r["cond"], target=r["method"], setting=r["setting"],
                     iter=r["iter"]) for r in ext if r["iter"] is not None and r["iter"] >= MAXITER]
    n_noiter = sum(1 for r in ext if r["iter"] is None)
    write_csv("iter_max.csv", it_rows, ["kind", "image", "cond", "target", "setting", "iter"])
    md += [f"### 最大反復 ({MAXITER}) で止まった実行", ""]
    md += [f"- {r['kind']} {r['image']} 条件{r['cond']} {r['target']} {r['setting']}" for r in it_rows] or ["- なし"]
    if n_noiter:
        md += [f"- 既存法 {n_noiter} 件は other_result.iteration を読めず未確認 (LRTDTV などの形式の違い, または GASSTV_result に無い)"]
    md += [""]

    # ---------------- 3.1 主表 ----------------
    def best_g(image, cond, v):
        d = G.get((image, cond, v))
        return max(d.values(), key=lambda r: r["mpsnr"]) if d else None

    def best_e(image, cond, m):
        d = E.get((image, cond, m))
        return max(d.values(), key=lambda r: r["mpsnr"]) if d else None

    main_rows = []
    md += ["## 3.1 主表", "",
           "各比較対象は omega を振った中で MPSNR が最良の実行, 既存法は設定を振った中で MPSNR が最良の実行 (rho=0.95)。"
           "MSSIM, SAM は MPSNR で選んだ実行の値。太字は各列の最良。由来 old は result/old_results から補った結果を含む。", ""]
    for image in IMAGES:
        md += [f"### {image}", ""]
        for cond in range(1, 10):
            ents = {}
            for v in VARIANTS:
                b = best_g(image, cond, v)
                ents[v] = None if b is None else dict(mpsnr=b["mpsnr"], mssim=b["mssim"], sam=b["sam"],
                                                      setting=f"ω={b['omega']:g}", source=b["source"],
                                                      n=len(G[(image, cond, v)]))
            for m in EXISTING:
                b = best_e(image, cond, m)
                ents[m] = None if b is None else dict(mpsnr=b["mpsnr"], mssim=b["mssim"], sam=b["sam"],
                                                      setting=b["setting"],
                                                      source="+".join(sorted({x["source"] for x in E[(image, cond, m)].values()})),
                                                      n=len(E[(image, cond, m)]))
            for name, e in ents.items():
                if e is not None:
                    main_rows.append(dict(image=image, cond=cond, noise=NOISE_NAME[cond], target=name, **e))
            if all(e is None for e in ents.values()):
                md += [f"**条件 {cond}: {NOISE_NAME[cond]}** — 結果なし", ""]
                continue
            bp = bold_best({k: e and e["mpsnr"] for k, e in ents.items()}, "mpsnr", True, 3)
            bs = bold_best({k: e and e["mssim"] for k, e in ents.items()}, "mssim", True, 4)
            bsa = bold_best({k: e and e["sam"] for k, e in ents.items()}, "sam", False, 3)
            body = []
            for k, e in ents.items():
                exp_n = 1 if k == "omega0" else (5 if k in VARIANTS else len(EXPECTED_SETTINGS[k]))
                note = "" if e is None else (e["source"] if e["n"] == exp_n else f"{e['source']}, {e['n']}/{exp_n} 設定")
                body.append([k, bp[k], bs[k], bsa[k], "–" if e is None else e["setting"], note])
            md += [f"**条件 {cond}: {NOISE_NAME[cond]}**", ""] + md_table(["手法", "MPSNR", "MSSIM", "SAM", "選んだ設定", "由来"], body) + [""]
    write_csv("main_table.csv", main_rows,
              ["image", "cond", "noise", "target", "mpsnr", "mssim", "sam", "setting", "source", "n"])

    # 一覧 (MPSNR / SAM, 行 = 手法, 列 = 条件)
    for metric, nd, higher in (("mpsnr", 3, True), ("sam", 3, False)):
        md += [f"### 一覧: {metric.upper()} (列ごとの最良を太字)", ""]
        for image in IMAGES:
            conds = [c for c in range(1, 10) if any(r["image"] == image and r["cond"] == c for r in main_rows)]
            cols = {}
            for c in conds:
                vals = {t: next((r[metric] for r in main_rows if r["image"] == image and r["cond"] == c and r["target"] == t), None)
                        for t in VARIANTS + EXISTING}
                cols[c] = bold_best(vals, metric, higher, nd)
            body = [[t] + [cols[c][t] for c in conds] for t in VARIANTS + EXISTING]
            md += [f"{image}", ""] + md_table(["手法"] + [str(c) for c in conds], body) + [""]

    # ---------------- 3.2 現行グラフとの差 ----------------
    def wtl(ds):
        ds = np.asarray(ds)
        return int((ds > TIE_DB).sum()), int((np.abs(ds) <= TIE_DB).sum()), int((ds < -TIE_DB).sum())

    delta_rows = []
    for image in IMAGES:
        for cond in range(1, 10):
            b = best_g(image, cond, "baseline")
            if b is None:
                continue
            for v in ["omega0", "p2_mknn", "p2_mknn_ns", "p2_eps"]:
                x = best_g(image, cond, v)
                if x is None:
                    continue
                delta_rows.append(dict(image=image, cond=cond, variant=v, d_mpsnr=x["mpsnr"] - b["mpsnr"],
                                       d_sam=x["sam"] - b["sam"], omega=x["omega"], omega_baseline=b["omega"]))
    write_csv("delta_best.csv", delta_rows, ["image", "cond", "variant", "d_mpsnr", "d_sam", "omega", "omega_baseline"])
    dsum = []
    for v in ["omega0", "p2_mknn", "p2_mknn_ns", "p2_eps"]:
        for image in IMAGES + ["両画像"]:
            sub = [r for r in delta_rows if r["variant"] == v and (image == "両画像" or r["image"] == image)]
            if not sub:
                continue
            d = np.array([r["d_mpsnr"] for r in sub])
            s = np.array([r["d_sam"] for r in sub])
            w, t, l = wtl(d)
            dsum.append(dict(variant=v, image=image, n=len(sub), mean=d.mean(), min=d.min(), max=d.max(),
                             win=w, tie=t, loss=l, sam_mean=s.mean(), sam_min=s.min(), sam_max=s.max()))
    write_csv("delta_best_summary.csv", dsum,
              ["variant", "image", "n", "mean", "min", "max", "win", "tie", "loss", "sam_mean", "sam_min", "sam_max"])

    md += ["## 3.2 現行グラフ (baseline) との差", "",
           "Δ = (比較対象の最良) − (baseline の最良)。MPSNR は正, SAM は負が改善。"
           f"勝ち / 引き分け / 負けは MPSNR で |Δ| ≤ {TIE_DB} dB を引き分け。", ""]
    body = []
    for image in IMAGES:
        for cond in sorted({r["cond"] for r in delta_rows if r["image"] == image}):
            cells = [image, cond]
            for v in ["omega0", "p2_mknn", "p2_mknn_ns", "p2_eps"]:
                r = next((r for r in delta_rows if r["image"] == image and r["cond"] == cond and r["variant"] == v), None)
                cells.append("–" if r is None else f"{f(r['d_mpsnr'], 3, True)} / {f(r['d_sam'], 3, True)}")
            body.append(cells)
    md += md_table(["画像", "条件", "omega0 (MPSNR / SAM)", "p2_mknn", "p2_mknn_ns", "p2_eps"], body) + [""]
    md += md_table(["比較対象", "画像", "n", "ΔMPSNR 平均 / 最小 / 最大", "勝 / 分 / 負", "ΔSAM 平均 / 最小 / 最大"],
                   [[r["variant"], r["image"], r["n"], f"{f(r['mean'],3,True)} / {f(r['min'],3,True)} / {f(r['max'],3,True)}",
                     f"{r['win']} / {r['tie']} / {r['loss']}",
                     f"{f(r['sam_mean'],3,True)} / {f(r['sam_min'],3,True)} / {f(r['sam_max'],3,True)}"] for r in dsum]) + [""]

    # 同じ omega どうし
    po_rows = []
    for image in IMAGES:
        for cond in range(1, 10):
            bd = G.get((image, cond, "baseline"), {})
            for v in ["p2_mknn", "p2_mknn_ns", "p2_eps"]:
                xd = G.get((image, cond, v), {})
                for o in OMEGAS:
                    if o in bd and o in xd:
                        po_rows.append(dict(image=image, cond=cond, variant=v, omega=o,
                                            d_mpsnr=xd[o]["mpsnr"] - bd[o]["mpsnr"], d_sam=xd[o]["sam"] - bd[o]["sam"]))
    write_csv("delta_per_omega.csv", po_rows, ["image", "cond", "variant", "omega", "d_mpsnr", "d_sam"])
    pos = []
    for v in ["p2_mknn", "p2_mknn_ns", "p2_eps"]:
        for o in OMEGAS:
            for image in IMAGES + ["両画像"]:
                sub = [r for r in po_rows if r["variant"] == v and r["omega"] == o and (image == "両画像" or r["image"] == image)]
                if not sub:
                    continue
                d = np.array([r["d_mpsnr"] for r in sub])
                s = np.array([r["d_sam"] for r in sub])
                w, t, l = wtl(d)
                pos.append(dict(variant=v, omega=o, image=image, n=len(sub), mean=d.mean(), min=d.min(), max=d.max(),
                                win=w, tie=t, loss=l, sam_mean=s.mean()))
    write_csv("delta_per_omega_summary.csv", pos,
              ["variant", "omega", "image", "n", "mean", "min", "max", "win", "tie", "loss", "sam_mean"])
    md += ["### 同じ omega どうしの差 (比較対象 − baseline)", ""]
    md += md_table(["比較対象", "ω", "画像", "n", "ΔMPSNR 平均 / 最小 / 最大", "勝 / 分 / 負", "ΔSAM 平均"],
                   [[r["variant"], f"{r['omega']:g}", r["image"], r["n"],
                     f"{f(r['mean'],3,True)} / {f(r['min'],3,True)} / {f(r['max'],3,True)}",
                     f"{r['win']} / {r['tie']} / {r['loss']}", f(r["sam_mean"], 3, True)] for r in pos]) + [""]

    # ---------------- 3.3 頑健性 ----------------
    rob = []
    for image in IMAGES:
        for cond in range(1, 10):
            for v in ["baseline", "p2_mknn", "p2_mknn_ns", "p2_eps"]:
                d = G.get((image, cond, v), {})
                if not d:
                    continue
                best = max(r["mpsnr"] for r in d.values())
                at1000 = d[1000.0]["mpsnr"] if 1000.0 in d else np.nan
                rob.append(dict(image=image, cond=cond, variant=v, best=best, at1000=at1000, drop=best - at1000,
                                n_good=sum(1 for r in d.values() if r["mpsnr"] >= best - GOOD_DB), n_omega=len(d),
                                best_omega=max(d.values(), key=lambda r: r["mpsnr"])["omega"]))
    write_csv("robustness.csv", rob, ["image", "cond", "variant", "best", "best_omega", "at1000", "drop", "n_good", "n_omega"])
    rsum = []
    for v in ["baseline", "p2_mknn", "p2_mknn_ns", "p2_eps"]:
        for image in IMAGES + ["両画像"]:
            sub = [r for r in rob if r["variant"] == v and (image == "両画像" or r["image"] == image)]
            if sub:
                rsum.append(dict(variant=v, image=image, n=len(sub), drop_mean=np.nanmean([r["drop"] for r in sub]),
                                 drop_max=np.nanmax([r["drop"] for r in sub]),
                                 n_good_mean=np.mean([r["n_good"] for r in sub])))
    write_csv("robustness_summary.csv", rsum, ["variant", "image", "n", "drop_mean", "drop_max", "n_good_mean"])
    md += ["## 3.3 omega への頑健性", "",
           f"劣化 = 最良 − (ω=1000 の値) [dB]。良好域 = 最良から {GOOD_DB} dB 以内に入る ω の個数 (5 個中)。", ""]
    body = []
    for image in IMAGES:
        for cond in sorted({r["cond"] for r in rob if r["image"] == image}):
            cells = [image, cond]
            for v in ["baseline", "p2_mknn", "p2_mknn_ns", "p2_eps"]:
                r = next((r for r in rob if r["image"] == image and r["cond"] == cond and r["variant"] == v), None)
                cells.append("–" if r is None else f"{f(r['drop'])} / {r['n_good']} (最良 ω={r['best_omega']:g})")
            body.append(cells)
    md += md_table(["画像", "条件", "baseline 劣化 / 良好域", "p2_mknn", "p2_mknn_ns", "p2_eps"], body) + [""]
    md += md_table(["比較対象", "画像", "n", "劣化 平均", "劣化 最大", "良好域 平均"],
                   [[r["variant"], r["image"], r["n"], f(r["drop_mean"]), f(r["drop_max"]), f"{r['n_good_mean']:.2f}"]
                    for r in rsum]) + [""]

    # ---------------- 3.4 符号 ----------------
    sg = []
    for image in IMAGES:
        for cond in range(1, 10):
            a, b = G.get((image, cond, "p2_mknn"), {}), G.get((image, cond, "p2_mknn_ns"), {})
            for o in OMEGAS:
                if o in a and o in b:
                    sg.append(dict(image=image, cond=cond, omega=o, d_mpsnr=a[o]["mpsnr"] - b[o]["mpsnr"],
                                   d_sam=a[o]["sam"] - b[o]["sam"]))
            if a and b:
                ba, bb = max(a.values(), key=lambda r: r["mpsnr"]), max(b.values(), key=lambda r: r["mpsnr"])
                sg.append(dict(image=image, cond=cond, omega="best", d_mpsnr=ba["mpsnr"] - bb["mpsnr"], d_sam=ba["sam"] - bb["sam"]))
    write_csv("sign_effect.csv", sg, ["image", "cond", "omega", "d_mpsnr", "d_sam"])
    md += ["## 3.4 符号の効果 (p2_mknn − p2_mknn_ns)", "", "各セル: ΔMPSNR / ΔSAM。best は各々の最良 ω どうし。", ""]
    body = []
    for image in IMAGES:
        for cond in sorted({r["cond"] for r in sg if r["image"] == image}):
            cells = [image, cond]
            for o in OMEGAS + ["best"]:
                r = next((r for r in sg if r["image"] == image and r["cond"] == cond and r["omega"] == o), None)
                cells.append("–" if r is None else f"{f(r['d_mpsnr'],3,True)} / {f(r['d_sam'],3,True)}")
            body.append(cells)
    md += md_table(["画像", "条件"] + [f"ω={o:g}" for o in OMEGAS] + ["best"], body) + [""]

    # ---------------- 3.5 ノイズの種類・強さ ----------------
    nt = []
    for cond in range(1, 10):
        g, ps, pt, _, pd = NOISE_CONDITIONS[cond - 1]
        for image in IMAGES:
            bb = {v: best_g(image, cond, v) for v in VARIANTS}
            if bb["baseline"] is None or bb["omega0"] is None:
                continue
            o0 = bb["omega0"]["mpsnr"]
            nt.append(dict(cond=cond, type=NOISE_TYPE[cond], sigma=g, image=image,
                           graph_baseline=bb["baseline"]["mpsnr"] - o0,
                           graph_p2_mknn=bb["p2_mknn"]["mpsnr"] - o0 if bb["p2_mknn"] else np.nan,
                           graph_p2_eps=bb["p2_eps"]["mpsnr"] - o0 if bb["p2_eps"] else np.nan,
                           p2_mknn_vs_base=bb["p2_mknn"]["mpsnr"] - bb["baseline"]["mpsnr"] if bb["p2_mknn"] else np.nan,
                           p2_eps_vs_base=bb["p2_eps"]["mpsnr"] - bb["baseline"]["mpsnr"] if bb["p2_eps"] else np.nan))
    write_csv("noise_type_effect.csv", nt, ["cond", "type", "sigma", "image", "graph_baseline", "graph_p2_mknn",
                                            "graph_p2_eps", "p2_mknn_vs_base", "p2_eps_vs_base"])
    md += ["## 3.5 ノイズの種類・強さとの関係", "",
           "グラフ項の効果 = (各グラフの最良) − omega0。案2 の効果 = (案2 の最良) − (baseline の最良)。単位 dB (MPSNR)。", ""]
    md += md_table(["条件", "種類", "σ", "画像", "baseline − omega0", "p2_mknn − omega0", "p2_eps − omega0",
                    "p2_mknn − baseline", "p2_eps − baseline"],
                   [[r["cond"], r["type"], r["sigma"], r["image"], f(r["graph_baseline"], 3, True), f(r["graph_p2_mknn"], 3, True),
                     f(r["graph_p2_eps"], 3, True), f(r["p2_mknn_vs_base"], 3, True), f(r["p2_eps_vs_base"], 3, True)]
                    for r in sorted(nt, key=lambda r: (r["image"], r["cond"]))]) + [""]

    # ---------------- 3.6 共通 omega ----------------
    pairs = sorted({(r["image"], r["cond"]) for r in gas})
    common = [p for p in pairs if all(G.get((p[0], p[1], v)) for v in VARIANTS)
              and all(o in G[(p[0], p[1], v)] for v in VARIANTS[1:] for o in OMEGAS)]
    fo = []
    for v in VARIANTS[1:]:
        for o in OMEGAS:
            vals = [G[(i, c, v)][o] for i, c in common]
            fo.append(dict(variant=v, omega=o, n=len(vals), mpsnr_mean=np.mean([r["mpsnr"] for r in vals]),
                           mssim_mean=np.mean([r["mssim"] for r in vals]), sam_mean=np.mean([r["sam"] for r in vals])))
    o0 = [G[(i, c, "omega0")][0.0] for i, c in common]
    fo.append(dict(variant="omega0", omega=0.0, n=len(o0), mpsnr_mean=np.mean([r["mpsnr"] for r in o0]),
                   mssim_mean=np.mean([r["mssim"] for r in o0]), sam_mean=np.mean([r["sam"] for r in o0])))
    write_csv("fixed_omega_mean.csv", fo, ["variant", "omega", "n", "mpsnr_mean", "mssim_mean", "sam_mean"])
    best_common = {v: max((r for r in fo if r["variant"] == v), key=lambda r: r["mpsnr_mean"])["omega"] for v in VARIANTS[1:]}
    best_common["omega0"] = 0.0
    md += ["## 3.6 共通の omega を使った比較", "",
           f"全条件・両画像で共通の 1 つの ω を使ったときの平均 (対象: {len(common)} 組 = "
           + ", ".join(f"{i[:1]}{c}" for i, c in common) + ")。太字は各比較対象で平均 MPSNR が最良の ω。", ""]
    body = []
    for v in VARIANTS[1:]:
        cells = [v]
        for o in OMEGAS:
            r = next(r for r in fo if r["variant"] == v and r["omega"] == o)
            s = f"{r['mpsnr_mean']:.3f}"
            cells.append(f"**{s}**" if o == best_common[v] else s)
        body.append(cells)
    r0 = next(r for r in fo if r["variant"] == "omega0")
    body.append(["omega0"] + [f"{r0['mpsnr_mean']:.3f}"] * len(OMEGAS))
    md += md_table(["比較対象"] + [f"ω={o:g}" for o in OMEGAS], body) + [""]
    md += md_table(["比較対象", "共通 ω", "平均 MPSNR", "平均 MSSIM", "平均 SAM"],
                   [[v, f"{best_common[v]:g}"] + [f"{x:.{nd}f}" for x, nd in
                    ((r["mpsnr_mean"], 3), (r["mssim_mean"], 4), (r["sam_mean"], 3))]
                    for v in VARIANTS for r in [next(r for r in fo if r["variant"] == v and r["omega"] == best_common[v])]]) + [""]

    ft = []
    for i, c in common:
        vals = {v: G[(i, c, v)][best_common[v]] for v in VARIANTS}
        for v in VARIANTS:
            ft.append(dict(image=i, cond=c, variant=v, omega=best_common[v], mpsnr=vals[v]["mpsnr"],
                           mssim=vals[v]["mssim"], sam=vals[v]["sam"], d_mpsnr=vals[v]["mpsnr"] - vals["baseline"]["mpsnr"]))
    write_csv("fixed_omega_table.csv", ft, ["image", "cond", "variant", "omega", "mpsnr", "mssim", "sam", "d_mpsnr"])
    md += ["### 共通 ω での条件ごとの値 (MPSNR, 括弧は baseline との差。各行の最良を太字)", ""]
    body = []
    for i, c in common:
        vals = {v: next(r for r in ft if r["image"] == i and r["cond"] == c and r["variant"] == v) for v in VARIANTS}
        bb = bold_best({v: vals[v]["mpsnr"] for v in VARIANTS}, "mpsnr")
        body.append([i, c] + [bb[v] + ("" if v in ("baseline", "omega0") else f" ({f(vals[v]['d_mpsnr'],3,True)})") for v in VARIANTS])
    md += md_table(["画像", "条件"] + [f"{v} (ω={best_common[v]:g})" for v in VARIANTS], body) + [""]
    dd = defaultdict(list)
    for r in ft:
        dd[r["variant"]].append(r["d_mpsnr"])
    md += md_table(["比較対象", "ΔMPSNR 平均 / 最小 / 最大 (共通 ω)", "勝 / 分 / 負"],
                   [[v, f"{f(np.mean(dd[v]),3,True)} / {f(np.min(dd[v]),3,True)} / {f(np.max(dd[v]),3,True)}",
                     "{} / {} / {}".format(*wtl(dd[v]))] for v in ["omega0", "p2_mknn", "p2_mknn_ns", "p2_eps"]]) + [""]

    # ---------------- 3.7 既存法との差 ----------------
    vx = []
    for image in IMAGES:
        for cond in range(1, 10):
            be = [(m, best_e(image, cond, m)) for m in EXISTING]
            be = [(m, b) for m, b in be if b is not None]
            if not be:
                continue
            m_best, b_best = max(be, key=lambda x: x[1]["mpsnr"])
            complete = all(len(E.get((image, cond, m), {})) == len(EXPECTED_SETTINGS[m]) for m in EXISTING)
            row = dict(image=image, cond=cond, best_existing=m_best, best_existing_setting=b_best["setting"],
                       best_existing_mpsnr=b_best["mpsnr"], best_existing_sam=b_best["sam"], existing_complete=complete,
                       existing_source="+".join(sorted({x["source"] for m in EXISTING for x in E.get((image, cond, m), {}).values()})))
            for v in VARIANTS:
                b = best_g(image, cond, v)
                row[f"d_{v}"] = b["mpsnr"] - b_best["mpsnr"] if b else np.nan
                if (image, cond) in common:
                    row[f"d_{v}_fixed"] = G[(image, cond, v)][best_common[v]]["mpsnr"] - b_best["mpsnr"]
            vx.append(row)
    write_csv("vs_existing.csv", vx, ["image", "cond", "best_existing", "best_existing_setting", "best_existing_mpsnr",
                                      "best_existing_sam", "existing_complete", "existing_source"]
              + [f"d_{v}" for v in VARIANTS] + [f"d_{v}_fixed" for v in VARIANTS])
    md += ["## 3.7 最良の既存法との差 (MPSNR)", "",
           "左: 各比較対象の最良 ω − 最良の既存法。右 (共通 ω): 3.6 の共通 ω での値 − 最良の既存法。"
           "既存法の設定が揃っていない条件は「既存法 揃い」が no。", ""]
    body = []
    for r in vx:
        body.append([r["image"], r["cond"], f"{r['best_existing']} ({r['best_existing_setting']})", f(r["best_existing_mpsnr"]),
                     "yes" if r["existing_complete"] else "no", r["existing_source"]]
                    + [f(r.get(f"d_{v}"), 3, True) for v in VARIANTS]
                    + [f(r.get(f"d_{v}_fixed"), 3, True) for v in ("baseline", "p2_mknn", "p2_eps")])
    md += md_table(["画像", "条件", "最良の既存法", "MPSNR", "既存法 揃い", "由来"] + [f"Δ {v}" for v in VARIANTS]
                   + ["Δ baseline (共通 ω)", "Δ p2_mknn (共通 ω)", "Δ p2_eps (共通 ω)"], body) + [""]

    # ---------------- 未確認 ----------------
    md += ["## 未確認の点", ""] + ([f"- {u}" for u in unconfirmed] or ["- なし"]) + [""]
    (HERE / "tables.md").write_text("\n".join(md), encoding="utf-8")

    # ---------------- 図 ----------------
    plot_omega_curves(G)
    plot_fixed_omega(fo)
    plot_delta_per_omega(pos)
    print(f"-> {HERE}")


# ---------------------------------------------------------------------------
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for name in ("Yu Gothic", "Meiryo", "MS Gothic", "Noto Sans CJK JP", "IPAexGothic"):
        if any(name == f_.name for f_ in font_manager.fontManager.ttflist):
            plt.rcParams["font.family"] = name
            break
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#b8b7b0",
                         "axes.labelcolor": "#3d3d3a", "xtick.color": "#5f5e5a", "ytick.color": "#5f5e5a",
                         "axes.grid": True, "grid.color": "#e8e7e2", "grid.linewidth": 0.6, "font.size": 9})
    return plt


def plot_omega_curves(G):
    plt = _mpl()
    conds = sorted({c for (i, c, v) in G})
    fig, axes = plt.subplots(len(IMAGES), len(conds), figsize=(2.6 * len(conds), 2.5 * len(IMAGES) + 0.8), squeeze=False)
    for ri, image in enumerate(IMAGES):
        for ci, cond in enumerate(conds):
            ax = axes[ri][ci]
            o0 = G.get((image, cond, "omega0"), {}).get(0.0)
            if o0:
                ax.axhline(o0["mpsnr"], color=COLOR["omega0"], lw=1.5, ls="--", label=LABEL["omega0"])
            for v in ["baseline", "p2_mknn", "p2_mknn_ns", "p2_eps"]:
                d = G.get((image, cond, v), {})
                if not d:
                    continue
                xs = sorted(d)
                ax.plot(xs, [d[x]["mpsnr"] for x in xs], color=COLOR[v], lw=2, marker=MARKER[v], ms=5,
                        markeredgecolor="white", markeredgewidth=1, label=LABEL[v])
            ax.set_xscale("log")
            ax.set_xticks(OMEGAS)
            ax.set_xticklabels([f"{o:g}" for o in OMEGAS])
            ax.minorticks_off()
            ax.set_title(f"{image}  条件{cond}\n{NOISE_NAME[cond]}", fontsize=8.5, color="#3d3d3a")
            if ci == 0:
                ax.set_ylabel("MPSNR [dB]")
            if ri == len(IMAGES) - 1:
                ax.set_xlabel("ω")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=5, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(HERE / "fig_omega_curves.png", dpi=150)
    plt.close(fig)


def plot_fixed_omega(fo):
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    r0 = next(r for r in fo if r["variant"] == "omega0")
    ax.axhline(r0["mpsnr_mean"], color=COLOR["omega0"], lw=1.5, ls="--", label=LABEL["omega0"])
    for v in ["baseline", "p2_mknn", "p2_mknn_ns", "p2_eps"]:
        rs = sorted((r for r in fo if r["variant"] == v), key=lambda r: r["omega"])
        ax.plot([r["omega"] for r in rs], [r["mpsnr_mean"] for r in rs], color=COLOR[v], lw=2, marker=MARKER[v], ms=6,
                markeredgecolor="white", markeredgewidth=1, label=LABEL[v])
    ax.set_xscale("log")
    ax.set_xticks(OMEGAS)
    ax.set_xticklabels([f"{o:g}" for o in OMEGAS])
    ax.minorticks_off()
    ax.set_xlabel("共通の ω")
    ax.set_ylabel("平均 MPSNR [dB]")
    ax.set_title(f"全条件・両画像で共通の ω を使ったときの平均 (n={r0['n']})", fontsize=9, color="#3d3d3a")
    ax.legend(frameon=False, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(HERE / "fig_fixed_omega.png", dpi=150)
    plt.close(fig)


def plot_delta_per_omega(pos):
    plt = _mpl()
    fig, axes = plt.subplots(1, len(IMAGES), figsize=(4.4 * len(IMAGES), 3.2), sharey=True)
    for ax, image in zip(axes, IMAGES):
        ax.axhline(0, color="#8a8a85", lw=1)
        for v in ["p2_mknn", "p2_mknn_ns", "p2_eps"]:
            rs = sorted((r for r in pos if r["variant"] == v and r["image"] == image), key=lambda r: r["omega"])
            xs = [r["omega"] for r in rs]
            ax.plot(xs, [r["mean"] for r in rs], color=COLOR[v], lw=2, marker=MARKER[v], ms=6,
                    markeredgecolor="white", markeredgewidth=1, label=LABEL[v])
            ax.fill_between(xs, [r["min"] for r in rs], [r["max"] for r in rs], color=COLOR[v], alpha=0.12, lw=0)
        ax.set_xscale("log")
        ax.set_xticks(OMEGAS)
        ax.set_xticklabels([f"{o:g}" for o in OMEGAS])
        ax.minorticks_off()
        ax.set_xlabel("ω")
        ax.set_title(f"{image}: 同じ ω での baseline との差 (線: 条件平均, 帯: 最小〜最大)", fontsize=8, color="#3d3d3a")
    axes[0].set_ylabel("ΔMPSNR [dB]")
    axes[0].legend(frameon=False, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(HERE / "fig_delta_per_omega.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
