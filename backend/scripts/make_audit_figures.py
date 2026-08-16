"""Figures for understanding the current results visually.

NOT report figures. These are drawn on the substrate as it stands on 2026-08-16,
before Phase 0.1 drops the four degenerate topology features and re-derives the
components. Regenerate after that lands.

Writes PNG + PDF into audit/figures/.

House palette, from tailwind.config.js: navy/steel/grey for structure, ochre for
the thing being argued about, green for CTCF. No purple.

Run in cd4env.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                   "cd4_cleaned/scripts_cleaned/audit/GW/scripts")

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import hdbscan
import diptest
from _shape import MAG_OVERALL

from app.store import get_store
from app import store_schema as S

warnings.filterwarnings("ignore")

OUT = Path(__file__).resolve().parents[2] / "audit" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

INK, STEEL, PALE = "#2b5070", "#5b89ae", "#c3d6e4"
OCHRE, GREEN, GREY = "#c2703d", "#4a7c59", "#8a99a8"
plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 200, "savefig.bbox": "tight",
})


def save(fig, name, caption):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}")
    plt.close(fig)
    print(f"  {name}.png   {caption}")


# ---------------------------------------------------------------------------
print("loading store and building substrates...")
s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
genes = s.genes.set_index("gene_id").loc[W.index]
sym = genes["symbol_key"].to_numpy()
grp = genes["group"].to_numpy()

rep = s.table("reproducibility_pairs")
rho = {}
for f, sub in rep.groupby("feature"):
    if len(sub) > 10 and sub["gw"].std() > 0 and sub["immune"].std() > 0:
        rho[f] = float(np.corrcoef(sub["gw"].rank(), sub["immune"].rank())[0, 1])

DEGEN = ["mean_degree", "mean_degree_raw", "frac_active_pairs", "n_isolates_raw"]
mag = [m for m in MAG_OVERALL if m in W.columns]
keep = [c for c in W.columns if c not in mag and c not in DEGEN]
A = np.column_stack([np.ones(len(W)), W[mag].to_numpy(float)])
Xc = W[keep].to_numpy(float)
beta, *_ = np.linalg.lstsq(A, Xc, rcond=None)
Rres = Xc - A @ beta
pca_c = PCA(n_components=18, random_state=0).fit(Rres)
D = pca_c.transform(Rres)
D = (D - D.mean(0)) / D.std(0)
n = len(D)

trust_c = []
for j in range(18):
    load = pd.Series(pca_c.components_[j], index=keep)
    m = load[[f for f in keep if f in rho]]
    w = (m ** 2) / (m ** 2).sum()
    trust_c.append(float((w * (pd.Series({f: rho[f] for f in m.index}) > 0.7)).sum()))

cm = s.table("cohort_membership")
sizes = cm.groupby("group").size()
sets = sorted(sizes[sizes >= S.MIN_GROUP_N].index)

print("computing displacements...")
rows = []
for g in sets:
    members = set(cm.loc[cm["group"] == g, "symbol_key"])
    hit = np.fromiter((v in members for v in sym), bool, n)
    k = int(hit.sum())
    if k < S.MIN_GROUP_N:
        continue
    rng = np.random.default_rng(0)
    obs = float(np.linalg.norm(D[hit].mean(0)))
    null = np.array([np.linalg.norm(D[rng.choice(n, size=k, replace=False)].mean(0))
                     for _ in range(1000)])
    j = int(np.abs(D[hit].mean(0)).argmax())
    a_, b_ = D[hit, j], D[~hit, j]
    sp = np.sqrt(((k - 1) * a_.var(ddof=1) + (n - k - 1) * b_.var(ddof=1)) / (n - 2))
    d = float((a_.mean() - b_.mean()) / sp) if sp else 0.0
    from math import erf, sqrt
    rows.append(dict(group=g, n=k, axis=j + 1, var=pca_c.explained_variance_ratio_[j] * 100,
                     trust=trust_c[j] * 100, d=d,
                     p=(np.sum(null >= obs) + 1) / 1001,
                     overlap=100 * (1 - erf(abs(d) / (2 * sqrt(2))))))
ext = pd.DataFrame(rows)
SE = "dbSUPER_CD4_SE_TSS_pm50kb"


# =========================================================================
# FIG 1: the continuum
# =========================================================================
print("fig 1: the continuum...")
fig, ax = plt.subplots(1, 3, figsize=(9.5, 2.9))

Xall = W.to_numpy(float)
Sall = PCA(n_components=19, random_state=0).fit_transform(Xall)
rngp = np.random.default_rng(0)
Xp = Xall.copy()
for j in range(Xp.shape[1]):
    rngp.shuffle(Xp[:, j])
Sp = PCA(n_components=19, random_state=0).fit_transform(Xp)

ks = range(2, 9)
sil = [silhouette_score(Sall, KMeans(k, n_init=10, random_state=0).fit_predict(Sall)) for k in ks]
sil_c = [silhouette_score(D, KMeans(k, n_init=10, random_state=0).fit_predict(D)) for k in ks]
sil_n = [silhouette_score(Sp, KMeans(k, n_init=10, random_state=0).fit_predict(Sp)) for k in ks]
ax[0].plot(list(ks), sil, "o-", color=INK, label="all 91 features", ms=4)
ax[0].plot(list(ks), sil_c, "s-", color=STEEL, label="amount-corrected", ms=4)
ax[0].plot(list(ks), sil_n, "^--", color=GREY, label="permuted null", ms=4)
ax[0].axhline(0.5, color=OCHRE, ls=":", lw=1)
ax[0].text(8, 0.52, "0.5, weak-structure\nconvention", color=OCHRE, fontsize=6, ha="right")
ax[0].set_xlabel("k"); ax[0].set_ylabel("silhouette")
ax[0].set_title("No k gives separated clusters", loc="left")
ax[0].set_ylim(0, 0.62); ax[0].legend(frameon=False, fontsize=6.5)

dips = [diptest.diptest(D[:, i])[1] for i in range(8)]
ax[1].bar(range(1, 9), dips, color=STEEL, width=.65)
ax[1].axhline(0.05, color=OCHRE, ls="--", lw=1)
ax[1].text(1.0, 0.10, "p = 0.05", color=OCHRE, fontsize=6.5)
ax[1].set_ylim(0, 1.05); ax[1].set_xlabel("component"); ax[1].set_ylabel("dip test p")
ax[1].set_title("Every component is unimodal", loc="left")

lab = hdbscan.HDBSCAN(min_cluster_size=25).fit_predict(D)
n_hdb = len(set(lab)) - (1 if -1 in lab else 0)
ax[2].axis("off")
ax[2].text(0.02, 0.80, "HDBSCAN, min_cluster_size 25", fontsize=8, weight="bold")
ax[2].text(0.02, 0.62, f"clusters found:  {n_hdb}", fontsize=11, color=OCHRE, family="monospace")
ax[2].text(0.02, 0.47, f"genes unassigned:  {(lab == -1).mean():.0%}", fontsize=11,
           color=OCHRE, family="monospace")
ax[2].text(0.02, 0.24, "0 clusters in all 16 conditions tested:\n"
                       "4 substrates x whole panel and amount tertiles,\n"
                       "under two definitions of amount.",
           fontsize=7, color="#333", va="top")
ax[2].set_title("Density clustering finds nothing", loc="left")
fig.suptitle("Aim 2: the landscape is continuous, not categorical", x=0.09, ha="left",
             fontsize=10, weight="bold")
save(fig, "fig1_continuum", "no k separates, all components unimodal, HDBSCAN empty")


# =========================================================================
# FIG 2: axis rank, the strongest form of the SE null
# =========================================================================
print("fig 2: axis rank...")
# y is |d|, not -log10 p. Seventeen sets share the permutation floor at
# p = 1/1001, so a p axis stacks them all on one line and no label fits. Effect
# size separates them and says more.
fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.0),
                       gridspec_kw={"width_ratios": [1.5, 1]})
fig.subplots_adjust(top=.82, wspace=.34, bottom=.16)

e = ext.copy()
e["ad"] = e["d"].abs()
is_se = e["group"] == SE
sig = e["p"] < 0.05

sc = ax[0].scatter(e["axis"], e["ad"],
                   s=25 + 55 * np.sqrt(e["n"] / e["n"].max()),
                   c=e["trust"], cmap="Blues", vmin=0, vmax=100,
                   edgecolor=[INK if v else GREY for v in sig],
                   linewidth=[.9 if v else .9 for v in sig],
                   alpha=[1.0 if v else .45 for v in sig], zorder=3)
ax[0].scatter(e.loc[is_se, "axis"], e.loc[is_se, "ad"], s=190, facecolor="none",
              edgecolor=OCHRE, linewidth=2.2, zorder=5)

ax[0].axvspan(0.4, 4.5, color=PALE, alpha=.35, zorder=0)
ax[0].text(2.4, .80, "top 4 components", ha="center", fontsize=7, color=STEEL)

for lab_, dx, dy, ha in [("Lambert_TF", .5, .015, "left"),
                         ("Roadmap_silenced", .5, .015, "left"),
                         ("ChromHMM_bivalent", .5, .01, "left"),
                         ("GWAS_immune_hot", .5, .012, "left"),
                         ("gene_desert_bottomQ_density", .5, -.045, "left")]:
    r = e[e["group"] == lab_]
    if len(r):
        ax[0].annotate(lab_.replace("_bottomQ_density", ""),
                       (r["axis"].iloc[0], r["ad"].iloc[0]),
                       xytext=(r["axis"].iloc[0] + dx, r["ad"].iloc[0] + dy),
                       fontsize=6.8, color=INK, ha=ha, va="center")
r = e[is_se]
ax[0].annotate("super-enhancers\nno preferred direction\nin the measured space",
               (r["axis"].iloc[0], r["ad"].iloc[0]),
               xytext=(r["axis"].iloc[0] + .8, r["ad"].iloc[0] + .20),
               fontsize=7.5, color=OCHRE, weight="bold",
               arrowprops=dict(arrowstyle="->", color=OCHRE, lw=1.4))

ax[0].set_xlabel("component carrying the set's largest displacement (rank)")
ax[0].set_ylabel("|Cohen's d| on that component")
ax[0].set_title("Real structure lives in the top components;\n"
                "super-enhancers lean on the 12th", loc="left", pad=6)
ax[0].set_xlim(0, 18.5); ax[0].set_ylim(0, .88)
ax[0].legend(handles=[
    Line2D([], [], marker="o", ls="", mfc=STEEL, mec=INK, ms=6, label="p < 0.05"),
    Line2D([], [], marker="o", ls="", mfc=STEEL, mec=GREY, ms=6, alpha=.45,
           label="not significant")], frameon=False, fontsize=6.5, loc="upper right")
cb = fig.colorbar(sc, ax=ax[0], pad=.015, fraction=.04)
cb.set_label("% of axis on reproducible features", fontsize=6.2)
cb.ax.tick_params(labelsize=6)

s_, ns_ = e[sig], e[~sig]
bp = ax[1].boxplot([s_["axis"], ns_["axis"]], widths=.5, patch_artist=True,
                   medianprops=dict(color=OCHRE, lw=2))
for patch, c in zip(bp["boxes"], [PALE, "#eeeeee"]):
    patch.set_facecolor(c); patch.set_edgecolor(INK)
for i, (dat, c) in enumerate([(s_, INK), (ns_, GREY)], start=1):
    ax[1].scatter(np.random.default_rng(1).normal(i, .065, len(dat)), dat["axis"],
                  s=16, color=c, zorder=3, alpha=.85)
ax[1].scatter([2], e.loc[is_se, "axis"], s=130, facecolor="none",
              edgecolor=OCHRE, linewidth=2.2, zorder=4)
ax[1].annotate("SE", (2, e.loc[is_se, "axis"].iloc[0]), xytext=(2.28, e.loc[is_se, "axis"].iloc[0]),
               fontsize=7.5, color=OCHRE, weight="bold", va="center")
ax[1].set_xticks([1, 2])
ax[1].set_xticklabels([f"displaced\n(p < 0.05, n={len(s_)})",
                       f"not displaced\n(n={len(ns_)})"], fontsize=7)
ax[1].set_ylabel("axis rank"); ax[1].invert_yaxis(); ax[1].set_xlim(.5, 2.8)
ax[1].set_title(f"median rank {s_['axis'].median():.0f} vs {ns_['axis'].median():.0f}",
                loc="left")
fig.suptitle("Aim 3: where a set leans says more than how far", x=0.02, y=.99,
             ha="left", fontsize=10.5, weight="bold")
save(fig, "fig2_axis_rank", "SE leans on sPC12; displaced sets lean on the top 4")


# =========================================================================
# FIG 3: displaced but not separated
# =========================================================================
print("fig 3: displacement forest...")
fig, ax = plt.subplots(1, 2, figsize=(9.5, 4.4),
                       gridspec_kw={"width_ratios": [1, 1]})
o = ext.sort_values("d")
cols = [OCHRE if g == SE else (GREEN if g == "gene_desert_bottomQ_density" else INK)
        for g in o["group"]]
ax[0].barh(range(len(o)), o["d"], color=cols, height=.68)
ax[0].axvline(0, color="#666", lw=.8)
for x, lab_ in [(-0.8, "large"), (-0.5, "medium"), (-0.2, "small")]:
    ax[0].axvline(x, color=GREY, ls=":", lw=.7)
    ax[0].axvline(-x, color=GREY, ls=":", lw=.7)
ax[0].set_yticks(range(len(o)))
ax[0].set_yticklabels([g.replace("_", " ")[:30] for g in o["group"]], fontsize=6.5)
ax[0].set_xlabel("Cohen's d on the set's strongest component")
ax[0].set_title("Displaced, but never separated", loc="left")
ax[0].set_xlim(-1.0, 1.0)

ax[1].scatter(o["d"].abs(), o["overlap"], s=20 + 45 * np.sqrt(o["n"] / o["n"].max()),
              c=cols, zorder=3)
for g in [SE, "Lambert_TF", "Roadmap_silenced", "GWAS_immune_hot"]:
    r = o[o["group"] == g]
    if len(r):
        ax[1].annotate(g.replace("_", " ")[:22],
                       (abs(r["d"].iloc[0]), r["overlap"].iloc[0]),
                       xytext=(abs(r["d"].iloc[0]) + .025, r["overlap"].iloc[0] + 1.6),
                       fontsize=6.5, color=OCHRE if g == SE else INK,
                       weight="bold" if g == SE else "normal")
ax[1].axhline(70, color=GREY, ls=":", lw=.8)
ax[1].text(.02, 71, "70% overlap", fontsize=6, color=GREY)
ax[1].set_xlabel("|d|"); ax[1].set_ylabel("% overlap of the two distributions")
ax[1].set_ylim(65, 100)
ax[1].set_title("Even the largest effect leaves 71% overlap", loc="left")
fig.suptitle("External categories are directions, not regions", x=0.06, ha="left",
             fontsize=10, weight="bold")
save(fig, "fig3_displaced_not_separated", "d and overlap, SE ochre, density control green")


# =========================================================================
# FIG 4: the dimensions and how much to trust them
# =========================================================================
print("fig 4: dimensions...")
scree = s.table("pc_scree")
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.0))
r20 = scree[scree["pc"] <= 20]
c = [INK if a else PALE for a in r20["above_noise"]]
ax[0].bar(r20["pc"], r20["variance_pct"], color=c, width=.7)
ax[0].plot(r20["pc"], r20["noise_pct"], color=OCHRE, lw=1.3, label="permuted noise ceiling")
ax[0].set_xlabel("component"); ax[0].set_ylabel("% variance")
ax[0].legend(frameon=False, fontsize=6.5)
ax[0].set_title(f"{int(scree['above_noise'].sum())} components above the noise ceiling", loc="left")

tr = np.array(trust_c) * 100
bars = ax[1].bar(range(1, 19), tr, color=[INK if t >= 50 else OCHRE for t in tr], width=.7)
ax[1].axhline(50, color=GREY, ls="--", lw=.9)
ax[1].set_xlabel("component (amount-corrected substrate)")
ax[1].set_ylabel("% of loading on features with rho > 0.7")
ax[1].set_title("Not all components are equally measurable", loc="left")
ax[1].set_ylim(0, 100)
ax[1].legend(handles=[Line2D([], [], color=INK, lw=6, label="trustworthy"),
                      Line2D([], [], color=OCHRE, lw=6, label="below 50%")],
             frameon=False, fontsize=6.5, loc="lower right")
fig.suptitle("There is a lot of real structure, and it is not all equally solid",
             x=0.06, ha="left", fontsize=10, weight="bold")
save(fig, "fig4_dimensions", "scree with noise ceiling, plus per-component trust")


# =========================================================================
# FIG 5: Aim 1, nested baselines
# =========================================================================
print("fig 5: nested baselines...")
nb = pd.read_csv("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/"
                 "scripts_cleaned/audit/continuous_methods/nested_baselines.tsv", sep="\t")
lin = nb[nb["model"] == "linear"]
piv = lin.pivot_table(index="target", columns="features", values="score")
order = piv["4. all features (91)"].sub(piv["3. magnitude (11)"]).sort_values().index
fig, ax = plt.subplots(figsize=(7.6, 3.6))
y = np.arange(len(order))
ax.barh(y - .21, piv.loc[order, "3. magnitude (11)"], height=.4, color=PALE,
        label="11-feature magnitude basis")
ax.barh(y + .21, piv.loc[order, "4. all features (91)"], height=.4, color=INK,
        label="all 91 features")
for i, t in enumerate(order):
    gap = piv.loc[t, "4. all features (91)"] - piv.loc[t, "3. magnitude (11)"]
    ax.text(piv.loc[t, "4. all features (91)"] + .012, i,
            f"+{gap:.3f}", va="center", fontsize=6.5,
            color=OCHRE if gap > .04 else GREY,
            weight="bold" if gap > .04 else "normal")
ax.set_yticks(y); ax.set_yticklabels(order, fontsize=7)
ax.set_xlabel("cross-validated score (AUC for binary, r for continuous)")
ax.legend(frameon=False, fontsize=7, loc="lower right")
ax.set_title("Aim 1: shape adds to amount on 9 of 9 targets, but substantially on two",
             loc="left", fontsize=9, weight="bold")
ax.set_xlim(0, .85)
save(fig, "fig5_aim1_baselines", "91 features vs the magnitude basis, margins annotated")


# =========================================================================
# FIG 6: the map itself, with the k=2 cut and arch-HK
# =========================================================================
print("fig 6: the map...")
emb = s.table("embeddings").set_index("gene_id").loc[W.index]
lab2 = KMeans(2, n_init=20, random_state=0).fit_predict(D)
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.3))
GC = {"arch-HK": "#3d6b91", "arch-ME-constitutive": INK, "arch-ME-effector": OCHRE,
      "arch-sparse": "#8badc9", "arch-off": PALE}
for a, (col, title) in zip(ax, [(None, "The panel, PC2 x PC3"),
                                ("group", "Imposed k=4 partition"),
                                ("k2", "The k=2 cut, amount-free")]):
    x, y_ = emb["pc2"].to_numpy(), emb["pc3"].to_numpy()
    if col is None:
        a.scatter(x, y_, s=3, color=STEEL, alpha=.5, lw=0)
    elif col == "group":
        for g_, c_ in GC.items():
            m = grp == g_
            a.scatter(x[m], y_[m], s=3, color=c_, alpha=.75, lw=0, label=g_.replace("arch-", ""))
        a.legend(frameon=False, fontsize=5.5, markerscale=2.2, loc="upper left")
    else:
        # Cluster index is arbitrary, so name them by their entropy, which is
        # what the split is actually made of.
        ent = W["signal_entropy"].to_numpy()
        hi = 1 if ent[lab2 == 1].mean() > ent[lab2 == 0].mean() else 0
        for v in (1 - hi, hi):
            m = lab2 == v
            l_ = "dispersed / long-range" if v == hi else "focal / local"
            a.scatter(x[m], y_[m], s=3, color=(INK if v == hi else "#8badc9"),
                      alpha=.75, lw=0, label=l_)
        a.legend(frameon=False, fontsize=6, markerscale=2.2, loc="upper left")
    a.set_title(title, loc="left"); a.set_xticks([]); a.set_yticks([])
    a.set_xlabel("PC2"); a.set_ylabel("PC3")
v = 0 if ((grp == "arch-HK") & (lab2 == 0)).sum() > ((grp == "arch-HK") & (lab2 == 1)).sum() else 1
agree = lab2 == v
n_hk = int((grp == "arch-HK").sum())
ov = int(((grp == "arch-HK") & agree).sum())
exp = n_hk * agree.sum() / len(grp)
fig.subplots_adjust(top=.74)
fig.suptitle(f"On the correct substrate the k=2 cut is NOT arch-HK: {ov} of {n_hk} "
             f"arch-HK genes against {exp:.0f} expected by chance ({ov / exp:.2f}x)",
             x=0.02, y=1.02, ha="left", fontsize=10, weight="bold")
fig.text(0.02, .90, "MAG_OVERALL-corrected, degenerate topology features dropped. An earlier\n"
                    "correction against total_mcc alone gave 819 of 844; that number is retracted.",
         fontsize=6.8, color=GREY, ha="left", va="top")
save(fig, "fig6_map_and_k2", "the map, the imposed partition, and the natural cut")

print(f"\nwrote {len(list(OUT.glob('*.png')))} figures to {OUT}")
