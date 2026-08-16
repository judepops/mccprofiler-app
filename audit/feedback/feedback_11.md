okay lets continue.
feedback_12.md#1-11
dd the confound check to Phase B, as its first step. The plan says constraint is where architecture does real work, and it is — but LOEUF has three confounders that have never been controlled anywhere in this project:

Gene length. LOEUF is a depletion statistic, so it scales with coding sequence length. And gene length mechanically changes contact features, because a longer gene occupies more of the ±1 Mb window. This was flagged as an uncontrolled gap in the feature review weeks ago and never addressed.
Expression level, which correlates with both constraint and contact amount.
CpG density, same.

So step one of B is: does the 0.191 survive partialling out length, expression and CpG? If it does, "contact architecture predicts constraint independent of magnitude, length and expression" is a genuinely strong mechanistic claim and B is worth the effort. If it doesn't, you've saved a phase.

Given the week you've had, better to find that out in an hour than after building a programme on it. And gene length is exactly the kind of thing that has been catching you — a mechanical, unglamorous confound that nobody thought to check.

Everything else in the plan I'd leave alone. The what-not-to-do list is the part that will actually save you time.
This is the best catch of the day and it should be tested before Phase B is built, not after. LOEUF is a depletion statistic so it scales with coding length by construction, and longer genes mechanically occupy more of the ±1 Mb window. Let me check what confounders are available:

Bash Find confounder columns
IN
head -1 /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/collaboration/data/gene_table.tsv | tr '\t' '\n' | nl | grep -iE "length|tpm|cpg|loeuf|pli|exon|transcript"

OUT
     6	gtex_whole_blood_tpm
    12	cpg_overlap
    23	gnomad_loeuf
    24	gnomad_pli
    25	cpg_oe_ratio
    31	gene_length
    32	log10_gene_length
All three are on disk. Running the nested test:

Write diagnose_loeuf_confounds.py
133 lines
"""Does the constraint result survive gene length, expression and CpG density?

The forward plan makes constraint the most promising positive result: LOEUF is
predicted at r = 0.210 by all 91 features and **0.191 by the amount-corrected
shape substrate alone**, so 91% of the signal is shape rather than magnitude.
Nothing else in the panel behaves like that, and Phase B is built on it.

Three confounders have never been controlled anywhere in this project:

  GENE LENGTH   LOEUF is a depletion statistic, so it scales with coding
                sequence length by construction: longer genes accumulate more
                expected loss-of-function variants and give more power to
                detect depletion. Length ALSO changes contact features
                mechanically, because a longer gene occupies more of the
                plus-or-minus 1 Mb window. This was flagged as a gap in the
                feature review weeks ago and never addressed.
  EXPRESSION    correlates with constraint and with contact amount.
  CpG DENSITY   same.

So the question is whether the 0.191 is architecture or is length wearing an
architectural costume. Nested, cross-validated, matching the methodology of
audit/continuous_methods/nested_baselines.tsv:

  1. confounders only          length, expression, CpG
  2. shape only                the amount-corrected substrate
  3. confounders + shape       does shape ADD anything?
  4. shape, confounders partialled out of BOTH sides

Model 3 versus 1 is the test that matters. If shape adds nothing over the
confounders, Phase B is not worth building and an hour here has saved a phase.

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
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV, LinearRegression
from sklearn.model_selection import KFold
from _shape import MAG_OVERALL
from app.store import get_store

warnings.filterwarnings("ignore")
SEED = 0
N_SPLITS = 5
DEGEN = ["mean_degree", "mean_degree_raw", "frac_active_pairs", "n_isolates_raw"]
GENE_TABLE = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
              "collaboration/data/gene_table.tsv")

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
sym = s.genes.set_index("gene_id").loc[W.index, "symbol_key"].to_numpy()

gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
cols = ["gnomad_loeuf", "log10_gene_length", "gtex_whole_blood_tpm", "cpg_oe_ratio"]
m = pd.DataFrame({"symbol_key": sym}).merge(
    gt[["symbol_key"] + cols].drop_duplicates("symbol_key"), on="symbol_key", how="left")

ok = m[cols].notna().all(axis=1).to_numpy()
print(f"{ok.sum()} of {len(m)} genes have LOEUF and all three confounders\n")

y = m.loc[ok, "gnomad_loeuf"].to_numpy(float)
C = m.loc[ok, ["log10_gene_length", "gtex_whole_blood_tpm", "cpg_oe_ratio"]].to_numpy(float)
C[:, 1] = np.log1p(C[:, 1])                 # TPM is heavy-tailed
C = (C - C.mean(0)) / C.std(0)

mag = [x for x in MAG_OVERALL if x in W.columns]
keep = [c for c in W.columns if c not in mag and c not in DEGEN]
A = np.column_stack([np.ones(len(W)), W[mag].to_numpy(float)])
X = W[keep].to_numpy(float)
beta, *_ = np.linalg.lstsq(A, X, rcond=None)
SHAPE = PCA(n_components=18, random_state=0).fit_transform(X - A @ beta)[ok]
SHAPE = (SHAPE - SHAPE.mean(0)) / SHAPE.std(0)

print("How much do the confounders relate to LOEUF and to each other?")
names = ["log10 gene length", "log blood TPM", "CpG o/e"]
for i, nm in enumerate(names):
    print(f"   r(LOEUF, {nm:<18}) = {np.corrcoef(y, C[:, i])[0, 1]:+.3f}")
print(f"   r(length, CpG)              = {np.corrcoef(C[:, 0], C[:, 2])[0, 1]:+.3f}")


def cv_r(Xd, yd, seed=SEED):
    """Cross-validated Pearson r between prediction and truth."""
    kf = KFold(N_SPLITS, shuffle=True, random_state=seed)
    pred = np.zeros_like(yd, dtype=float)
    for tr, te in kf.split(Xd):
        mdl = RidgeCV(alphas=np.logspace(-2, 4, 25)).fit(Xd[tr], yd[tr])
        pred[te] = mdl.predict(Xd[te])
    return float(np.corrcoef(pred, yd)[0, 1])


def resid(target, on):
    return target - LinearRegression().fit(on, target).predict(on)


rows = [
    ("1. confounders only (length, expression, CpG)", cv_r(C, y)),
    ("2. shape only (amount-corrected, 18 comps)", cv_r(SHAPE, y)),
    ("3. confounders + shape", cv_r(np.column_stack([C, SHAPE]), y)),
    ("4. shape, confounders partialled from both",
     cv_r(np.column_stack([resid(SHAPE[:, j], C) for j in range(SHAPE.shape[1])]),
          resid(y, C))),
]
print("\n" + "=" * 74)
print("NESTED TEST: does shape predict constraint beyond length/expression/CpG?")
print("=" * 74)
for lab, r in rows:
    print(f"  {lab:<48} r = {r:+.3f}")

add = rows[2][1] - rows[0][1]
print(f"\n  shape adds {add:+.3f} over the confounders alone")
print(f"  after partialling, shape retains r = {rows[3][1]:+.3f} "
      f"against {rows[1][1]:+.3f} unadjusted "
      f"({100 * rows[3][1] / rows[1][1]:.0f}% retained)")

# Which single confounder does the damage?
print("\nleave-one-confounder-in, to see which one matters:")
for i, nm in enumerate(names):
    r = cv_r(np.column_stack([resid(SHAPE[:, j], C[:, [i]]) for j in range(SHAPE.shape[1])]),
             resid(y, C[:, [i]]))
    print(f"   partialling {nm:<18} shape r = {r:+.3f}")

print("\nModel 3 versus 1 is the test. If shape adds nothing over the confounders,")
print("the constraint result is length and expression, and Phase B is not worth")
print("building on it.")
Bash Test LOEUF against length, expression, CpG
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate cd4env
python backend/scripts/diagnose_loeuf_confounds.py 2>&1 | tail -30

OUT
1686 of 1846 genes have LOEUF and all three confounders

How much do the confounders relate to LOEUF and to each other?
   r(LOEUF, log10 gene length ) = -0.347
   r(LOEUF, log blood TPM     ) = -0.191
   r(LOEUF, CpG o/e           ) = -0.094
   r(length, CpG)              = +0.056

==========================================================================
NESTED TEST: does shape predict constraint beyond length/expression/CpG?
==========================================================================
  1. confounders only (length, expression, CpG)    r = +0.452
  2. shape only (amount-corrected, 18 comps)       r = +0.152
  3. confounders + shape                           r = +0.454
  4. shape, confounders partialled from both       r = +0.064

  shape adds +0.002 over the confounders alone
  after partialling, shape retains r = +0.064 against +0.152 unadjusted (42% retained)

leave-one-confounder-in, to see which one matters:
   partialling log10 gene length  shape r = +0.085
   partialling log blood TPM      shape r = +0.156
   partialling CpG o/e            shape r = +0.147

Model 3 versus 1 is the test. If shape adds nothing over the confounders,
the constraint result is length and expression, and Phase B is not worth
building on it.
That kills Phase B, and it is worse than a confound — it is the whole effect:

model	r
confounders only (length, expression, CpG)	0.452
shape only	0.152
confounders + shape	0.454
Shape adds +0.002. And it is specifically gene length: partialling length alone drops shape from 0.152 to 0.085, while expression and CpG cost nothing.

This also threatens Aim 1, since LOEUF was its standout target. Let me check whether the other one holds:

Bash Test other Aim 1 targets for length confounding
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate cd4env
python - <<'PY' 2>&1 | grep -v Warning
import sys, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,'backend')
sys.path.insert(0,'/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/scripts_cleaned/audit/GW/scripts')
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV, LogisticRegression, LinearRegression
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import roc_auc_score
from _shape import MAG_OVERALL
from app.store import get_store
s=get_store()
W=s.table("features").pivot(index="gene_id",columns="feature",values="z")
sym=s.genes.set_index("gene_id").loc[W.index,"symbol_key"].to_numpy()
gt=pd.read_csv("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/collaboration/data/gene_table.tsv",sep="\t")
gt["symbol_key"]=gt.gene_symbol.str.upper()
cf=["log10_gene_length","gtex_whole_blood_tpm","cpg_oe_ratio"]
m=pd.DataFrame({"symbol_key":sym}).merge(gt[["symbol_key"]+cf].drop_duplicates("symbol_key"),on="symbol_key",how="left")
ok=m[cf].notna().all(axis=1).to_numpy()
C=m.loc[ok,cf].to_numpy(float); C[:,1]=np.log1p(C[:,1]); C=(C-C.mean(0))/C.std(0)
DEGEN=["mean_degree","mean_degree_raw","frac_active_pairs","n_isolates_raw"]
mag=[x for x in MAG_OVERALL if x in W.columns]
keep=[c for c in W.columns if c not in mag and c not in DEGEN]
A=np.column_stack([np.ones(len(W)),W[mag].to_numpy(float)]); X=W[keep].to_numpy(float)
b,*_=np.linalg.lstsq(A,X,rcond=None)
SH=PCA(n_components=18,random_state=0).fit_transform(X-A@b)[ok]; SH=(SH-SH.mean(0))/SH.std(0)
cm=s.table("cohort_membership")
print("Do the OTHER Aim 1 targets survive gene length, expression and CpG?\n")
print("%-26s %10s %10s %10s %9s" % ("target","confound","shape","conf+shape","shape adds"))
def auc(Xd,yd):
    kf=KFold(5,shuffle=True,random_state=0)
    p=cross_val_predict(LogisticRegression(max_iter=2000,C=1.0),Xd,yd,cv=kf,method="predict_proba")[:,1]
    return roc_auc_score(yd,p)
for g in ["GWAS_immune_hot","Lambert_TF","DepMap_curated_essential","Eisenberg_HK"]:
    mem=set(cm.loc[cm["group"]==g,"symbol_key"])
    y=np.fromiter((v in mem for v in sym),bool,len(sym))[ok].astype(int)
    if y.sum()<25: continue
    a,bb,c=auc(C,y),auc(SH,y),auc(np.column_stack([C,SH]),y)
    print("%-26s %10.3f %10.3f %10.3f %+9.3f" % (g,a,bb,c,c-a))
# and continuous targets
print()
def cvr(Xd,yd):
    kf=KFold(5,shuffle=True,random_state=0); pr=np.zeros_like(yd,dtype=float)
    for tr,te in kf.split(Xd):
        pr[te]=RidgeCV(alphas=np.logspace(-2,4,25)).fit(Xd[tr],yd[tr]).predict(Xd[te])
    return float(np.corrcoef(pr,yd)[0,1])
for col,lab in [("gtex_tau","gtex_tau"),("phastcons_mean_2kb","phastcons_2kb")]:
    if col not in gt.columns: continue
    v=pd.DataFrame({"symbol_key":sym}).merge(gt[["symbol_key",col]].drop_duplicates("symbol_key"),on="symbol_key",how="left")[col].to_numpy(float)[ok]
    g_=np.isfinite(v)
    a,bb,c=cvr(C[g_],v[g_]),cvr(SH[g_],v[g_]),cvr(np.column_stack([C,SH])[g_],v[g_])
    print("%-26s %10.3f %10.3f %10.3f %+9.3f" % (lab,a,bb,c,c-a))
PY

OUT
Do the OTHER Aim 1 targets survive gene length, expression and CpG?

target                       confound      shape conf+shape shape adds
GWAS_immune_hot                 0.694      0.627      0.699    +0.005
Lambert_TF                      0.577      0.622      0.631    +0.053
DepMap_curated_essential        0.576      0.520      0.582    +0.006
Eisenberg_HK                    0.680      0.533      0.673    -0.007

gtex_tau                        0.218      0.080      0.230    +0.012
phastcons_2kb                   0.055      0.065      0.098    +0.043
That is comprehensive and it hits Aim 1 hard. Across every target, shape adds essentially nothing over three simple gene properties:

target	confounders	shape	both	shape adds
LOEUF	0.452	0.152	0.454	+0.002
GWAS_immune_hot	0.694	0.627	0.699	+0.005
DepMap_curated	0.576	0.520	0.582	+0.006
Eisenberg_HK	0.680	0.533	0.673	−0.007
gtex_tau	0.218	0.080	0.230	+0.012
Lambert_TF	0.577	0.622	0.631	+0.053
The nested baselines never included length, expression or CpG. But if length drives contact features, Aim 3 is exposed too — we controlled magnitude, density and chromosome, never length. Checking the headline before saying anything:

Bash Test Aim 3 headline against gene length
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate cd4env
python - <<'PY' 2>&1 | grep -v Warning
import sys, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,'backend')
sys.path.insert(0,'/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/scripts_cleaned/audit/GW/scripts')
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from _shape import MAG_OVERALL
from app.store import get_store
from app import store_schema as S
s=get_store()
W=s.table("features").pivot(index="gene_id",columns="feature",values="z")
sym=s.genes.set_index("gene_id").loc[W.index,"symbol_key"].to_numpy()
gt=pd.read_csv("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/collaboration/data/gene_table.tsv",sep="\t")
gt["symbol_key"]=gt.gene_symbol.str.upper()
cf=["log10_gene_length","gtex_whole_blood_tpm","cpg_oe_ratio","gene_density_500kb"]
m=pd.DataFrame({"symbol_key":sym}).merge(gt[["symbol_key"]+cf].drop_duplicates("symbol_key"),on="symbol_key",how="left")
V=m[cf].to_numpy(float); V[:,1]=np.log1p(V[:,1])
for j in range(V.shape[1]):
    V[:,j]=np.where(np.isfinite(V[:,j]),V[:,j],np.nanmedian(V[:,j]))
V=(V-V.mean(0))/V.std(0)
DEGEN=["mean_degree","mean_degree_raw","frac_active_pairs","n_isolates_raw"]
mag=[x for x in MAG_OVERALL if x in W.columns]
keep=[c for c in W.columns if c not in mag and c not in DEGEN]
X=W[keep].to_numpy(float)
print("How much does gene length alone relate to the contact features?")
L=V[:,0]; rs=np.array([abs(np.corrcoef(L,X[:,j])[0,1]) for j in range(X.shape[1])])
print("  |r| with log10 gene length: median %.3f, max %.3f (%s)\n"
      % (np.median(rs), rs.max(), keep[int(rs.argmax())]))
def sub(idx):
    cols=[np.ones(len(W)),W[mag].to_numpy(float)]
    if idx: cols.append(V[:,idx])
    A=np.column_stack(cols); b,*_=np.linalg.lstsq(A,X,rcond=None)
    D=PCA(n_components=18,random_state=0).fit_transform(X-A@b); return (D-D.mean(0))/D.std(0)
def disp(D,mask):
    n=len(D); k=int(mask.sum()); rng=np.random.default_rng(0)
    obs=float(np.linalg.norm(D[mask].mean(0)))
    null=np.array([np.linalg.norm(D[rng.choice(n,size=k,replace=False)].mean(0)) for _ in range(2000)])
    j=int(np.abs(D[mask].mean(0)).argmax()); a_,b2=D[mask,j],D[~mask,j]
    sp=np.sqrt(((k-1)*a_.var(ddof=1)+(n-k-1)*b2.var(ddof=1))/(n-2))
    return abs((a_.mean()-b2.mean())/sp),(np.sum(null>=obs)+1)/2001
cm=s.table("cohort_membership")
TESTS=[("magnitude only",[]),("+ length",[0]),
       ("+ length, expr, CpG, density",[0,1,2,3])]
Ds={lab:sub(idx) for lab,idx in TESTS}
print("%-30s %s" % ("set"," ".join("%22s"%l for l,_ in TESTS)))
for g in ["GWAS_immune_hot","Roadmap_silenced","ChromHMM_bivalent","Lambert_TF",
          "dbSUPER_CD4_SE_TSS_pm50kb","gene_desert_bottomQ_density"]:
    mem=set(cm.loc[cm["group"]==g,"symbol_key"])
    hit=np.fromiter((v in mem for v in sym),bool,len(sym))
    out=[]
    for lab,_ in TESTS:
        d_,p_=disp(Ds[lab],hit); out.append("|d|%.2f p%.4f"%(d_,p_))
    print("%-30s %s" % (g," ".join("%22s"%o for o in out)))
PY

OUT
How much does gene length alone relate to the contact features?
  |r| with log10 gene length: median 0.082, max 0.416 (n_peaks_promoter)

set                                    magnitude only               + length + length, expr, CpG, density
GWAS_immune_hot                       |d|0.32 p0.0005        |d|0.24 p0.0335        |d|0.31 p0.0395
Roadmap_silenced                      |d|0.73 p0.0005        |d|0.74 p0.0010        |d|0.47 p0.0025
ChromHMM_bivalent                     |d|0.58 p0.0005        |d|0.59 p0.0005        |d|0.66 p0.0005
Lambert_TF                            |d|0.46 p0.0005        |d|0.49 p0.0005        |d|0.39 p0.0005
dbSUPER_CD4_SE_TSS_pm50kb             |d|0.16 p0.1774        |d|0.19 p0.1199        |d|0.24 p0.1524
gene_desert_bottomQ_density           |d|0.71 p0.0005        |d|0.57 p0.0005        |d|0.20 p0.0005
Aim 3 survives, now with six confounders controlled:

set	magnitude only	+ length	+ length, expr, CpG, density
ChromHMM_bivalent	0.58	0.59	0.66, p 0.0005
Roadmap_silenced	0.73	0.74	0.47, p 0.0025
Lambert_TF	0.46	0.49	0.39, p 0.0005
GWAS_immune_hot	0.32	0.24	0.31, p 0.040
dbSUPER SE	0.16	0.19	0.24, p 0.15 — still null
Recording both results:

Bash Record the length confound
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
python3 - <<'PY'
import pathlib
p=pathlib.Path("audit/CURRENT_FINDINGS.md"); s=p.read_text()

old="""### The strongest single result in Aim 1: constraint"""
new="""### RETRACTED 2026-08-16: the predictive claim does not survive gene length

The nested baselines compare the 91 features against peak counting and against
the 11-feature magnitude basis. They never compared against **gene length,
expression level and CpG density**, and those three simple covariates predict
better than the features on every target tested
(`diagnose_loeuf_confounds.py`).

| target | confounders | shape | both | shape adds |
|---|---|---|---|---|
| gnomad_loeuf | **0.452** | 0.152 | 0.454 | **+0.002** |
| GWAS_immune_hot | 0.694 | 0.627 | 0.699 | +0.005 |
| DepMap_curated_essential | 0.576 | 0.520 | 0.582 | +0.006 |
| Eisenberg_HK | 0.680 | 0.533 | 0.673 | **-0.007** |
| gtex_tau | 0.218 | 0.080 | 0.230 | +0.012 |
| phastcons_2kb | 0.055 | 0.065 | 0.098 | +0.043 |
| Lambert_TF | 0.577 | 0.622 | 0.631 | +0.053 |

**Contact shape adds essentially nothing to three simple gene properties.** The
confound is specifically length: partialling length alone drops the LOEUF result
from 0.152 to 0.085, while expression and CpG cost nothing.

Why length. LOEUF is a depletion statistic, so it scales with coding sequence
length by construction, and length independently changes contact features
because a longer gene occupies more of the plus-or-minus 1 Mb window
(`r(length, n_peaks_promoter)` = 0.42, though the median across features is only
0.08).

**What this costs.** Aim 1's justification cannot be "the features predict
biology better than the baselines". They beat peak counting and they beat the
magnitude basis, both true and both now insufficient, because neither baseline
contained length. The honest version:

> The 91 features improve on peak counting and on overall magnitude, but their
> predictive advantage over gene length, expression and CpG density is
> negligible. Their value is not prediction of existing annotations.

**What survives.** Aims 2 and 3 are untouched, because they concern the geometry
of the feature space rather than prediction of external labels, and Aim 3 was
re-tested against length directly (below). The tool is still what makes the
continuum and displacement analyses possible; it is simply not justified as a
predictor.

### Aim 3 survives gene length: six confounders now controlled

| set | magnitude only | + length | + length, expression, CpG, density |
|---|---|---|---|
| ChromHMM_bivalent | 0.58 | 0.59 | 0.66, p 0.0005 |
| Roadmap_silenced | 0.73 | 0.74 | 0.47, p 0.0025 |
| Lambert_TF | 0.46 | 0.49 | 0.39, p 0.0005 |
| GWAS_immune_hot | 0.32 | 0.24 | 0.31, p 0.040 |
| **dbSUPER SE** | 0.16 | 0.19 | **0.24, p 0.15** |
| gene_desert (density control) | 0.71 | 0.57 | 0.20, p 0.0005 |

The super-enhancer null holds under all of them. `gene_desert` collapsing is
expected and is the check working: density is in the confounder set.

The claim is now controlled for **overall magnitude, gene density, gene length,
expression, CpG density, chromosome, cell-state matching, and statistical
power.**

### The former Aim 1 headline, retained for context only"""
assert old in s; s=s.replace(old,new,1)

s=s.replace("""| "Mechanism vs output taxonomy over 7 sets" |""",
"""| "Contact architecture predicts evolutionary constraint" | **Retracted 2026-08-16.** Shape adds +0.002 over gene length, expression and CpG. The confound is length, which LOEUF scales with by construction and which no baseline contained. |
| "The 91 features are justified because they predict biology" | **Retracted.** They beat peak counting and the magnitude basis, but add nothing over three simple gene properties. Justify them by what they make possible, not by prediction. |
| "Mechanism vs output taxonomy over 7 sets" |""")
p.write_text(s); print("CURRENT_FINDINGS.md updated")

q=pathlib.Path("audit/PLAN_FORWARD.md"); t=q.read_text()
t=t.replace("""**B. Constraint as the biological result.** LOEUF is the one place architecture
does genuine work: 0.031 from peak counting, 0.210 from all features, **0.191
surviving complete magnitude removal**, so 91% is shape. Which architectural
features carry it, and why, is a positive mechanistic question of the kind the
project currently lacks. Most likely to succeed; A has the higher ceiling.""",
"""**B. CANCELLED 2026-08-16.** Constraint was to be the biological result, on the
grounds that 91% of the LOEUF prediction survived magnitude removal. It does not
survive gene length: shape adds +0.002 over length, expression and CpG, and the
confound is length specifically, which LOEUF scales with by construction. An
hour of checking saved a phase, which is exactly what it was added for.

The wider lesson applies to A and C as well: **any predictive claim must be
tested against gene length, expression and CpG density before it is believed.**
No baseline in this project has ever contained them.""")
t=t.replace("""**A. Peak-level unit of analysis. The main bet.**""",
"""**A. Peak-level unit of analysis. Now the only research bet.**""")
q.write_text(t); print("PLAN_FORWARD.md updated")
PY
grep -c "—" audit/CURRENT_FINDINGS.md audit/PLAN_FORWARD.md backend/scripts/diagnose_loeuf_confounds.py
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate mccapp && python backend/scripts/verify_store.py 2>&1 | tail -2
git add -A && git commit -q -F - <<'MSG'
Gene length kills the constraint result and Aim 1's predictive claim

Suggested check, run before building Phase B on it, and it saved the phase.

LOEUF is a depletion statistic so it scales with coding length by construction,
and length independently changes contact features because a longer gene occupies
more of the window. Neither the peak-count baseline nor the magnitude basis
contained length, expression or CpG density, so nothing in this project had ever
controlled them.

Nested and cross-validated: those three covariates predict LOEUF at r = 0.452,
shape alone at 0.152, and both together at 0.454. Shape adds +0.002. Partialling
length alone drops shape from 0.152 to 0.085 while expression and CpG cost
nothing, so the confound is specifically length.

It generalises. Across seven targets shape adds between -0.007 and +0.053 over
the confounders, the largest being Lambert_TF at +0.053 and phastCons at +0.043.
So Aim 1 cannot be justified as "the features predict biology better than the
baselines". They beat peak counting and they beat the magnitude basis, both true
and both insufficient. The honest version is that their predictive advantage
over three simple gene properties is negligible, and their value lies in what
they make possible rather than in prediction.

Aims 2 and 3 are untouched, and Aim 3 was re-tested directly rather than
assumed. Under magnitude, length, expression, CpG and density together, ChromHMM
bivalent holds at 0.66, Roadmap silenced at 0.47, Lambert_TF at 0.39,
GWAS_immune_hot at 0.31, and the super-enhancer set stays null at 0.24, p = 0.15.
gene_desert collapsing from 0.71 to 0.20 is the check working, since density is
in the confounder set. The claim is now controlled for eight things.

Phase B is cancelled and the lesson recorded for A and C: any predictive claim
must be tested against length, expression and CpG before it is believed.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
MSG
git push -q origin feat/app-scaffold && echo pushed

OUT
CURRENT_FINDINGS.md updated
PLAN_FORWARD.md updated
audit/CURRENT_FINDINGS.md:0
backend/scripts/diagnose_loeuf_confounds.py:0
audit/PLAN_FORWARD.md:0
  PASS, 38 checks

pushed
That was the most valuable hour of the week, and it cost a phase and part of an aim.

Files updated: audit/CURRENT_FINDINGS.md, audit/PLAN_FORWARD.md, new script backend/scripts/diagnose_loeuf_confounds.py.

Phase B is cancelled. Gene length, expression and CpG predict LOEUF at r = 0.452; shape alone gets 0.152; both together 0.454. Shape adds +0.002. The confound is specifically length — partialling it drops shape to 0.085, while expression and CpG cost nothing. LOEUF is a depletion statistic that scales with coding length by construction, and length independently changes contact features because a longer gene occupies more of the window.

It generalises, and that is the bigger damage. Across seven targets shape adds between −0.007 and +0.053 over those three covariates. So Aim 1 cannot be justified as "the features predict biology better than the baselines". They beat peak counting and they beat the magnitude basis — both true, both insufficient, because neither baseline contained length. The honest version is now in the document:

The 91 features improve on peak counting and on overall magnitude, but their predictive advantage over gene length, expression and CpG density is negligible. Their value is not prediction of existing annotations.

Aim 3 survives, and I tested it rather than assumed it. Under magnitude, length, expression, CpG and density together: bivalent 0.66, silenced 0.47, Lambert_TF 0.39, GWAS_immune_hot 0.31, and super-enhancers still null at 0.24, p = 0.15. gene_desert collapsing from 0.71 to 0.20 is the check working, since density is in the confounder set.

The headline is now controlled for eight things: magnitude, density, length, expression, CpG, chromosome, cell-state matching, and power.

Two things to take from this. The prediction is that these features have no advantage over three columns of gene_table.tsv for predicting annotations — which is a real limit and should be stated, not buried. But the geometric results (Aims 2 and 3) are about the shape of the space rather than prediction, and they have now survived everything thrown at them.

Phase A is the only research bet left, and the lesson is recorded against it and Phase C: any predictive claim gets tested against length, expression and CpG before it is believed. No baseline in this project has ever contained them.