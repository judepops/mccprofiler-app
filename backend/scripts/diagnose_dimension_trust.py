"""Which dimensions rest on measurements that reproduce?

The naming diagnostics ask whether an axis is ABOUT one thing. This asks the
prior question: is the axis built on things we can measure twice and get the
same answer? A perfectly nameable axis carried by features with rho 0.3 is a
well-described artefact.

The number is the share of each component's squared loading carried by features
with cross-panel rho > 0.7, restricted to the 63 features measured in both
panels. Coverage is reported alongside, because a component whose loading mass
sits mostly on the 28 unmeasured features cannot be judged either way.

This exists because varimax rotation, which nearly doubles nameability for free,
concentrates its cleanest factors on the oe_asymmetry family, and that family is
the least reproducible thing in the feature set (rho 0.27 to 0.42). Rotation
should be weighted by this before it is adopted.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from app.store import get_store

TRUST = 0.7

s = get_store()
rep = s.table("reproducibility_pairs")
rho = {}
for f, sub in rep.groupby("feature"):
    if len(sub) > 10 and sub["gw"].std() > 0 and sub["immune"].std() > 0:
        rho[f] = float(np.corrcoef(sub["gw"].rank(), sub["immune"].rank())[0, 1])

L = s.table("pc_loadings")
scree = s.table("pc_scree")
print(f"{len(rho)} features measured in both panels, "
      f"{sum(v > TRUST for v in rho.values())} above rho {TRUST}\n")

print(f"{'PC':<5}{'var%':>8}{'coverage':>10}{'trusted':>9}   top loading (rho)")
rows = []
for r in scree.itertuples():
    pc = int(r.pc)
    sub = L[L["pc"] == pc].copy()
    if sub.empty:
        continue
    sub["rho"] = sub["feature"].map(rho)
    w_all = sub["loading"] ** 2
    m = sub.dropna(subset=["rho"])
    if m.empty:
        continue
    coverage = float(w_all[sub["feature"].isin(rho)].sum() / w_all.sum())
    w = (m["loading"] ** 2) / (m["loading"] ** 2).sum()
    trusted = float((w * (m["rho"] > TRUST)).sum())
    top = m.loc[m["loading"].abs().idxmax()]
    rows.append(dict(pc=pc, variance_pct=r.variance_pct, above_noise=bool(r.above_noise),
                     coverage=coverage, trusted=trusted,
                     top_feature=top["feature"], top_rho=top["rho"]))
    if pc <= 12:
        print(f"{pc:<5}{r.variance_pct:>7.2f}%{coverage:>9.0%}{trusted:>9.0%}"
              f"   {top['feature'][:32]} ({top['rho']:.2f})")

df = pd.DataFrame(rows)
real = df[df["above_noise"]]
print(f"\nAcross the {len(real)} components above the noise ceiling:")
print(f"  weighted by variance, {(real['trusted'] * real['variance_pct']).sum() / real['variance_pct'].sum():.0%}"
      f" of the retained structure rests on reproducible features")
weak = real[real["trusted"] < 0.5].sort_values("variance_pct", ascending=False)
if len(weak):
    print(f"\n  {len(weak)} component(s) below 50% trusted, carrying "
          f"{weak['variance_pct'].sum():.1f}% of variance:")
    for x in weak.itertuples():
        print(f"    PC{x.pc:<3} {x.variance_pct:5.2f}%  trusted {x.trusted:.0%}  "
              f"top {x.top_feature} (rho {x.top_rho:.2f})")

print("\nLeast reproducible features, which is where those components live:")
for f, v in sorted(rho.items(), key=lambda kv: kv[1])[:6]:
    print(f"  {v:.3f}  {f}")
