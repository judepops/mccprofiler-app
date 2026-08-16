Plan's good. One ordering error worth catching before it costs a day.

S1 should run after the feature removals, not before.

Current order is S1 → S2 (including both removals) → figures. But both removals change the substrate — component count, variance shares, every displacement value. So the core-promoter displacement test would run on the old feature set and then need re-running on the new one. You'd do S1 twice.

Revised:

S0  orthogonality               DONE
S2.6 PROJECT_STATUS redirect    (one hour, outward-facing risk)
S2.1 + S2.1b feature removals   (both)
S2.3 rebuild store              
S1  core promoters              ← now on the final substrate
S2.2, S2.4, S2.5                remaining corrections
S3  figures
S4  write

The only thing that must stay at the very front is S2.6, and for that one, make the redirect the choice rather than the alternative — point PROJECT_STATUS.md at CURRENT_FINDINGS.md with a dated note saying the numbers are being re-derived. Writing real numbers into it now guarantees writing them twice.

One note on S2.1b. "Re-support or retire" leaves the decision open, and it should probably be decided before you spend two days on it. The re-supporting route (stacked peaks per gene, or per element class) is worth trying because contact_asymmetry at ρ 0.900 proves the concept works with adequate support. But if the aggregated version doesn't clear ρ 0.7, retire them rather than iterating — 14.9% of variance from features you can't measure reliably is worth less than a smaller, trustworthy substrate.

Set that threshold now, before running it, for the same reason every other pre-registration in this project has earned its keep.

Everything else I'd leave alone. The stopping-rule clarification and the Phase A pre-registration are both exactly right.