# audit/scripts

Every analysis behind `audit/CURRENT_FINDINGS.md`. All of them read the **built
store**, none of them touch the pipeline, and all are safe to re-run.

## Re-running after the feature space changes

This is the reason the folder exists. Changing `mccprofiler`'s features moves the
component count, the variance shares, the trust weighting and every effect size,
so nothing in `CURRENT_FINDINGS.md` carries over.

```bash
conda activate cd4env
python backend/scripts/build_store.py          # rebuild from the new features
conda activate mccapp
python backend/scripts/verify_store.py         # WILL fail, see below
conda activate cd4env
bash audit/scripts/run_all.sh                  # ~30 min, writes to audit/scripts/output/
```

`verify_store.py` failing after a feature change is correct behaviour, not a
problem. Its `EXPECTED` constants encode the current numbers. Update them
deliberately once you have read why each moved; never silently re-baseline.

Then re-read `CURRENT_FINDINGS.md` against the new outputs. Section 7, the
retraction table, is the part to preserve: those claims were wrong for reasons
that do not change when the features do.

## Order, and why it is this order

**1. The substrate.** What is in the feature space and how far it can be
trusted. Run first because everything downstream is conditional on it.

| script | question |
|---|---|
| `diagnose_feature_novelty` | do gene length, expression and CpG explain the features? (2.6%) |
| `diagnose_dimension_trust` | which components rest on features that reproduce |
| `diagnose_pc_names` | how much of each axis name is actually true |

**2. Aim 2, is the structure categorical.** Depends on the substrate only.

| script | question |
|---|---|
| `experiment_cluster_search` | 16 conditions, four cluster tests each |
| `experiment_rotate_axes` | would varimax make the axes nameable, and at what cost |

**3. Aim 3, external sets.** Ordered by the confound each one adds. Run in
sequence: each assumes the previous correction is in place.

| script | adds |
|---|---|
| `diagnose_external_structure_corrected` | magnitude (`MAG_OVERALL`) |
| `diagnose_external_density_stratified` | gene density |
| `diagnose_chromosome_confound` | chromosome |
| `diagnose_external_power` | minimum detectable effect, so a null is interpretable |
| `diagnose_external_axis_rank` | which component each set leans on |
| `diagnose_definition_type` | does the KIND of definition predict displacement |
| `experiment_feature_search` | can ANY feature subset separate them, scored against permuted labels |
| `diagnose_promoter_assembly` | splits promoter contacts into own-TSS vs other genes', tests the HK assembly hypothesis |

`diagnose_external_structure.py` is the superseded version that ran on raw
components. Kept for provenance, not for quoting.

**4. Predictive claims.** `diagnose_loeuf_confounds` exists because no baseline
in this project ever contained gene length, and a whole phase was nearly built
on a length artefact. **Any new predictive claim goes through it first.**

**5. Formerly ad hoc.** `diagnose_adhoc_checks` captures four analyses that were
run interactively and each changed a conclusion: panel coverage, taxonomy
robustness, the Lambert DBD split, and whether chr19 is really gene density.

**6. Figures.** `make_audit_figures` writes six composites into `audit/figures/`.
Draw them before writing prose: figure 6 falsified a claim four rounds of written
review had missed, and figure 3 exposed an arithmetic error.

## Environments

| | |
|---|---|
| `cd4env` | everything here (needs sklearn, scipy, hdbscan, diptest, matplotlib) |
| `mccapp` | `verify_store.py` only, and the server |

`mccapp` is deliberately thin and has no sklearn. If a script fails on an
import, check which environment is active before anything else.

## External dependencies

Beyond the store, these scripts read:

- `audit/GW/scripts/_shape.py` for `MAG_OVERALL`, the 11-feature magnitude basis
- `collaboration/data/gene_table.tsv` for length, expression, CpG, density, chromosome
- `chromprofiler/benchmarking/src` for the Lambert DBD lookup
- `audit/continuous_methods/*.tsv` for nested baselines and summit precision

All are absolute paths into `scripts_cleaned`, by project convention.
