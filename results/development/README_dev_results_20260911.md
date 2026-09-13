# Development results, 11 Sep 2026 (not part of the frozen matrix)

All three runs use development seeds (31.5M–31.7M), repaired engine, REF nuisance unless stated.
They inform design and interpretation; any cell the manuscript makes a claim about is to be
frozen and confirmed separately.

## 1. Studentised statistic pilot — `runs/v2/dev/studentized_pilot/` (R=100, B=199, full budget)

$T^*_{gh} = d_K(\hat S_g,\hat S_h)/(\hat v_g+\hat v_h)$, $\hat v_g$ = mean Kendall distance between the
fitted ordering and 8 consensus-level bootstrap refits (mixture fixed).

| scheme | statistic | FWER (le, α/3) | 95% Wilson | mean p |
|---|---|---|---|---|
| unrestricted | plain | 0.180 | [0.117, 0.267] | 0.307 |
| unrestricted | studentised | 0.170 | [0.109, 0.255] | 0.318 |
| diagnosis | plain | 0.040 | [0.016, 0.098] | 0.476 |
| diagnosis | studentised | 0.030 | [0.010, 0.085] | 0.482 |

Paired, dataset-level: studentised − plain under unrestricted permutation **+0.012 [−0.017, +0.041]**;
D−U shift +0.170 (plain), +0.164 (studentised). Observed $\hat v$: e2 0.079, e33 0.042, e4 0.030.

**Reading.** This candidate does not restore calibration under unrestricted permutation; it changes
nothing. Per the audit's caution: this shows one candidate failing in one configuration, not that
conditioning is necessary. It does say that a scale correction of the Behrens–Fisher type — the
remedy that works when the nuisance is purely a variance — leaves the composition effect intact,
which points at a nuisance that is not only precision (see §3). Chunks 50–99 not submitted; the
paired comparison is already tight and the queue is fair-share-limited.

## 2. Heterogeneity the model allows — `runs/v2/hpc_results/dev_heterogeneity_a/` (R=300, B=199)

Common ordering; one factor differs by group. At B=199 with the le rule a valid pair test has
attainable size 3/200 = 0.015, so the Bonferroni family ceiling is ≈ 0.045 and the max-statistic
ceiling is exactly 0.050.

| cell | what differs (e4) | U-all FWER | D-all FWER [95% Wilson] | D-all max-rej | D−U shift (n=60) |
|---|---|---|---|---|---|
| PDF_H0 | amyloid +0.45σ, p-tau +0.25σ components | 0.227 | **0.067** [0.044, 0.101] | 0.060 | +0.154 |
| MISS_CSF_H0 | CSF availability by DX, opposite in e2 | 0.223 | **0.077** [0.052, 0.113] | 0.080 | +0.149 |
| MISS_ALL_H0 | same, all modalities | 0.221 | **0.060** [0.038, 0.093] | 0.087 | +0.128 |

Reference at the same B and rule: REF D-all 0.040 [0.016, 0.098] (pilot, R=100). No failed fits.

**Reading.** Group-specific densities and group-by-diagnosis missingness push the stratified test
above its ceiling, mildly (≈1.3–1.7× of 0.045; MISS_CSF's lower bound clears the ceiling, PDF's
sits on it, MISS_ALL's does not), while the unrestricted test sits at ≈0.22 (≈5×). This is the
boundary the audit asked for: diagnosis blocks absorb the composition effect but not differences
in what the model estimates within a diagnosis. R=300 gives MCSE ≈ 0.014 at these rates; a
frozen confirmation at R=1000 is needed before any of the three numbers is quoted.

## 3. Mechanism, fixed true ordering — `runs/v2/hpc_results/dev_mechanism_fixedS_a/` (R=200, no permutations)

See `scripts/v2/dev/analyze_mechanism_fixedS.py` output, recorded below when run.

Output of `analyze_mechanism_fixedS.py` (fixed $S$ = [5,9,1,12,11,10,3,8,13,0,6,7,2,4]):

| cell | R | e2: dist / maj-d / pairs>0.5 / max | e33 | e4 |
|---|---|---|---|---|
| IID_S | 200 | 0.313 / 0.110 / 10 / 0.55 | 0.205 / 0.022 / 2 / 0.57 | 0.195 / 0.033 / 3 / 0.53 |
| REF_S | 200 | 0.384 / **0.264** / **24** / 0.69 | 0.261 / 0.176 / 16 / 0.69 | 0.176 / 0.022 / 2 / 0.62 |
| BAL_S | 200 | 0.314 / 0.066 / 6 / 0.57 | 0.205 / 0.033 / 3 / 0.54 | 0.208 / 0.033 / 3 / 0.58 |
| IID_S_N4 | 200 | 0.207 / 0.044 / 4 / 0.54 | 0.117 / 0.000 / 0 / 0.50 | 0.114 / 0.011 / 1 / 0.53 |
| REF_S_N4 | 200 | 0.306 / **0.242** / **22** / **0.84** | 0.186 / 0.099 / 9 / 0.73 | 0.096 / 0.022 / 2 / 0.70 |
| BAL_S_N4 | 200 | 0.185 / 0.011 / 1 / 0.52 | 0.114 / 0.000 / 0 / 0.49 | 0.105 / 0.000 / 0 / 0.45 |

dist = mean $d_K$(fitted, $S$); maj-d = $d_K$(pairwise-majority ordering, $S$); pairs>0.5 = pairs
reversed in the majority; max = largest per-pair reversal probability.

**Reading.** Under IID and BALANCED composition the fitted ordering converges to $S$ as $N$ grows
(majority distance → 0): a precision effect. Under REF the CN-heavy e2 group's majority ordering
stays ≈0.24 from $S$ at $4N$ and its strongest reversals *sharpen* (max 0.69 → 0.84): the
estimator concentrates on a **different ordering** — a composition-dependent estimation target,
not noise. e33 (59 % CN) shows it at about half strength; e4 (23 % CN) does not.

### Cross-check on the frozen core run (random $S$ per dataset, R=816, repaired engine)

Mean positional bias (fitted position − true position; negative = fitted earlier):

| biomarker | IID e2 / e33 / e4 | REF e2 / e33 / e4 |
|---|---|---|
| ABETA | +0.37 / +0.33 / +0.11 | **+1.60 / +1.23 / −0.84** |
| PTAU | −0.08 / +0.20 / +0.08 | +0.78 / +0.83 / −0.57 |
| NG | −0.23 / −0.13 / −0.08 | **−1.63 / −1.51 / +0.74** |
| NFL | −0.16 / −0.25 / −0.09 | −0.99 / −1.21 / +0.72 |
| ADAS13 | +0.40 / +0.21 / +0.20 | **+2.77 / +1.62 / −1.08** |
| MMSE | +0.06 / +0.19 / +0.18 | +1.26 / +0.98 / −0.66 |
| WholeBrain | −0.09 / −0.06 / −0.16 | −1.18 / −0.62 / +0.61 |
| Ventricles | −0.34 / −0.46 / −0.27 | **−1.75 / −1.26 / +0.71** |

Under a shared true ordering, CN-heavy and AD-heavy groups are pulled in **opposite** directions
along the separability / availability axis: high-separability markers (ABETA, ADAS13, MMSE,
PTAU) come out *late* in CN-heavy groups and *early* in AD-heavy ones; low-separability or
low-availability markers (NG, NFL, Ventricles, WholeBrain) the reverse. Between e2 and e4 the
induced difference in target position is 2.4 positions for ABETA and NG, 3.9 for ADAS13 — with
no true difference.

**Consequences.** (i) The nuisance is a target shift, not only a precision loss, which is why the
studentised statistic (a scale correction) leaves the inflation intact. (ii) Diagnosis-stratified
permutation remains valid because pseudo-groups keep their composition and hence their targets;
what it tests is equality of the composition-specific targets. (iii) The published ε2 finding
("neurogranin among the earliest events; Aβ42 after p-tau") lies on this axis and in this
direction for a CN-heavy group; whether it is composition-induced is testable on ADNI by
subsampling the ε4 group to ε2's composition and refitting. This is a development result with
placeholder separabilities; it is not a statement about the real cohort until that test is run.

### Localising the bias: engines and oracles (frozen mechanism_b, REF_H0 and IID_H0, R=200)

Mean |e2 − e4| target-position difference over the 14 biomarkers, same true ordering:

| cell | original | repaired | oracle equal-prior | **oracle known density + estimated prior** |
|---|---|---|---|---|
| IID_H0 | 0.24 | 0.24 | 0.50 | 0.17 |
| REF_H0 | **1.60** | **1.59** | 2.45 | **0.21** |

Core_a REF, original vs repaired (R=816): ABETA e2 +1.63 / +1.60, NG e2 −1.68 / −1.63 — identical.

**Reading.** The composition-dependent target is (i) not the optimiser defect (original = repaired),
(ii) not removed by a fixed equal prior (worse), and (iii) **removed by supplying the true
component densities** (REF gap 0.21 ≈ IID). It arises in the per-group mixture estimation of
co-init DEBM: in a CN-heavy group the abnormal component of low-separability / low-availability
biomarkers is poorly identified, normal subjects receive inflated abnormality posteriors, and
those events are placed early; in an AD-heavy group the normal component is the poorly
identified one and the pattern reverses.

**Estimator-level implication.** Sharing one mixture across groups (fit on the pooled sample)
and estimating only the ordering per group would remove the composition dependence by
construction — the "coupled" DEBM variant that the 2021 paper set aside in favour of co-init on
the basis of simulations in which groups *did* differ in densities. When groups share densities
but differ in composition, the trade goes the other way. This is a development question for a
new run, not a claim.

## 5. Shared-mixture DEBM — `runs/v2/dev/shared_mixture/` (dev seeds 31.8M; analysis in `results/v2_dev/shared_mixture_analysis_20260911T1612.txt`)

Pooled mixture fitted once on all subjects, replicated to the three groups; only the Mallows
consensus uses labels. Mechanism: fixed $S$, R=200 per cell, all 1,200 jobs ok. Calibration:
REF_H0, R=300 planned, **87 done at 16:20 UTC 09-11**, B=199 per scheme, cache-consistency
check passed on every dataset.

### Mechanism — same table as §3, per-group vs shared mixture

| cell | per-group e2: maj-d / pairs>0.5 / max | shared e2 | per-group \|e2−e4\| gap | shared gap |
|---|---|---|---|---|
| IID_S | 0.110 / 10 / 0.55 | 0.011 / 1 / 0.52 | 0.98 | 0.10 |
| REF_S | **0.264 / 24 / 0.69** | **0.132 / 12 / 0.61** | **2.05** | **0.88** |
| BAL_S | 0.066 / 6 / 0.57 | 0.011 / 1 / 0.51 | 0.86 | 0.10 |
| IID_S_N4 | 0.044 / 4 / 0.54 | 0.011 / 1 / 0.53 | 0.61 | 0.05 |
| REF_S_N4 | **0.242 / 22 / 0.84** | **0.110 / 10 / 0.76** | **2.22** | **1.17** |
| BAL_S_N4 | 0.011 / 1 / 0.52 | 0.000 / 0 / 0.47 | 0.48 | 0.07 |

(gap = mean over 14 biomarkers of |mean fitted position e2 − mean fitted position e4|, fixed $S$.)
Mean $d_K$ to $S$, REF_S: e2 0.384 → 0.262 (better), e33 0.261 → 0.237, **e4 0.176 → 0.217 (worse)**;
at 4N e4 0.096 → 0.139.

**Reading.** Sharing the mixture removes the composition-dependent target **under IID and
BALANCED entirely** and **about half of it under REF** (gap 2.05 → 0.88; majority-reversed pairs
24 → 12). The remaining half is a target, not noise: it persists and grows at 4N (gap 1.17,
max reversal 0.76). It therefore lives in the **consensus step** — the weighted Mallows
aggregate of subject orderings targets a different ordering when the subjects are mostly CN
than when they are mostly AD, even with identical posterior models. Shared mixture improves
the CN-heavy group and *harms* the AD-heavy group.

### Calibration, REF_H0 (87/300, B=199, dev rule: pair rejects iff #{T_b ≥ T_obs} ≤ 2)

| scheme | pair rej e2-e33 / e2-e4 / e33-e4 | family any-pair | max-stat | ceiling |
|---|---|---|---|---|
| unrestricted (published) | 0.069 / **0.483** / **0.632** | **0.782 [0.695, 0.868]** | 0.632 | 0.045 / 0.050 |
| diagnosis-stratified | 0.000 / 0.011 / 0.000 | 0.011 [0.000, 0.034] | 0.023 | 0.045 / 0.050 |

p-shift (D − U): e2-e33 +0.12, e2-e4 +0.49, e33-e4 +0.56. Mean observed $d_K$: 0.095 / 0.159 / 0.102.

**Reading.** With a shared mixture the published (unrestricted) permutation becomes **far more
anti-conservative** (family FWER ≈ 0.78 vs ≈ 0.11 with per-group mixtures): the pooled mixture
removes the refit noise that inflated the reference distribution, and the residual
composition-dependent target difference (0.88 positions between e2 and e4) is now detected as
"different" in most no-difference cohorts. In the per-group design that noise was partly
masking the target shift. Diagnosis-stratified permutation stays calibrated (0.011 / 0.023).

**Consequence for an estimator-level fix.** Shared mixture alone is *not* a fix; it is a
sharper estimator that exposes the consensus-step target shift. A complete estimator-level fix
must also standardise the consensus to a common composition (e.g. weight subjects so every
group has the reference diagnostic composition inside the weighted Mallows step, or subsample
to a matched composition). Untested; proposed as the next dev run.

## 6. Composition-standardised consensus — `runs/v2/dev/stdcons/` (dev seeds 31.8M; analysis in `results/v2_dev/stdcons_analysis_20260911T1637.txt`)

Subject weights $w_i = \pi_\text{ref}(d_i)/\pi_g(d_i)$ enter the Mallows objective as a weighted
mean over subjects (`scripts/v2/dev/standardized_consensus_dev.py`); $\pi_\text{ref}$ = pooled
composition ("pooled") or renormalised element-wise minimum over groups ("min"). Weights are
recomputed from the (pseudo-)labels at every fit. Uniform weights reproduce the plain fit exactly.
Mechanism: fixed $S$, 6 cells × 200 × 3 variants, all 3,600 jobs ok (job 10356899, 16:26–17:10 UTC).
Calibration (job 10356900) running.

### Mechanism — mean |e2−e4| target-position gap (14 biomarkers, fixed $S$)

| cell | per-group mixture | shared mixture | per-group + std(pooled) | **shared + std(pooled)** | **shared + std(min)** | IID reference |
|---|---|---|---|---|---|---|
| REF_S | 2.05 | 0.88 | 1.55 | **0.16** | **0.12** | 0.10–0.12 |
| REF_S_N4 | 2.22 | 1.17 | 1.45 | **0.09** | **0.08** | 0.05 |

REF_S e2: majority-reversed pairs 24 → 12 (shared) → 16 (std only) → **2 / 2** (shared + std);
at 4N 22 → 10 → 11 → **1 / 1**; max per-pair reversal 0.84 → 0.76 → 0.70 → **0.55 / 0.60** at 4N.
Mean $d_K$ to $S$, REF_S (e2 / e33 / e4): per-group 0.384 / 0.261 / 0.176; shared + std(pooled)
**0.245 / 0.224 / 0.225**; shared + std(min) 0.249 / 0.232 / 0.233. ESS (pooled ref) e2 45 of 75,
e33 367 of 411, e4 397 of 485; (min ref) 53 / 388 / 348.

**Reading.** Neither step alone suffices: sharing the mixture removes about half of the REF
target gap, standardising the consensus alone removes about a quarter; **together they remove
it entirely** — the REF gap falls to the IID level and stays there at 4N, the majority ordering
of every group coincides with $S$ up to 1–2 pairs, and the three groups reach the same accuracy
(≈0.23–0.25) regardless of composition. Cost: the AD-heavy group loses some accuracy relative
to the per-group fit (0.176 → 0.225 at N, 0.096 → 0.142 at 4N), inherited from the shared
mixture, and the CN-heavy group's effective sample size drops by ~40 %. The two reference
compositions are equivalent in target; "min" has lower weight variance.

**Status.** Composition-invariant estimation is established at the mechanism level in the
development design (identical component densities across groups by construction). Whether the
resulting group orderings support valid inference under unrestricted and stratified permutation
is the calibration run (REF_H0, R=300, B=199), pending. Expected: D calibrated; U still inflated
by the remaining precision (ESS) difference.

## 7. SA-EBM as an alternative estimator — complete-data comparison (`runs/v2/dev/stdcons_nomiss/`, `runs/v2/dev/saebm/`; analysis `results/v2_dev/saebm_comparison_20260911T1835.txt`)

pysaebm (Hao et al., ML4H 2025) needs a complete participant × biomarker matrix, so all three
estimators were run on the **same cohorts simulated with `missing=False`** (fixed $S$, dev seeds,
R=200 per cell, all jobs ok). SA-EBM: `conjugate_priors`, n_iter 10,000, burn-in 2,500, one fit per
group, `diseased` = not CN. Jobs 10360492 and 10360746.

| estimator | cell | e2 dist / maj-d / pairs>0.5 | e33 dist | e4 dist | \|e2−e4\| gap | s/dataset |
|---|---|---|---|---|---|---|
| per-group DEBM | REF_S | 0.360 / 0.253 / 23 | 0.233 | 0.161 | 1.98 | 83 |
| per-group DEBM | REF_S_N4 | 0.280 / 0.231 / 21 | 0.165 | 0.077 | 2.12 | 76 |
| **shared + std (ours)** | REF_S | 0.218 / 0.022 / 2 | 0.201 | 0.200 | **0.14** | 16 |
| **shared + std (ours)** | REF_S_N4 | 0.140 / 0.011 / 1 | 0.130 | 0.130 | **0.05** | 26 |
| SA-EBM | IID_S | 0.287 / 0.077 / 7 | 0.084 | 0.073 | 1.50 | 55 |
| SA-EBM | REF_S | 0.372 / 0.165 / 15 | 0.113 | 0.062 | 2.27 | 54 |
| SA-EBM | REF_S_N4 | 0.210 / 0.011 / 1 | 0.034 | 0.014 | 0.94 | 209 |

**Reading.**
1. **SA-EBM is a much better single-group estimator when the group has enough progressing
   subjects**: e33/e4 at REF reach 0.11/0.06 (DEBM 0.23/0.16; ours 0.20/0.20), and 0.03/0.01 at 4N.
2. **It does not make separately fitted orderings comparable at realistic sizes.** The CN-heavy
   e2 group (75 subjects, 6 MCI) is as poorly recovered as under DEBM (0.372 vs 0.360), its majority
   ordering is 0.165 from $S$ with 15 majority-reversed pairs, and the e2–e4 gap is 2.27 — larger
   than DEBM's. Unlike DEBM, this is mostly a **precision** effect: at 4N the e2 majority ordering
   converges (0.011) and the gap falls to 0.94, so SA-EBM is consistent where co-init DEBM is not.
   But the between-group nuisance at ADNI-like N is then exactly the Behrens–Fisher situation
   (very unequal precision), which the unrestricted permutation test cannot absorb.
3. **Ours equalises groups rather than maximising per-group accuracy**: all three groups land at
   ≈0.20–0.22 (N) and ≈0.13–0.14 (4N) with gap ≤0.14; the pooled mixture lends the small CN-heavy
   group the strength it lacks, at a cost on the large groups relative to SA-EBM.
4. Practical: SA-EBM cannot take missing data (ADNI CSF availability ≈ 17 % missing in the REF
   design), and is ~4–8× slower per fit at 4N.

**Position for the paper.** Cite SA-EBM as convergent evidence that composition degrades EBM
orderings and as the best available single-group estimator; show that it leaves the
between-group problem in place at realistic sizes (gap 2.27 under REF); recommend our
composition-invariant estimator **for comparison purposes** and conditional inference in either
case. A natural extension — SA-EBM's stage-aware likelihood with a pooled measurement model and
standardised aggregation — is future work.

### §6 addendum — calibration final (300/300, 21:0x UTC 09-11; job 10356900; shared-mixture job 10353054 also 300/300)

REF_H0, B=199 per scheme, dev rule (pair rejects iff #{T_b ≥ T_obs} ≤ 2; family = any pair; max iff ≤ 9):

| estimator | scheme | pair rej e2-e33 / e2-e4 / e33-e4 | family any-pair [95 % CI] | max-stat |
|---|---|---|---|---|
| shared mixture only | unrestricted | 0.083 / 0.547 / 0.633 | **0.783 [0.737, 0.830]** | 0.677 |
| shared mixture only | diagnosis | 0.007 / 0.010 / 0.007 | 0.017 [0.002, 0.031] | 0.050 |
| shared + std (pooled) | unrestricted | 0.043 / 0.050 / 0.023 | **0.097 [0.063, 0.130]** | 0.100 |
| shared + std (pooled) | diagnosis | 0.003 / 0.010 / 0.010 | 0.020 [0.004, 0.036] | 0.027 |
| shared + std (min) | unrestricted | 0.027 / 0.043 / 0.017 | **0.077 [0.047, 0.107]** | 0.077 |
| shared + std (min) | diagnosis | 0.007 / 0.007 / 0.003 | 0.013 [0.000, 0.026] | 0.023 |

Ceilings: pair 0.015, family 0.045, max 0.050. p-shift (D − U): shared-only +0.12 / +0.46 / +0.52;
shared + std ≈ +0.13 / +0.14 / +0.05–0.07. Mean observed d_K (e2-e4): shared-only 0.170; shared + std 0.093.
Cache-consistency check passed on all 900 datasets. ESS (pooled) e2 45 / e33 367 / e4 397; (min) 53 / 388 / 348.

**Reading (final).** (i) The composition-invariant estimator brings the published (unrestricted)
procedure from 0.78 (shared-only, target shift fully exposed) down to 0.08–0.10 — the residual is
the precision (ESS) difference between groups, and the 'min' reference, which has the smaller
weight variance, sits closer to the ceiling. It is still above 0.045, so the estimator alone does
not make the published test valid. (ii) Diagnosis-stratified permutation is calibrated for every
estimator (family 0.013–0.020; max 0.023–0.050). (iii) Recommended default: shared mixture +
standardised consensus with the 'min' reference (lower weight variance, same target), followed by
diagnosis-stratified permutation. Both layers are necessary: see the 0.78 row for what happens
with the estimator half-fixed, and the 0.08–0.10 rows for what remains after the estimator is fixed.
