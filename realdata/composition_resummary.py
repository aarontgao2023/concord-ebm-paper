#!/usr/bin/env python3
"""Six-event panel, step 6a: mean-position gap and event shifts of the controlled composition
experiment (Fig. 4b and 4d).

Standard library only. No model is fitted and no participant-level file is read. The inputs are
the fit-level event positions written by summarize_common_panel.py:
composition_event_positions.csv and composition_fits.csv (200 draws x 2 compositions x 4 models
x 2 groups, per cohort).

For each cohort, model and composition (matched, different) and each event, the shift is the
position in group B minus the position in group A (1-based positions), averaged over the 200
draws. In the different composition A is the CN-heavy group (112/32/16) and B the AD-heavy group
(48/64/48); Fig. 4d plots the opposite sign (CN-heavy minus AD-heavy). The mean-position gap is
the mean over the six events of the absolute mean shift. Intervals use 10,000 paired bootstrap
resamples of the 200 draws with Python's random.Random(seed), seed 20260924 (ADNI) and 20260925
(NACC); within a cohort the same resampled draw indices are used for every model, composition
and event, and the paired increment gap(different) - gap(matched) is recomputed in every
resample. The sign-flip columns of real_comp_gap.csv (10,000 flips of each draw's shift vector,
seeds 20260924 + 1000 and 20260925 + 1000) were an internal check and are not reported in the
paper. These summaries were defined after the fits.

Model labels: likelihood_ebm = likelihood EBM, separate = separately fitted DEBM,
shared = pooled-score DEBM, concord = CONCORD.

  python composition_resummary.py --summary-dir $CONCORD_WORK_DIR/summary --out $CONCORD_WORK_DIR/resummary

Outputs: real_comp_event_shifts.csv, real_comp_gap.csv, real_comp_resummary_meta.json.
"""
import argparse, csv, json, math, random, operator, os, time
from pathlib import Path

B = 10000
SEED = {'adni': 20260924, 'nacc': 20260925}
FLIP_B = 10000
EVENTS = ['ABETA', 'TAU', 'PTAU', 'MEM', 'EXF', 'LAN']
ENGINES = ['likelihood_ebm', 'separate', 'shared', 'concord']
DESIGNS = ['matched', 'different']
COHORTS = ['adni', 'nacc']
R = 200


def pct(sorted_vals, q):
    # linear-interpolated percentile (numpy default 'linear')
    n = len(sorted_vals)
    h = (n - 1) * q
    lo = math.floor(h); hi = math.ceil(h)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (h - lo)


def ci(vals):
    s = sorted(vals)
    return pct(s, 0.025), pct(s, 0.975)


def mean(v):
    return sum(v) / len(v)


def sd(v):
    m = mean(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--summary-dir', type=Path, required=True,
                    help='folder with composition_event_positions.csv and composition_fits.csv')
    ap.add_argument('--out', type=Path, required=True, help='output folder')
    args = ap.parse_args(argv)
    POS = args.summary_dir / 'composition_event_positions.csv'
    FITS = args.summary_dir / 'composition_fits.csv'
    OUT = args.out
    OUT.mkdir(parents=True, exist_ok=True)

    # ---------- load ----------
    pos = {}  # (cohort, rep, design, engine, group) -> {event: pos}
    with open(POS) as f:
        for row in csv.DictReader(f):
            k = (row['cohort'], int(row['replicate']), row['design'], row['engine'], row['group'])
            pos.setdefault(k, {})[row['event']] = int(row['position_1_based'])
    assert len(pos) == 2 * R * 2 * 4 * 2, len(pos)
    for k, d in pos.items():
        assert sorted(d) == sorted(EVENTS), k
        assert sorted(d.values()) == [1, 2, 3, 4, 5, 6], k

    # cross-check against composition_fits orders and normalized Kendall distances
    checked = 0
    with open(FITS) as f:
        for row in csv.DictReader(f):
            c, r, dsg, eng = row['cohort'], int(row['replicate']), row['design'], row['engine']
            assert row['status'] == 'ok'
            for g, col in (('A', 'group_A_event_order'), ('B', 'group_B_event_order')):
                order = json.loads(row[col])
                p = pos[(c, r, dsg, eng, g)]
                assert all(p[e] == i + 1 for i, e in enumerate(order)), (c, r, dsg, eng, g)
            pa, pb = pos[(c, r, dsg, eng, 'A')], pos[(c, r, dsg, eng, 'B')]
            disc = sum(1 for i in range(6) for j in range(i + 1, 6)
                       if (pa[EVENTS[i]] - pa[EVENTS[j]]) * (pb[EVENTS[i]] - pb[EVENTS[j]]) < 0)
            assert abs(disc / 15 - float(row['normalized_kendall'])) < 1e-12
            checked += 1
    assert checked == 3200
    print('input checks pass: 3200 fits, 38400 positions, orders and Kendall distances reconcile')

    # ---------- per-draw series ----------
    # series[(cohort, engine, design, name)] = list of R values indexed by draw
    series = {}
    for c in COHORTS:
        for eng in ENGINES:
            for dsg in DESIGNS:
                A = [pos[(c, r, dsg, eng, 'A')] for r in range(R)]
                Bg = [pos[(c, r, dsg, eng, 'B')] for r in range(R)]
                for e in EVENTS:
                    series[(c, eng, dsg, 'shift_' + e)] = [Bg[r][e] - A[r][e] for r in range(R)]
                    series[(c, eng, dsg, 'posA_' + e)] = [A[r][e] for r in range(R)]
                    series[(c, eng, dsg, 'posB_' + e)] = [Bg[r][e] for r in range(R)]

    # ---------- paired bootstrap of draw means ----------
    keys_by_cohort = {c: [k for k in series if k[0] == c] for c in COHORTS}
    boot = {}  # key -> list of B bootstrap means
    t0 = time.time()
    for c in COHORTS:
        rng = random.Random(SEED[c])
        keys = keys_by_cohort[c]
        mat = [series[k] for k in keys]
        out = [[0.0] * B for _ in keys]
        for b in range(B):
            cnt = [0] * R
            for _ in range(R):
                cnt[rng.randrange(R)] += 1
            for i, row in enumerate(mat):
                out[i][b] = sum(map(operator.mul, row, cnt)) / R
        for i, k in enumerate(keys):
            boot[k] = out[i]
    print(f'bootstrap done in {time.time() - t0:.1f}s')

    # ---------- sign-flip reference for the gap (internal check) ----------
    flip = {}
    for c in COHORTS:
        rng = random.Random(SEED[c] + 1000)
        signs = [[rng.choice((-1, 1)) for _ in range(R)] for _ in range(FLIP_B)]
        for eng in ENGINES:
            for dsg in DESIGNS:
                sh = [series[(c, eng, dsg, 'shift_' + e)] for e in EVENTS]
                gaps = []
                for s in signs:
                    ms = [sum(map(operator.mul, v, s)) / R for v in sh]
                    gaps.append(mean([abs(x) for x in ms]))
                flip[(c, eng, dsg)] = gaps

    # ---------- table 1: per-event signed shifts (B - A) ----------
    rows = []
    for c in COHORTS:
        for eng in ENGINES:
            for dsg in DESIGNS:
                for e in EVENTS:
                    v = series[(c, eng, dsg, 'shift_' + e)]
                    m = mean(v)
                    lo, hi = ci(boot[(c, eng, dsg, 'shift_' + e)])
                    npos = sum(x > 0 for x in v); nneg = sum(x < 0 for x in v); nz = R - npos - nneg
                    same = npos if m > 0 else (nneg if m < 0 else nz)
                    opp = nneg if m > 0 else (npos if m < 0 else 0)
                    rows.append(dict(cohort=c, engine=eng, design=dsg, event=e, n_draws=R,
                                     mean_pos_A=round(mean(series[(c, eng, dsg, 'posA_' + e)]), 4),
                                     mean_pos_B=round(mean(series[(c, eng, dsg, 'posB_' + e)]), 4),
                                     mean_shift_B_minus_A=round(m, 4), sd_shift=round(sd(v), 4),
                                     ci95_low=round(lo, 4), ci95_high=round(hi, 4),
                                     ci_excludes_0=int(lo > 0 or hi < 0),
                                     frac_later_in_B=round(npos / R, 4), frac_earlier_in_B=round(nneg / R, 4), frac_no_change=round(nz / R, 4),
                                     frac_same_sign_as_mean=round(same / R, 4), frac_opposite_sign=round(opp / R, 4),
                                     same_over_nonzero=round(same / (npos + nneg), 4) if npos + nneg else ''))
    with open(OUT / 'real_comp_event_shifts.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    # ---------- table 2: gap and different-minus-matched increment ----------
    def gap_from(means):
        return mean([abs(x) for x in means])

    def rms_debiased(c, eng, dsg):
        # mean over events of (mean^2 - var/n): unbiased for mean squared true shift
        vals = []
        for e in EVENTS:
            v = series[(c, eng, dsg, 'shift_' + e)]
            vals.append(mean(v) ** 2 - sd(v) ** 2 / R)
        return math.sqrt(max(0.0, mean(vals)))

    rows = []
    gap_boot = {}
    for c in COHORTS:
        for eng in ENGINES:
            for dsg in DESIGNS:
                bm = [boot[(c, eng, dsg, 'shift_' + e)] for e in EVENTS]
                gb = [mean([abs(bm[i][b]) for i in range(6)]) for b in range(B)]
                gap_boot[(c, eng, dsg)] = gb
        for eng in ENGINES:
            for dsg in DESIGNS:
                g = gap_from([mean(series[(c, eng, dsg, 'shift_' + e)]) for e in EVENTS])
                lo, hi = ci(gap_boot[(c, eng, dsg)])
                fg = flip[(c, eng, dsg)]
                fs = sorted(fg)
                p_flip = (sum(x >= g - 1e-12 for x in fg) + 1) / (FLIP_B + 1)
                per_rep = mean([mean([abs(series[(c, eng, dsg, 'shift_' + e)][r]) for e in EVENTS]) for r in range(R)])
                row = dict(cohort=c, engine=eng, design=dsg, gap_mean_abs_mean_shift=round(g, 4),
                           boot_ci95_low=round(lo, 4), boot_ci95_high=round(hi, 4),
                           boot_mean=round(mean(gap_boot[(c, eng, dsg)]), 4),
                           signflip_null_median=round(pct(fs, 0.5), 4), signflip_null_p95=round(pct(fs, 0.95), 4),
                           signflip_p_one_sided=round(p_flip, 5),
                           noise_corrected_rms_shift=round(rms_debiased(c, eng, dsg), 4),
                           mean_per_draw_abs_shift=round(per_rep, 4))
                if dsg == 'different':
                    gm = gap_from([mean(series[(c, eng, 'matched', 'shift_' + e)]) for e in EVENTS])
                    diffb = [a - b for a, b in zip(gap_boot[(c, eng, 'different')], gap_boot[(c, eng, 'matched')])]
                    dlo, dhi = ci(diffb)
                    row.update(diff_minus_matched_gap=round(g - gm, 4), dmm_ci95_low=round(dlo, 4), dmm_ci95_high=round(dhi, 4),
                               dmm_frac_boot_gt0=round(sum(x > 0 for x in diffb) / B, 4))
                else:
                    row.update(diff_minus_matched_gap='', dmm_ci95_low='', dmm_ci95_high='', dmm_frac_boot_gt0='')
                rows.append(row)
    with open(OUT / 'real_comp_gap.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    meta = dict(source=POS.name, source_fits=FITS.name, bootstrap_resamples=B, bootstrap_seeds=SEED, signflip_resamples=FLIP_B,
                signflip_seeds={k: v + 1000 for k, v in SEED.items()}, percentile='linear interpolation, pointwise 2.5/97.5',
                pairing='same draw indices for all models/compositions/events/statistics within cohort',
                group_labels='different composition: A=112/32/16 (CN-heavy), B=48/64/48 (AD-heavy); matched: both 80/48/32',
                status='descriptive summary of stored positions, defined after the fits; no refit')
    with open(OUT / 'real_comp_resummary_meta.json', 'w') as f:
        json.dump(meta, f, indent=1)
    print('done')


if __name__ == '__main__':
    main()
