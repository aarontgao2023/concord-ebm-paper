#!/usr/bin/env python3
"""MRI-cognition populations: covariate adjustment, paired draws and campaign.json (no fitting).

For ADNI (1,025) and NACC (1,886) this script copies the prepared input (prepare_adni.py,
prepare_nacc.py) to inputs/<cohort>/raw_deid.csv next to this script, adjusts all eight events
once (run_shared_panel.correct: OLS in CN participants; cognition on age, sex and education, MRI
volumes on age, sex and ICV) and draws the 200 pairs of groups of Fig. 4e from the e3/e3
participants: numpy.random.default_rng(SeedSequence([20260920, cohort, 9001, replicate])),
cohort 1 = ADNI, 2 = NACC. The group compositions are those of the ADNI 12-event population
(69/5/12 and 56/58/86) for ADNI, with counts scaled by 93/98 (largest-remainder rounding), and
those of the NACC amyloid PET population (73/34/10 and 133/78/79) for NACC. It then writes
campaign.json with the seeds (fits of draw i: 2026092000 + 10000 * cohort + i; APOE fits:
2026099000 + cohort) and the SHA-256 of every input file, of the vendored concord package
(vendor/concord, version 0.1.1 as run) and of run_shared_panel.py.

inputs/ and campaign.json are written next to this script (both ignored by git; inputs/ holds
participant-level data). The script refuses to run if campaign.json already exists.

  python freeze_campaign.py --prepared-dir $CONCORD_WORK_DIR/mri_cognition
"""
import argparse
import os
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import shutil
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
VENDOR = HERE / 'vendor/concord'
from run_shared_panel import correct, sha, save, EVENTS, PANELS, GROUPS

def hamilton(old, factor):
    targets = [factor * x for x in old]
    result = [v.numerator // v.denominator for v in targets]
    total = (factor * sum(old)).numerator // (factor * sum(old)).denominator
    order = sorted(range(len(old)), key=lambda k: (-(targets[k] - result[k]), k))
    for k in order[:total - sum(result)]:
        result[k] += 1
    return result

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--prepared-dir', type=Path, default=None,
                    help='folder with adni/ and nacc/ written by prepare_adni.py and prepare_nacc.py '
                         '(default: $CONCORD_WORK_DIR/mri_cognition)')
    args = ap.parse_args(argv)
    if args.prepared_dir is None and os.environ.get('CONCORD_WORK_DIR'):
        args.prepared_dir = Path(os.environ['CONCORD_WORK_DIR']) / 'mri_cognition'
    if args.prepared_dir is None:
        ap.error('--prepared-dir is required (or set CONCORD_WORK_DIR)')
    config_path = HERE / 'campaign.json'
    if config_path.exists():
        raise FileExistsError('Campaign already frozen; explicit amendment required')
    config = dict(name='nc_shared_panel_20260920', version=1, R=200, B=599, alpha=.05,
                  stage='additional analysis of earlier populations, not specified in advance',
                  panels=PANELS, actual_APOE_panel='MRIcog8', consensus='repaired',
                  stop_when_decided=False, stability_resamples=0, cohorts={}, frozen_files={})
    source_files = {'adni':'raw_input.csv','nacc':'nacc_common_raw_deid.csv'}
    expected = {'adni':1025, 'nacc':1886}
    old = {'adni':([69,5,12],[56,58,86]), 'nacc':([73,34,10],[133,78,79])}
    for code, cohort in enumerate(['adni','nacc'], 1):
        origin = args.prepared_dir / cohort
        target = HERE / 'inputs' / cohort
        target.mkdir(parents=True, exist_ok=True)
        rawpath = target / 'raw_deid.csv'
        shutil.copy2(origin / source_files[cohort], rawpath)
        shutil.copy2(origin / 'manifest.json', target / 'preparation_manifest.json')
        raw = pd.read_csv(rawpath, float_precision='round_trip')
        assert len(raw) == expected[cohort]
        labels = ['CN', 'MCIc' if cohort == 'adni' else 'MCI', 'AD']
        assert set(raw.Diagnosis) == set(labels) and set(raw.APOE) == set(GROUPS)
        adjusted, correction = correct(raw)
        adjusted_path = target / 'adjusted_deid.csv'
        adjusted.to_csv(adjusted_path, index=False, float_format='%.17g')
        loaded = pd.read_csv(adjusted_path, float_precision='round_trip')
        assert np.array_equal(adjusted[EVENTS].to_numpy(),loaded[EVENTS].to_numpy())
        save(target/'correction.json', dict(formula='y - (X - mean_full_cohort(X)) @ beta_CN_slopes',
             applies_to='entire complete8 frame before all draws; K4 is a column subset', events=correction))
        pools = [np.flatnonzero((adjusted.APOE.eq('e33') & adjusted.Diagnosis.eq(d)).to_numpy()) for d in labels]
        aold, bold = old[cohort]
        factor = min([Fraction(1)] + [Fraction(len(pool),a+b) for pool,a,b in zip(pools,aold,bold)])
        a,b = hamilton(aold,factor),hamilton(bold,factor)
        assert all(x>0 for x in a+b) and all(x+y<=len(pool) for x,y,pool in zip(a,b,pools))
        all_a,all_b = [],[]
        for rep in range(config['R']):
            rng = np.random.default_rng(np.random.SeedSequence([20260920,code,9001,rep]))
            aa,bb = [],[]
            for pool,na,nb in zip(pools,a,b):
                picks = rng.permutation(pool)[:na+nb]
                aa.extend(picks[:na]);bb.extend(picks[na:])
            assert not set(aa)&set(bb)
            all_a.append(aa);all_b.append(bb)
        drawpath = target / 'n2_draws.npz'
        np.savez_compressed(drawpath,a=np.asarray(all_a,dtype=np.int64),b=np.asarray(all_b,dtype=np.int64))
        counts = pd.crosstab(adjusted.APOE,adjusted.Diagnosis).reindex(index=GROUPS,columns=labels)
        c = dict(n=len(raw),labels=labels,counts=counts.to_dict(orient='index'),
                 raw_file=str(rawpath.relative_to(HERE)),adjusted_file=str(adjusted_path.relative_to(HERE)),
                 adjusted_sha256=sha(adjusted_path),draw_file=str(drawpath.relative_to(HERE)),
                 n2_counts=dict(a=a,b=b),original_n2_counts=dict(a=aold,b=bold),
                 capacity_factor=str(factor), allocation='floor scaled arm total; Hamilton by fractional remainder; ties diagnosis order',
                 source_pool_counts=[len(x) for x in pools], n2_draw_seed_components=[20260920,code,9001,'replicate'],
                 n2_fit_seed_base=2026092000+code*10000, apoe_seed=2026099000+code,
                 correction_file=str((target/'correction.json').relative_to(HERE)))
        config['cohorts'][cohort]=c
    for path in sorted((HERE/'inputs').rglob('*')):
        if path.is_file():
            path.chmod(0o600)
            config['frozen_files'][str(path.relative_to(HERE))] = sha(path)
    for path in sorted(VENDOR.glob('*.py')):
        config['frozen_files'][str(path.relative_to(HERE))] = sha(path)
    config['frozen_files']['run_shared_panel.py'] = sha(HERE/'run_shared_panel.py')
    save(config_path,config)
    print(json.dumps({k:{name:v[name] for name in ['n','counts','source_pool_counts','n2_counts','capacity_factor']} for k,v in config['cohorts'].items()},indent=2))

if __name__=='__main__':
    main()
