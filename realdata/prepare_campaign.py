"""Six-event panel, step 3: covariate adjustment, common proportions, draws and campaign.json.

No model is fitted. For each cohort this script
  * adjusts every event once, by OLS in all CN participants of the cohort (CSF events on age and
    sex; cognitive scores on age, sex and education): y - (z - mean z) . slopes;
  * computes the common CN/MCI/AD proportions: the smallest proportion of each diagnosis over the
    six cohort-by-genotype groups, rescaled to sum to one (48.9/28.6/22.5%);
  * fixes the group size of the controlled composition experiment (160; Fig. 4a) and draws its 200
    rosters from the e3/e3 participants with numpy.random.default_rng(SeedSequence([2026092300 +
    cohort, replicate])), cohort 1 = ADNI, 2 = NACC; the same 320 participants form the matched
    (80/48/32 vs 80/48/32) and the different (112/32/16 vs 48/64/48) split;
  * writes <runtime>/inputs/<cohort>/ (participant-level; ignored by git) and <runtime>/campaign.json.

The asserted genotype-by-diagnosis counts are those of Table 2. The published runtime/campaign.json
is the file that was used. With the same data and software versions this script writes it again
unchanged; it refuses to overwrite a frozen file with different content.

  python prepare_campaign.py --work-dir $CONCORD_WORK_DIR

reads <work-dir>/audit/adni/raw_deid.csv and <work-dir>/audit/nacc/prepared/nacc_csfcog6_raw_deid.csv
(from cohort/adni_export_core6.py and cohort/nacc_prepare_core6.py) unless --adni-input and
--nacc-input are given, and writes preparation_check.json to <work-dir>.
"""
import argparse
import os
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = None   # set in main()
RUNTIME = HERE / 'runtime'
EVENTS = ['ABETA', 'TAU', 'PTAU', 'MEM', 'EXF', 'LAN']
GROUPS = ['e2', 'e33', 'e4']
DX = ['CN', 'MCI', 'AD']
EXPECTED = {
    'adni': [[50, 39, 8], [222, 273, 65], [98, 299, 149]],
    'nacc': [[106, 15, 22], [480, 104, 175], [281, 94, 307]],
}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()
    if path.exists() and path.read_bytes() != payload:
        raise RuntimeError(f'Refusing changed frozen content: {path}')
    path.write_bytes(payload)

def main(argv=None):
    global RESULTS, RUNTIME
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--work-dir', type=Path, default=os.environ.get('CONCORD_WORK_DIR'),
                    help='working folder outside the repository (default: $CONCORD_WORK_DIR)')
    ap.add_argument('--adni-input', type=Path, help='default: <work-dir>/audit/adni/raw_deid.csv')
    ap.add_argument('--nacc-input', type=Path, help='default: <work-dir>/audit/nacc/prepared/nacc_csfcog6_raw_deid.csv')
    ap.add_argument('--runtime-dir', type=Path, default=HERE / 'runtime',
                    help='runtime folder that receives inputs/ and campaign.json (default: realdata/runtime)')
    args = ap.parse_args(argv)
    if args.work_dir is None:
        ap.error('--work-dir is required (or set CONCORD_WORK_DIR)')
    RESULTS, RUNTIME = args.work_dir.resolve(), args.runtime_dir.resolve()
    sources = {'adni': args.adni_input or RESULTS / 'audit/adni/raw_deid.csv',
               'nacc': args.nacc_input or RESULTS / 'audit/nacc/prepared/nacc_csfcog6_raw_deid.csv'}
    frames, counts, cohorts = {}, {}, {}
    for cohort in ['adni', 'nacc']:
        source = sources[cohort]
        x = pd.read_csv(source, float_precision='round_trip')
        wanted = ['PTID', 'APOE', 'Diagnosis', 'Age', 'Sex', 'Education'] + EVENTS
        assert set(x.columns) == set(wanted), (cohort, x.columns.tolist())
        x = x[wanted]
        assert x.PTID.is_unique and set(x.APOE) == set(GROUPS)
        assert set(x.Diagnosis) == set(DX) and set(x.Sex) == {'Male', 'Female'}
        assert np.isfinite(x[EVENTS + ['Age', 'Education']].to_numpy(float)).all()
        assert x.Age.between(18, 120).all() and x.Education.between(0, 40).all()
        c = pd.crosstab(x.APOE, x.Diagnosis).reindex(index=GROUPS, columns=DX, fill_value=0)
        assert c.to_numpy().tolist() == EXPECTED[cohort], (cohort, c.to_numpy().tolist())
        frames[cohort], counts[cohort] = x, c
    proportions = np.concatenate([c.to_numpy() / c.sum(axis=1).to_numpy()[:, None] for c in counts.values()])
    reference = proportions.min(axis=0)
    reference = reference / reference.sum()
    # Same N in both cohorts; rule uses coverage only, before any fit.
    source_counts = np.stack([c.loc['e33'].to_numpy() for c in counts.values()])
    n = int(min(200, np.min(source_counts / np.array([1., .6, .4]))) // 10 * 10)
    assert n == 160
    allocation = {'matched_a': [n//2, 3*n//10, n//5], 'matched_b': [n//2, 3*n//10, n//5],
                  'imbalanced_a': [7*n//10, n//5, n//10], 'imbalanced_b': [3*n//10, 2*n//5, 3*n//10]}
    for ci, (cohort, raw) in enumerate(frames.items(), start=1):
        dest = RUNTIME / 'inputs' / cohort
        dest.mkdir(parents=True, exist_ok=True)
        raw.to_csv(dest/'raw_deid.csv', index=False, float_format='%.17g')
        x = raw.copy()
        x['SexMale'] = x.Sex.eq('Male').astype(float)
        cn = x.Diagnosis.eq('CN').to_numpy()
        correction = {'cohort': cohort, 'fit_population': 'All CN in the frozen complete-six-event cohort',
                      'group_labels_used': False, 'performed_once_before_sampling': True, 'events': {}}
        for event in EVENTS:
            factors = ['Age', 'SexMale'] + (['Education'] if event in EVENTS[3:] else [])
            z = x[factors].to_numpy(float)
            design = np.column_stack([np.ones(cn.sum()), z[cn]])
            assert np.linalg.matrix_rank(design) == len(factors)+1
            beta = np.linalg.lstsq(design, x.loc[cn, event].to_numpy(float), rcond=None)[0]
            center = z.mean(axis=0)
            x[event] = x[event].to_numpy(float) - (z-center) @ beta[1:]
            correction['events'][event] = {'factors': factors, 'CN_n': int(cn.sum()), 'intercept': float(beta[0]),
                                            'slopes': beta[1:].tolist(), 'full_cohort_center': center.tolist()}
        adjusted = x[['PTID', 'APOE', 'Diagnosis'] + EVENTS]
        adjusted.to_csv(dest/'adjusted_deid.csv', index=False, float_format='%.17g')
        reloaded = pd.read_csv(dest/'adjusted_deid.csv', float_precision='round_trip')
        pd.testing.assert_frame_equal(adjusted, reloaded, check_dtype=False, check_exact=True)
        write(dest/'correction.json', correction)
        pools = [np.flatnonzero((x.APOE.eq('e33') & x.Diagnosis.eq(d)).to_numpy()) for d in DX]
        plans = {k: [] for k in allocation}
        roster_hashes = []
        for rep in range(200):
            rng = np.random.default_rng(np.random.SeedSequence([2026092300+ci, rep]))
            selected = [rng.permutation(p)[:allocation['matched_a'][j]+allocation['matched_b'][j]] for j,p in enumerate(pools)]
            for design in ['matched', 'imbalanced']:
                a = np.concatenate([p[:allocation[design+'_a'][j]] for j,p in enumerate(selected)])
                b = np.concatenate([p[allocation[design+'_a'][j]:] for j,p in enumerate(selected)])
                assert len(a) == len(b) == n and not set(a) & set(b)
                assert len(set(a) | set(b)) == 2*n
                plans[design+'_a'].append(a); plans[design+'_b'].append(b)
            left = np.sort(np.r_[plans['matched_a'][-1], plans['matched_b'][-1]])
            right = np.sort(np.r_[plans['imbalanced_a'][-1], plans['imbalanced_b'][-1]])
            assert np.array_equal(left, right)
            roster_hashes.append(hashlib.sha256(left.astype('<i8').tobytes()).hexdigest())
        np.savez_compressed(dest/'composition_draws.npz', **{k: np.asarray(v, dtype='<i8') for k,v in plans.items()})
        write(dest/'draw_audit.json', {'R':200, 'allocation':allocation, 'seed_components':[2026092300+ci,'replicate'],
              'same_roster_across_mixes':True, 'within_mix_disjoint':True, 'source_APOE':'e33', 'roster_sha256':roster_hashes})
        cells = counts[cohort]
        ess = {}
        for g in GROUPS:
            cg = cells.loc[g].to_numpy(float); pg = cg/cg.sum()
            w = reference/pg
            ess[g] = {'n':int(cg.sum()), 'ESS':float(cg.sum()**2/np.sum(cg*w*w)),
                      'weights_CN_MCI_AD':w.tolist(), 'cells_CN_MCI_AD':cg.astype(int).tolist()}
        write(dest/'support.json', {'counts': {g: cells.loc[g].astype(int).tolist() for g in GROUPS}, 'reference':reference.tolist(),
              'ESS':ess, 'all_event_pairs_complete':True, 'min_group_diagnosis_cell':int(cells.to_numpy().min())})
        cohorts[cohort] = {'n':len(raw), 'raw_file':f'inputs/{cohort}/raw_deid.csv',
             'adjusted_file':f'inputs/{cohort}/adjusted_deid.csv','draw_file':f'inputs/{cohort}/composition_draws.npz',
             'draw_seed_base':2026092300+ci,'fit_seed_base':2026230000+1000*ci,
             'apoe_seed':2026239000+ci,'bootstrap_seed':2026238000+ci,
             'counts':{g:cells.loc[g].astype(int).tolist() for g in GROUPS}, 'labels':DX}
    config = {'name':'nc_common_panel_20260923','version':1,'R':200,'B':599,'alpha':.05,'bootstrap_R':200,
        'n_per_group':n,'biomarkers':EVENTS,'panels':{'CSFcog6':EVENTS,'omit_ptau':[e for e in EVENTS if e!='PTAU']},
        'primary_panel':'CSFcog6','reference':reference.tolist(),'reference_rule':'Normalize per-diagnosis minimum across both cohorts and all three APOE groups; compute from frozen complete6 counts, hold fixed for all fits.',
        'diagnosis_order':DX,'group_order':GROUPS,'composition_allocation':allocation,'fit_timeout_s':900,
        'cohorts':cohorts,'frozen_files':{str(p.relative_to(RUNTIME)):sha(p) for p in sorted((RUNTIME/'inputs').rglob('*')) if p.is_file()},
        'scope':'New coverage-selected complete common CSF-cognition campaign; no historical outputs relabelled; no simulations.',
        'selection_basis':'Biomarker definitions and aggregate coverage only; no fitted orderings or test outcomes used.',
        'primary_inference':'Three APOE D-pair comparisons per cohort on CSFcog6, Bonferroni within cohort; full 599 permutations.',
        'sensitivity':'Omit p-tau on the identical participants, corrections and reference; report alongside primary, never select by P value.',
        'bootstrap':'200 resamples within APOE x diagnosis, refit shared scoring and ordering on fixed corrected inputs; reference fixed; conditional rank stability, not full preprocessing uncertainty.'}
    write(RUNTIME/'campaign.json', config)
    write(RESULTS/'preparation_check.json', {'passed':True, 'fitting_performed':False,'counts':{k:v['counts'] for k,v in cohorts.items()},
         'reference':reference.tolist(),'n_per_N2_group':n,'planned_observed_composition_fits':3200,
         'planned_APOE_fits':7192,'planned_bootstrap_fits':400,'campaign_sha256':sha(RUNTIME/'campaign.json')})
    print(json.dumps({'passed':True,'reference':reference.tolist(),'n_per_group':n,'campaign_sha256':sha(RUNTIME/'campaign.json')}))

if __name__ == '__main__':
    main()
