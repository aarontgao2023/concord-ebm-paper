"""Shared, fit-free I/O and deterministic draw/operator definitions."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

LABELS = ('CN', 'MCI', 'AD')
GROUPS = ('e2', 'e33', 'e4')
ENGINES = ('likelihood_ebm', 'separate', 'shared', 'concord')
DESIGNS = ('matched', 'different')
PAIR_CODES = ((0, 1), (0, 2), (1, 2))


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def file_sha(path):
    return sha_bytes(Path(path).read_bytes())


def canonical(value):
    return sha_bytes(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def frame_sha(frame):
    return sha_bytes(frame.to_csv(index=False, float_format='%.17g').encode())


def label_sha(labels):
    return sha_bytes(np.asarray(labels, dtype='<i8').tobytes())


def kendall(left, right):
    k = len(left)
    if k < 2 or sorted(left) != list(range(k)) or sorted(right) != list(range(k)):
        raise ValueError('Invalid or mismatched event permutation')
    a, b = np.argsort(left), np.argsort(right)
    return float(sum((a[i]-a[j])*(b[i]-b[j]) < 0 for i in range(k) for j in range(i+1,k)) / (k*(k-1)/2))


def reference_check(reference):
    p = np.asarray(reference, dtype=float)
    if p.shape != (3,) or not np.isfinite(p).all() or (p <= 0).any() or not np.isclose(p.sum(), 1, rtol=0, atol=1e-12):
        raise ValueError('Reference must contain three positive CN/MCI/AD proportions summing to one')
    return p


def read_campaign(path):
    path = Path(path).resolve()
    cfg = json.loads(path.read_text())
    if cfg['R'] != 200 or cfg['B'] != 599 or cfg['alpha'] != .05:
        raise ValueError('Expected the authorized R200, B599, alpha=.05 campaign')
    if cfg['n_per_group'] % 10 or not 10 <= cfg['n_per_group'] <= 200:
        raise ValueError('Common N must be a positive multiple of 10, at most 200')
    if len(cfg['biomarkers']) < 2 or len(set(cfg['biomarkers'])) != len(cfg['biomarkers']):
        raise ValueError('Invalid biomarker panel')
    reference_check(cfg['reference'])
    for relative, expected in cfg['frozen_files'].items():
        if file_sha(path.parent / relative) != expected:
            raise ValueError(f'Frozen input changed: {relative}')
    return cfg, path.parent


def read_frame(cfg, base, cohort):
    frame = pd.read_csv(base / cfg['cohorts'][cohort]['adjusted_file'], float_precision='round_trip', dtype={'PTID': str})
    columns = ['PTID', 'APOE', 'Diagnosis'] + cfg['biomarkers']
    frame = frame[columns].copy()
    if not frame.PTID.is_unique or frame[columns].isna().any().any():
        raise ValueError('Unique complete participant rows are required')
    if not np.isfinite(frame[cfg['biomarkers']].to_numpy(float)).all():
        raise ValueError('Nonfinite biomarker')
    if set(frame.APOE) != set(GROUPS) or set(frame.Diagnosis) != set(LABELS):
        raise ValueError('Expected all three APOE and diagnosis groups')
    table = pd.crosstab(frame.APOE, frame.Diagnosis).reindex(index=GROUPS, columns=LABELS, fill_value=0)
    if (table.to_numpy() == 0).any():
        raise ValueError('Every cohort by APOE by diagnosis cell needs support')
    return frame


def draw_counts(n):
    return {
        'matched': (np.array([5, 3, 2])*n//10, np.array([5, 3, 2])*n//10),
        'different': (np.array([7, 2, 1])*n//10, np.array([3, 4, 3])*n//10),
    }


def selected_indices(frame, n, seed):
    """One roster per diagnosis; both designs use exactly the same participants."""
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), 20260923, 71000000]))
    required = np.array([10, 6, 4])*n//10
    rosters = []
    for diagnosis, count in zip(LABELS, required):
        available = np.flatnonzero(frame.APOE.eq('e33') & frame.Diagnosis.eq(diagnosis))
        if len(available) < count:
            raise ValueError(f'Insufficient e33/{diagnosis}: need {count}, have {len(available)}')
        rosters.append(rng.permutation(available)[:count])
    return rosters


def n2_frame(frame, cfg, cohort, replicate, design, base=None):
    seed = cfg['cohorts'][cohort]['draw_seed_base'] + replicate
    entry = cfg['cohorts'][cohort]
    if 'draw_file' in entry:
        if base is None:
            raise ValueError('A base path is required for frozen draw files')
        with np.load(base / entry['draw_file'], allow_pickle=False) as plans:
            prefix = 'matched' if design == 'matched' else 'imbalanced'
            a, b = plans[prefix+'_a'][replicate], plans[prefix+'_b'][replicate]
            other = 'imbalanced' if design == 'matched' else 'matched'
            other_rows = np.r_[plans[other+'_a'][replicate], plans[other+'_b'][replicate]]
        rows = np.r_[a,b]
        if not np.array_equal(np.sort(rows), np.sort(other_rows)):
            raise ValueError('The two composition designs do not use the same selected participants')
        if len(a) != cfg['n_per_group'] or len(b) != cfg['n_per_group'] or len(set(rows.tolist())) != len(rows):
            raise ValueError('N2 arms overlap or have wrong size')
        if not frame.iloc[rows].APOE.eq('e33').all():
            raise ValueError('N2 source is not exclusively e33')
        a_counts,b_counts = draw_counts(cfg['n_per_group'])[design]
        for part, counts in ((a,a_counts),(b,b_counts)):
            if frame.iloc[part].Diagnosis.value_counts().reindex(LABELS,fill_value=0).tolist() != counts.tolist():
                raise ValueError('Frozen N2 count mismatch')
        rows = np.sort(rows)
        data = frame.iloc[rows].copy().reset_index(drop=True)
        data['APOE'] = np.isin(rows,b).astype(int)
        roster_hash = sha_bytes(np.sort(rows).astype('<i8').tobytes() + frame_sha(frame).encode())
        return data, {'selected_indices':rows.tolist(),'selected_roster_sha256':roster_hash,
                      'draw_seed_components':[entry['draw_seed_base'],replicate],
                      'rng_scheme':'numpy.default_rng(SeedSequence([draw_seed_base, replicate])); permute each CN/MCI/AD source pool in that order',
                      'group_counts':[a_counts.tolist(),b_counts.tolist()], 'labels_sha256':label_sha(data.APOE),
                      'frame_sha256':frame_sha(data),'frozen_draw_file_sha256':file_sha(base/entry['draw_file'])}
    rosters = selected_indices(frame, cfg['n_per_group'], seed)
    a_counts, b_counts = draw_counts(cfg['n_per_group'])[design]
    rows = np.concatenate(rosters)
    data = frame.iloc[rows].copy().reset_index(drop=True)
    group = np.concatenate([np.r_[np.zeros(a, dtype=int), np.ones(b, dtype=int)] for a,b in zip(a_counts,b_counts)])
    data['APOE'] = group
    if len(set(rows.tolist())) != len(rows) or np.bincount(group).tolist() != [cfg['n_per_group']]*2:
        raise ValueError('N2 draw is not disjoint with fixed equal group size')
    roster_hash = sha_bytes(np.sort(rows).astype('<i8').tobytes() + frame_sha(frame).encode())
    return data, {'selected_indices': rows.tolist(), 'selected_roster_sha256': roster_hash,
                  'draw_seed': int(seed), 'group_counts': [a_counts.tolist(), b_counts.tolist()],
                  'labels_sha256': label_sha(group), 'frame_sha256': frame_sha(data)}


def pair_shuffle(frame, seed, pair_index, permutation_id):
    a, b = PAIR_CODES[pair_index]
    name = f'diagnosis_pair_{GROUPS[a]}-{GROUPS[b]}'
    tag = int.from_bytes(hashlib.sha256(name.encode()).digest()[:4], 'little')
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), tag, int(permutation_id), 61000000]))
    result = frame.copy()
    groups, diagnosis = frame.APOE.to_numpy(), frame.Diagnosis.to_numpy()
    shuffled = groups.copy()
    for d in np.unique(diagnosis):
        ix = np.flatnonzero((diagnosis == d) & np.isin(groups, [a, b]))
        shuffled[ix] = rng.permutation(shuffled[ix])
    result['APOE'] = shuffled
    if not np.array_equal(shuffled[~np.isin(groups, [a,b])], groups[~np.isin(groups, [a,b])]):
        raise ValueError('D-pair changed the third group')
    if not pd.crosstab(result.APOE, result.Diagnosis).equals(pd.crosstab(frame.APOE, frame.Diagnosis)):
        raise ValueError('D-pair changed diagnosis margins')
    return result


def bootstrap_frame(frame, seed, biomarkers):
    rng = np.random.default_rng(np.random.SeedSequence([int(seed),20260923,62000000]))
    rows = []
    for group in range(3):
        for diagnosis in LABELS:
            ix = np.flatnonzero(frame.APOE.eq(group) & frame.Diagnosis.eq(diagnosis))
            if not len(ix):
                raise ValueError('Empty bootstrap stratum')
            rows.extend(rng.choice(ix, size=len(ix), replace=True).tolist())
    data = frame.iloc[rows][['PTID','APOE','Diagnosis']+list(biomarkers)].copy().reset_index(drop=True)
    data['PTID'] = [f'boot_{i:07d}' for i in range(len(data))]
    return data, {'bootstrap_indices':rows,'bootstrap_indices_sha256':sha_bytes(np.asarray(rows,dtype='<i8').tobytes()),
                  'bootstrap_seed':int(seed),'strata':'APOE x Diagnosis','n':len(data)}
