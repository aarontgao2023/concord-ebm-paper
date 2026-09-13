#!/usr/bin/env python3
"""Catalog saved audit certificates; never reconstruct data or read fit-log bodies."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def chunk_key(value):
    run, cell, chunk = value.get('run_id'), value.get('cell_index'), value.get('chunk')
    if not isinstance(run, str) or not re.fullmatch(r'[A-Za-z0-9_-]+', run):
        raise ValueError('Invalid run identity')
    if any(type(x) is not int or x < 0 for x in (cell, chunk)):
        raise ValueError('Invalid cell/chunk identity')
    return run, cell, chunk


def path_key(path):
    p = Path(path)
    if not re.fullmatch(r'cell_[0-9]+', p.parent.name) or not re.fullmatch(r'chunk_[0-9]+', p.name):
        raise ValueError('Invalid cell/chunk path')
    return chunk_key({'run_id': p.parent.parent.name,
                      'cell_index': int(p.parent.name.removeprefix('cell_')),
                      'chunk': int(p.name.removeprefix('chunk_'))})


def build_catalog(audit_roots, raw_root, protocol_path):
    started = utc()
    audit_roots = [Path(p).resolve() for p in audit_roots]
    raw_root, protocol_path = Path(raw_root).resolve(), Path(protocol_path).resolve()
    reads, problems, metadata_notes = {}, [], []

    def read(path):
        path = Path(path).resolve()
        if str(path) not in reads:
            data = path.read_bytes()
            reads[str(path)] = {'path': str(path), 'sha256': sha(data), 'bytes': len(data),
                                'read_utc': utc(), 'value': json.loads(data)}
        return reads[str(path)]['value']

    protocol = read(protocol_path)
    protocol_hash = reads[str(protocol_path)]['sha256']
    plans = {x['run_id']: x for x in protocol['runs']}
    json_paths = sorted({p.resolve() for root in audit_roots for p in root.rglob('*.json')})
    manifest_paths = sorted(raw_root.glob('*/cell_*/chunk_*/manifest.json'))
    certificates, summaries = {}, []
    for path in json_paths:
        try:
            value = read(path)
            if value.get('schema_version') == 'confirmation_event_audit_v1':
                chunk_key(value)
                certificates[path] = value
            elif value.get('schema_version') == 'confirmation-event-batch-v1' and 'outcomes' in value:
                summaries.append((path, value))
        except Exception as exc:
            problems.append({'path': str(path), 'error': str(exc)})

    bindings = defaultdict(list)
    batch_reports = []
    for path, summary in summaries:
        batch = {k: summary.get(k) for k in ('finished_utc', 'status', 'source_stable',
                 'source_sha256', 'changed_source_paths', 'changed_selected_manifest_paths',
                 'selection_sha256', 'scopes', 'selected_chunks')}
        batch.update(path=str(path), sha256=reads[str(path)]['sha256'], qualified_outcomes=0)
        batch_reports.append(batch)
        try:
            selection_path = path.parent/'selection.json'
            selection = read(selection_path)
            batch['selection_path'] = str(selection_path)
            batch['actual_selection_sha256'] = reads[str(selection_path)]['sha256']
            common_ok = (summary.get('source_stable') is True
                         and summary.get('changed_source_paths') == []
                         and batch['actual_selection_sha256'] == summary.get('selection_sha256')
                         and selection.get('source_sha256') == summary.get('source_sha256'))
            selected_items = [x for x in selection['items'] if x.get('selected') is True]
            selected = {path_key(x['chunk_path']): x for x in selected_items}
            if len(selected) != len(selected_items):
                raise ValueError('Duplicate selected chunk identity')
            counts = Counter(path_key(x['chunk_path']) for x in summary['outcomes'])
            for outcome in summary['outcomes']:
                if not outcome.get('audit_json'):
                    continue
                name = outcome['audit_json']
                if not isinstance(name, str) or Path(name).name != name:
                    raise ValueError('Audit outcome must reference a local filename')
                cert_path = (path.parent/name).resolve()
                reasons = []
                if not common_ok:
                    reasons.append('summary_source_or_selection_not_verified')
                key = path_key(outcome['chunk_path'])
                cert = certificates.get(cert_path)
                item = selected.get(key)
                if counts[key] != 1:
                    reasons.append('duplicate_outcome_identity')
                if not cert or reads[str(cert_path)]['sha256'] != outcome.get('audit_json_sha256'):
                    reasons.append('certificate_missing_or_hash_mismatch')
                if (outcome.get('selected') is not True or outcome.get('checkpoint_certified') is not True
                        or outcome.get('manifest_matches_fixed_selection_at_finalize') is not True):
                    reasons.append('outcome_not_certified')
                if not item or not cert:
                    reasons.append('selection_or_certificate_missing')
                else:
                    if (chunk_key(cert) != key or outcome.get('run_id') != key[0]
                            or item.get('run_id') != key[0]
                            or item.get('cell_index') != key[1] or item.get('chunk_index') != key[2]
                            or item.get('planned_seeds') != cert.get('planned_seeds')):
                        reasons.append('selection_identity_mismatch')
                    evidence = cert.get('batch_evidence', {})
                    if (evidence.get('selection_sha256') != summary.get('selection_sha256')
                            or evidence.get('batch_source_sha256') != summary.get('source_sha256')
                            or evidence.get('candidate_index') != item.get('candidate_index')):
                        reasons.append('certificate_batch_binding_mismatch')
                    manifest = cert.get('inputs', {}).get('manifest.json', {})
                    if item.get('inputs', {}).get('manifest.json', {}).get('sha256') != manifest.get('sha256'):
                        reasons.append('fixed_manifest_mismatch')
                bindings[cert_path].append({'summary_path': str(path), 'summary_sha256': batch['sha256'],
                    'outcome_category': outcome.get('category'), 'outcome_done': outcome.get('certified_done_chunk'),
                    'binding_verified': not reasons, 'reasons': reasons})
        except Exception as exc:
            batch['catalog_error'] = str(exc)
            problems.append({'path': str(path), 'error': str(exc)})

    entries = {}

    def entry(key):
        return entries.setdefault(key, {'run_id': key[0], 'cell_index': key[1], 'chunk': key[2],
            'in_confirmation_protocol': key[0] in plans, 'audit_versions': [], 'local_manifest': None})

    for path, cert in certificates.items():
        key = chunk_key(cert)
        reasons = []
        if path_key(cert.get('chunk_path', '')) != key:
            reasons.append('certificate_path_identity_mismatch')
        plan = plans.get(key[0], {})
        inputs = cert.get('inputs', {})
        if (cert.get('checkpoint_certified') is not True
                or cert.get('findings') != [] or cert.get('deferred_reasons') != []):
            reasons.append('individual_certificate_not_clean_checkpoint')
        if (cert.get('protocol', {}).get('sha256') != protocol_hash
                or cert.get('protocol', {}).get('snapshot_sha256') != plan.get('snapshot_sha256')
                or cert.get('provenance', {}).get('config_sha256') != plan.get('config_sha256')):
            reasons.append('frozen_protocol_identity_mismatch')
        for name in ('manifest.json', 'progress.json', 'rows.jsonl', 'fit_records.jsonl'):
            inp = inputs.get(name, {})
            if (inp.get('exists') is not True or not re.fullmatch(r'[0-9a-f]{64}', inp.get('sha256', ''))
                    or not isinstance(inp.get('before'), dict) or inp.get('before') != inp.get('after')):
                reasons.append(f'audited_input_not_stable_or_hashed:{name}')
        related = bindings.get(path, [])
        qualified = not reasons and any(x['binding_verified'] for x in related)
        version = {'path': str(path), 'sha256': reads[str(path)]['sha256'],
            'inspected_utc': cert.get('inspected_utc'), 'status': cert.get('status'),
            'nchunks': cert.get('nchunks'), 'standalone_record': 'batch_evidence' not in cert,
            'audited_fit_record_count': cert.get('raw_recomputed_progress', {}).get('fit_records'),
            'individual_checkpoint_claim_verified': not reasons,
            'batch_checkpoint_qualified': qualified,
            'batch_done_qualified': qualified and cert.get('certified_done_chunk') is True
                and any(x['binding_verified'] and x['outcome_done'] is True for x in related),
            'qualification_reasons': reasons, 'batch_bindings': related,
            'audited_inputs': {name: {k: val.get(k) for k in ('path', 'sha256', 'before', 'after')}
                               for name, val in inputs.items()},
            'audit_source_sha256': cert.get('audit_source_sha256'),
            'protocol': cert.get('protocol'), 'provenance': cert.get('provenance')}
        entry(key)['audit_versions'].append(version)
        if qualified:
            for binding in related:
                if binding['binding_verified']:
                    next(b for b in batch_reports if b['path'] == binding['summary_path'])['qualified_outcomes'] += 1

    for path in manifest_paths:
        try:
            manifest = read(path)
            key = chunk_key({'run_id': manifest['config']['run_id'], 'cell_index': manifest['cell_index'], 'chunk': manifest['chunk']})
            disk_key = path_key(path.parent)
            if key != disk_key:
                if key[0] in plans or disk_key[0] in plans:
                    raise ValueError('Confirmation manifest path and content identities differ')
                metadata_notes.append({'path': str(path), 'disk_identity': list(disk_key),
                    'declared_identity': list(key), 'classification': 'outside_confirmation_reference_copy'})
                key = disk_key
            entry(key)['local_manifest'] = {'path': str(path), 'sha256': reads[str(path.resolve())]['sha256'],
                'declared_run_id': manifest['config']['run_id'], 'nchunks': manifest['nchunks'],
                'read_utc': reads[str(path.resolve())]['read_utc']}
        except Exception as exc:
            problems.append({'path': str(path), 'error': str(exc)})

    for value in entries.values():
        value['audit_versions'].sort(key=lambda x: (x.get('inspected_utc') or '', x['path']))
        qualified = [x for x in value['audit_versions'] if x['batch_checkpoint_qualified']]
        value['ever_batch_checkpoint_qualified'] = bool(qualified)
        value['ever_batch_done_qualified'] = any(x['batch_done_qualified'] for x in qualified)
        value['latest_historical_batch_certificate_path'] = qualified[-1]['path'] if qualified else None
        value['latest_historical_batch_audited_fit_record_count'] = qualified[-1]['audited_fit_record_count'] if qualified else None
        local = value['local_manifest']
        matching = [x for x in qualified if local and x['audited_inputs']['manifest.json']['sha256'] == local['sha256']]
        value['local_manifest_has_matching_batch_certificate'] = bool(matching) if local else None
        value['metadata_relation'] = ('certificate_only_pending_metadata_sync' if qualified and not local else
            'same_manifest_historically_batch_certified' if matching else
            'local_manifest_differs_from_historical_certificate' if qualified and local else
            'local_manifest_without_batch_certificate' if local else 'audit_record_only_not_batch_certified')
        latest = matching[-1] if matching else qualified[-1] if qualified else None
        value['comparison_certificate_path'] = latest['path'] if latest else None
        value['local_input_presence'] = {}
        if local:
            for name in ('progress.json', 'rows.jsonl', 'fit_records.jsonl'):
                p = Path(local['path']).parent/name
                exists = p.is_file()
                value['local_input_presence'][name] = {'exists': exists,
                    'bytes': p.stat().st_size if exists else None, 'sha256_recomputed': False,
                    'matches_audited_bytes': None}
        value['current_all_audited_input_bytes_match'] = None
        value['remote_current_bytes_match'] = None

    changed = []
    for path, record in reads.items():
        try:
            if sha(Path(path).read_bytes()) != record['sha256']:
                changed.append(path)
        except OSError:
            changed.append(path)
    final_manifests = sorted(raw_root.glob('*/cell_*/chunk_*/manifest.json'))
    final_audits = sorted({p.resolve() for root in audit_roots for p in root.rglob('*.json')})
    scope_changed = manifest_paths != final_manifests or json_paths != final_audits
    per_run = []
    for run in sorted(set(plans) | {k[0] for k in entries}):
        rows = [v for k, v in entries.items() if k[0] == run]
        counts = [x['latest_historical_batch_audited_fit_record_count'] for x in rows if x['ever_batch_checkpoint_qualified']]
        per_run.append({'run_id': run, 'in_confirmation_protocol': run in plans,
            'local_manifests': sum(x['local_manifest'] is not None for x in rows),
            'local_manifest_matching_batch_certificate': sum(x['local_manifest_has_matching_batch_certificate'] is True for x in rows),
            'local_manifest_without_matching_batch_certificate': sum(x['local_manifest'] is not None and x['local_manifest_has_matching_batch_certificate'] is False for x in rows),
            'unique_historical_batch_certified_chunks': sum(x['ever_batch_checkpoint_qualified'] for x in rows),
            'unique_historical_batch_done_chunks': sum(x['ever_batch_done_qualified'] for x in rows),
            'certified_chunks_without_local_manifest': sum(x['metadata_relation'] == 'certificate_only_pending_metadata_sync' for x in rows),
            'latest_qualified_record_count': sum(counts) if all(type(x) is int and x >= 0 for x in counts) else None,
            'scientific_planned_chunks': None, 'scientific_run_completion_assessed': False})
    return {'schema_version': 'raw_event_coverage_catalog_v1', 'scan_started_utc': started, 'scan_finished_utc': utc(),
        'catalog_source_sha256': sha(Path(__file__).read_bytes()), 'protocol_sha256': protocol_hash,
        'audit_roots': list(map(str, audit_roots)), 'local_raw_root': str(raw_root),
        'status': 'fixed_scan_scope_catalogued' if not changed and not scope_changed and not problems else 'scan_requires_attention',
        'source_bytes_stable_during_scan': not changed, 'changed_source_paths': changed,
        'directory_scope_changed_during_scan': scope_changed, 'problems': problems,
        'outside_confirmation_metadata_notes': metadata_notes,
        'counts_usable_for_fixed_scan_scope': not changed and not scope_changed and not problems,
        'per_run': per_run, 'batches': batch_reports,
        'entries': [entries[k] for k in sorted(entries)],
        'input_file_catalog': [{k: v for k, v in x.items() if k != 'value'} for x in reads.values()],
        'raw_fit_log_bodies_read': False, 'rows_bodies_read': False, 'hpc_calls': 0,
        'scientific_whole_campaign_completion_assessed': False,
        'limitations': ['Counts refer only to the explicitly scanned local manifest/certificate scope, never all planned scientific chunks.',
            'Only source-stable batch summaries with individually certified, hash-bound outcomes qualify as batch passes; old standalone records remain history.',
            'Audited remote bytes are historical. This catalog does not hash current fitlogs/rows, so current byte equivalence is unknown even when a local file exists.',
            'A certificate without a local manifest can reflect later remote execution or delayed synchronization; it is not invalid for that reason.',
            'Development manifests outside the confirmation protocol are listed separately and are not covered by the confirmation audit selection.']}


def markdown(report):
    lines = ['# 原始事件审计覆盖目录', '', f"扫描完成：{report['scan_finished_utc']}；状态：{report['status']}。", '',
        '计数按(run_id, cell_index, chunk)去重；每版证书及summary/selection的路径、SHA、审计时间保留在JSON中。', '',
        '| Run | 本地manifest | 同manifest批次已审 | 本地尚无匹配批次证书 | 历史批次通过chunk | 有证书但本地无manifest |',
        '|---|---:|---:|---:|---:|---:|']
    for r in report['per_run']:
        lines.append(f"| {r['run_id']}{'' if r['in_confirmation_protocol'] else '（development／范围外）'} | {r['local_manifests']} | {r['local_manifest_matching_batch_certificate']} | {r['local_manifest_without_matching_batch_certificate']} | {r['unique_historical_batch_certified_chunks']} | {r['certified_chunks_without_local_manifest']} |")
    lines += ['', '“历史批次通过”要求batch summary声明source稳定、selection与证书SHA绑定吻合，而且该outcome与证书均通过。单chunk早版不增加批次通过数。', '',
        '“同manifest已审”仅表示清单身份相同；本目录未重新读取fitlog或rows正文，不能声称本地或远端当前全部输入字节仍等于审计时。远端仍可能变化。', '',
        '本地缺manifest的晚近证书单独列为待同步；本地现有manifest不是全科学计划分母。本目录不判断FWER、功效、科学cohort或全模拟完成。', '',
        f"扫描期间源变化：{len(report['changed_source_paths'])}；目录范围变化：{report['directory_scope_changed_during_scan']}；问题：{len(report['problems'])}。若counts_usable_for_fixed_scan_scope为false，不使用本表作最终覆盖计数。", '']
    return '\n'.join(lines)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[2])
    p.add_argument('--output-dir', type=Path, required=True)
    args = p.parse_args(argv)
    root = args.project.resolve()
    audit_roots = [root/'results/v2_hpc_diagnostics/raw_events', root/'results/v2_hpc_diagnostics/raw_event_batches',
                   root/'results/v2_confirmation/mechanism_raw_event_batch_20260908']
    raw = root/'runs/v2/hpc_results'
    output = args.output_dir.resolve()
    if output.exists() or not output.is_relative_to(root/'results') or any(output.is_relative_to(x) for x in [raw, *audit_roots]):
        p.error('Output must be a new results directory outside all input roots')
    report = build_catalog(audit_roots, raw, root/'runs/v2/confirmation_protocol.json')
    output.mkdir(parents=True, exist_ok=False)
    (output/'coverage_catalog.json').write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    (output/'coverage_catalog.md').write_text(markdown(report))
    print(json.dumps({'status': report['status'], 'output': str(output), 'problems': len(report['problems'])}))


if __name__ == '__main__':
    main()
