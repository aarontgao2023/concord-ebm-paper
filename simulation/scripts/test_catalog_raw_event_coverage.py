import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


spec = importlib.util.spec_from_file_location('coverage_catalog', Path(__file__).with_name('catalog_raw_event_coverage.py'))
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)


class CoverageCatalogTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.audits, self.raw = self.root/'audits', self.root/'raw'
        self.audits.mkdir()
        self.raw.mkdir()
        self.protocol = self.root/'protocol.json'
        self.plan = {'run_id': 'confirm_toy', 'snapshot_sha256': 'a'*64, 'config_sha256': 'b'*64}
        self.write(self.protocol, {'runs': [self.plan]})
        self.local = self.raw/'confirm_toy/cell_0/chunk_0'
        self.manifest = {'config': {'run_id': 'confirm_toy'}, 'cell_index': 0, 'chunk': 0, 'nchunks': 2}
        self.write(self.local/'manifest.json', self.manifest)
        self.manifest_hash = M.sha((self.local/'manifest.json').read_bytes())
        self.batch = self.audits/'batch1'
        self.batch.mkdir()
        remote = '/remote/confirm_toy/cell_0/chunk_0'
        self.sources = {'/remote/auditor.py': 'c'*64}
        item = {'run_id': 'confirm_toy', 'cell_index': 0, 'chunk_index': 0, 'chunk_path': remote,
                'selected': True, 'candidate_index': 0, 'planned_seeds': [1],
                'inputs': {'manifest.json': {'sha256': self.manifest_hash}}}
        self.selection = {'items': [item], 'source_sha256': self.sources}
        self.write(self.batch/'selection.json', self.selection)
        selection_hash = M.sha((self.batch/'selection.json').read_bytes())
        self.cert = {'schema_version': 'confirmation_event_audit_v1', 'run_id': 'confirm_toy',
            'cell_index': 0, 'chunk': 0, 'chunk_path': remote, 'nchunks': 2, 'planned_seeds': [1],
            'inspected_utc': '2026-09-08T18:00:00Z', 'status': 'certified_done',
            'checkpoint_certified': True, 'certified_done_chunk': True, 'findings': [], 'deferred_reasons': [],
            'protocol': {'sha256': M.sha(self.protocol.read_bytes()), 'snapshot_sha256': 'a'*64},
            'provenance': {'config_sha256': 'b'*64},
            'inputs': {name: {'exists': True, 'sha256': self.manifest_hash if name == 'manifest.json' else 'd'*64,
                'before': {'bytes': 1}, 'after': {'bytes': 1}}
                for name in ('manifest.json', 'progress.json', 'rows.jsonl', 'fit_records.jsonl')},
            'batch_evidence': {'selection_sha256': selection_hash, 'candidate_index': 0,
                               'batch_source_sha256': self.sources}}
        self.cert_path = self.batch/'cert.json'
        self.write(self.cert_path, self.cert)
        self.summary = {'schema_version': 'confirmation-event-batch-v1', 'source_stable': True,
            'changed_source_paths': [], 'source_sha256': self.sources, 'selection_sha256': selection_hash,
            'finished_utc': '2026-09-08T18:01:00Z', 'outcomes': [{'run_id': 'confirm_toy', 'chunk_path': remote,
                'selected': True, 'checkpoint_certified': True, 'certified_done_chunk': True,
                'manifest_matches_fixed_selection_at_finalize': True, 'audit_json': 'cert.json',
                'audit_json_sha256': M.sha(self.cert_path.read_bytes()), 'category': 'certified_done'}]}
        self.write(self.batch/'summary.json', self.summary)

    def write(self, path, obj):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(obj))

    def report(self):
        return M.build_catalog([self.audits], self.raw, self.protocol)

    def save_cert(self):
        self.write(self.cert_path, self.cert)
        self.summary['outcomes'][0]['audit_json_sha256'] = M.sha(self.cert_path.read_bytes())
        self.write(self.batch/'summary.json', self.summary)

    def test_deduplicate_history_without_reading_raw_bodies(self):
        old = copy.deepcopy(self.cert)
        old.pop('batch_evidence')
        old['inspected_utc'] = '2026-09-08T17:00:00Z'
        self.write(self.audits/'early.json', old)
        (self.local/'fit_records.jsonl').write_text('DO NOT PARSE OR HASH THIS RAW BODY')
        original = Path.read_bytes
        def guarded(path):
            if path.name in ('fit_records.jsonl', 'rows.jsonl'):
                raise AssertionError('raw body read')
            return original(path)
        with patch.object(Path, 'read_bytes', guarded):
            r = self.report()
        self.assertTrue(r['counts_usable_for_fixed_scan_scope'])
        self.assertEqual(r['per_run'][0]['unique_historical_batch_certified_chunks'], 1)
        e = r['entries'][0]
        self.assertEqual(len(e['audit_versions']), 2)
        self.assertTrue(e['local_manifest_has_matching_batch_certificate'])
        self.assertIsNone(e['current_all_audited_input_bytes_match'])

    def test_unstable_batch_cannot_inherit_individual_pass(self):
        self.summary['source_stable'] = False
        self.write(self.batch/'summary.json', self.summary)
        self.assertEqual(self.report()['per_run'][0]['unique_historical_batch_certified_chunks'], 0)

    def test_certificate_hash_mismatch_is_not_qualified(self):
        self.summary['outcomes'][0]['audit_json_sha256'] = '0'*64
        self.write(self.batch/'summary.json', self.summary)
        self.assertFalse(self.report()['entries'][0]['ever_batch_checkpoint_qualified'])

    def test_selection_hash_mismatch_is_not_qualified(self):
        self.selection['extra'] = 1
        self.write(self.batch/'selection.json', self.selection)
        self.assertFalse(self.report()['entries'][0]['ever_batch_checkpoint_qualified'])

    def test_unapproved_outcome_or_changed_manifest_blocks_qualification(self):
        for field in ('checkpoint_certified', 'manifest_matches_fixed_selection_at_finalize'):
            with self.subTest(field=field):
                self.summary['outcomes'][0][field] = False
                self.write(self.batch/'summary.json', self.summary)
                self.assertFalse(self.report()['entries'][0]['ever_batch_checkpoint_qualified'])
                self.summary['outcomes'][0][field] = True

    def test_certificate_without_local_metadata_is_not_invalid(self):
        (self.local/'manifest.json').unlink()
        r = self.report()
        self.assertEqual(r['entries'][0]['metadata_relation'], 'certificate_only_pending_metadata_sync')
        self.assertEqual(r['per_run'][0]['certified_chunks_without_local_manifest'], 1)
        self.assertEqual(r['problems'], [])

    def test_changed_local_manifest_does_not_revoke_historical_pass(self):
        self.manifest['created_utc'] = 'later'
        self.write(self.local/'manifest.json', self.manifest)
        e = self.report()['entries'][0]
        self.assertTrue(e['ever_batch_checkpoint_qualified'])
        self.assertFalse(e['local_manifest_has_matching_batch_certificate'])
        self.assertEqual(e['metadata_relation'], 'local_manifest_differs_from_historical_certificate')

    def test_other_runs_are_outside_confirmation_scope(self):
        self.write(self.raw/'dev_toy/cell_0/chunk_0/manifest.json',
                   {'config': {'run_id': 'dev_toy'}, 'cell_index': 0, 'chunk': 0, 'nchunks': 1})
        r = self.report()
        dev = next(x for x in r['per_run'] if x['run_id'] == 'dev_toy')
        self.assertFalse(dev['in_confirmation_protocol'])
        self.assertEqual(dev['local_manifest_without_matching_batch_certificate'], 1)
        self.assertFalse(r['scientific_whole_campaign_completion_assessed'])

    def test_development_reference_alias_is_kept_outside_confirmation(self):
        self.write(self.raw/'dev_reference_copy/cell_0/chunk_0/manifest.json',
                   {'config': {'run_id': 'dev_toy'}, 'cell_index': 0, 'chunk': 0, 'nchunks': 1})
        r = self.report()
        self.assertEqual(r['problems'], [])
        self.assertTrue(r['counts_usable_for_fixed_scan_scope'])
        self.assertEqual(len(r['outside_confirmation_metadata_notes']), 1)
        alias = next(x for x in r['per_run'] if x['run_id'] == 'dev_reference_copy')
        self.assertEqual(alias['local_manifests'], 1)
        self.assertFalse(alias['in_confirmation_protocol'])

    def test_duplicate_selection_is_not_counted(self):
        self.selection['items'].append(copy.deepcopy(self.selection['items'][0]))
        self.write(self.batch/'selection.json', self.selection)
        selection_hash = M.sha((self.batch/'selection.json').read_bytes())
        self.summary['selection_sha256'] = selection_hash
        self.cert['batch_evidence']['selection_sha256'] = selection_hash
        self.save_cert()
        r = self.report()
        self.assertFalse(r['entries'][0]['ever_batch_checkpoint_qualified'])
        self.assertFalse(r['counts_usable_for_fixed_scan_scope'])

    def test_wrong_frozen_protocol_cannot_count(self):
        self.cert['protocol']['sha256'] = 'e'*64
        self.save_cert()
        self.assertFalse(self.report()['entries'][0]['ever_batch_checkpoint_qualified'])

    def test_duplicate_outcome_does_not_double_count_or_qualify(self):
        self.summary['outcomes'].append(copy.deepcopy(self.summary['outcomes'][0]))
        self.write(self.batch/'summary.json', self.summary)
        self.assertFalse(self.report()['entries'][0]['ever_batch_checkpoint_qualified'])

    def test_done_requires_both_outcome_and_certificate(self):
        self.summary['outcomes'][0]['certified_done_chunk'] = False
        self.write(self.batch/'summary.json', self.summary)
        e = self.report()['entries'][0]
        self.assertTrue(e['ever_batch_checkpoint_qualified'])
        self.assertFalse(e['ever_batch_done_qualified'])


if __name__ == '__main__':
    unittest.main()
