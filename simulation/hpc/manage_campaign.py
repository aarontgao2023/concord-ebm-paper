#!/usr/bin/env python3
"""Plan frozen-run continuation without submitting, cancelling, or changing jobs.

The campaign file registers immutable snapshots and chunk layouts. A fresh,
complete normalized queue snapshot and a submission-intent ledger prevent an
absent queue row from being mistaken for permission to submit again. Outputs are
JSON/Markdown. This program has no subprocess or network submission capability.
See MANAGE_CAMPAIGN.md for the schemas and the separate dispatch transaction.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shlex

import accounted_members as A
import inventory_run as I

ACTIVE_STATES = {"RUNNING", "PENDING", "CONFIGURING", "COMPLETING", "SUSPENDED",
                 "REQUEUED", "RESIZING", "SIGNALING", "STAGE_OUT", "MIXED"}
TERMINAL_STATES = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
                   "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE"}
LEDGER_OPEN = {"intent", "submission_unknown", "submitted"}
LEDGER_CLOSED = {"terminal_verified", "not_submitted_verified"}

# This is an execution-only exception for eight already audited immutable runs.
# A new scientific snapshot requires a new review, not a campaign-file override.
SINGLE_ENGINE_1H_POLICY = "frozen_single_engine_signal_1h_v1"
SINGLE_ENGINE_1H_PINS = {
    "confirm_stage_a": (
        "64dd4730cbcda4049bc956f5d6208fa9fec8f394b2e701fe33718e2b8677dda1",
        "d5b1d08322b8fe6e5ace3b0abaf9780354744e400bfa3266c2389bb2d4a81d92", 1000),
    "confirm_power_global_a": (
        "590d47752d16117ff545e0ad1c75f4ac4912dfe21061faf49392356d4f349aee",
        "8bc7fb417bc94e7e5c694122abb71d9e3ba610822c969a9b1b3664585be1eedb", 250),
    "confirm_pair_complete_ref_a": (
        "e5d909b7707d32f2e7509838670eb8a41ba4e146ed004b4499d31b44350886f3",
        "c8f97c11e0ffee91d27445ea9b5de7bd74184fc76d76dec46fa01113f3db809f", 500),
    "confirm_pair_power_k27_a": (
        "a9297e4fe381abc324bcb80ab02c53e9ec1751faf133c1095da686d92a9c8113",
        "696fdc7592a7a22ccc9999a95fc5065352229d3ba75612d334d128e53b43c1f6", 250),
    "confirm_pair_partial_e2_topup_a": (
        "4e813cd51fc094dcf1be35fa35d40a7c5ba7f8e0a592dd9ae66674f239393138",
        "b964bc198576e22558ca5a5d1954a45e1ee5bdf380d0e4b2a705d493e991c73b", 250),
    "confirm_pair_partial_e4_topup_a": (
        "8a1ca0c56df15472a1392d43a9df6fa986fb900715d174fa7916692523a7225b",
        "58895277b39e474b96ea0333b94fb1b09309894129a8b04e0db9328f91403f81", 250),
    "confirm_pair_partial_e33_a": (
        "57bb259d28fdd02a83d3ef9a0c383ec151320deb9743b3153ccc442d3aafe76b",
        "7ed869eae99ec985927ae215f18c692d7c57c15ddb3e72d5d6c9c03cce6a66d2", 500),
    "confirm_power_global_k14_k55_a": (
        "82985238d18c9ef324ed3f0dca124a514d4c71defd4f48e121f746b4fb587439",
        "6c59b0b2810277590869e39989464ba4066db9c7d713c7f7ed0f4652f7c9ed8c", 250),
}
SINGLE_ENGINE_1H_SOURCE_PINS = {
    "hpc/v2/run_array.sbatch": "1ff9a7ef3fdb8f2de0c4e4d00a5eb1382406bd2e1e8e819ce28b677c7c026add",
    "scripts/v2/run_v2.py": "5334c1c9b85de70af52e81f06292fc98c20bfd55dee16bfd57ea56c944d22926",
}
DEV_CALIBRATION_1H_POLICY = "frozen_dev_calibration_signal_1h_v1"
DEV_CALIBRATION_1H_PIN = (
    "e5a587f1af8959a0bbf4b6be731cc30dbb494b6fe9fb3c8ff2bfa025894e47ba",
    "c4a04bedbc872422ff5739c4dc4c3692460e14f4e3373ecbd9f857af0cbdcfdb",
)
DEV_CALIBRATION_1H_SOURCE_PINS = {
    "hpc/v2/run_array.sbatch": "a6b1b43684abbfaceb51bdb14f5aac2d4ba9e974d7095c6a59ad2bb75f1ec16b",
    "scripts/v2/run_v2.py": "b43085b2f88503cda40daab88c4b9176cb8dfac4c3c77dd4c825f76e11723f4d",
}
CORE_PAIRED_1H_POLICY = "frozen_core_paired_signal_1h_v1"
CORE_PAIRED_1H_PIN = (
    "e4701a864230fbc6faf64709defb4dec5fa5a3d07ac3625b95cc36cb06847464",
    "11b716e46bb5055f8b24a3ee5b756683508f66165fe3c361b359b477edf30a52",
)
CORE_PAIRED_1H_SOURCE_PINS = {
    "hpc/v2/run_array.sbatch": "1ff9a7ef3fdb8f2de0c4e4d00a5eb1382406bd2e1e8e819ce28b677c7c026add",
    "scripts/v2/run_v2.py": "5334c1c9b85de70af52e81f06292fc98c20bfd55dee16bfd57ea56c944d22926",
    "scripts/v2/paired_engine_v2.py": "797ab4c745544f762de025c95f1f0d4d0d58b787289f140ed48d62df0829c3f5",
}
# No campaign/config entry can choose another signal or lead time.
CORE_PAIRED_1H_SIGNAL_OPTION = "--signal=B:USR1@2100"


def utc():
    return datetime.now(timezone.utc)


def parse_utc(value):
    result = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if result.tzinfo is None:
        raise ValueError("Timestamp must include its UTC offset")
    return result.astimezone(timezone.utc)


def positive(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def read_json(path):
    value = json.loads(Path(path).read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    with temporary.open("w") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def snapshot_info(snapshot, config_name, expected_digest=None):
    """Verify all listed bytes; reject traversal and an unlisted runtime/config."""
    snapshot = Path(snapshot).resolve()
    checks = snapshot / "SHA256SUMS"
    digest = hashlib.sha256(checks.read_bytes()).hexdigest()
    if expected_digest and digest != expected_digest:
        raise ValueError("Snapshot SHA256SUMS digest differs from campaign registration")
    names = set()
    for line in checks.read_text().splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if not match:
            raise ValueError("Malformed SHA256SUMS line")
        declared, relative = match.groups()
        target = (snapshot / relative).resolve()
        if target == snapshot or not target.is_relative_to(snapshot) or relative in names:
            raise ValueError("Duplicate or unsafe SHA256SUMS path")
        names.add(relative)
        if hashlib.sha256(target.read_bytes()).hexdigest() != declared:
            raise ValueError(f"Snapshot byte mismatch: {relative}")
    config_relative = f"configs/v2/{config_name}.json"
    script_relative = "hpc/v2/run_array.sbatch"
    required = {config_relative, script_relative}
    required.update(f"scripts/v2/{name}" for name in I.RUNTIME_FILES)
    runtime = {f"scripts/v2/{name}" for name in I.RUNTIME_FILES + I.OPTIONAL_RUNTIME_FILES
               if (snapshot / "scripts/v2" / name).exists()}
    if not required | runtime <= names:
        raise ValueError("Required config, batch script or runtime file is not checksummed")
    script = (snapshot / script_relative).read_text()
    compact = re.sub(r"\s+", "", script)
    if "cell=$((index/nchunks))" in compact and "chunk=$((index%nchunks))" in compact:
        mapping = "cell_major"
    elif "cell=$((index%ncells))" in compact and "chunk=$((index/ncells))" in compact:
        mapping = "chunk_major"
    else:
        raise ValueError("Unknown frozen array mapping; do not guess indices")
    runtime_match = re.search(r"--time-budget-s\s+([0-9]+)", script)
    signal_match = re.search(r"#SBATCH\s+--signal=B:USR1@([0-9]+)", script)
    if not runtime_match or not signal_match:
        raise ValueError("Unrecognized runtime budget or stop-signal lead time")
    root_match = re.search(r"^root=(/[^\n]+)$", script, flags=re.M)
    if not root_match:
        raise ValueError("Cannot verify the batch script's output root")
    return {"sha256_manifest": digest, "mapping": mapping,
            "config_path": snapshot / config_relative,
            "script_path": snapshot / script_relative,
            "runtime_s": int(runtime_match[1]), "signal_lead_s": int(signal_match[1]),
            "remote_root": root_match[1]}


def validate_execution_walltime(frozen, config, cpus, wall, policy=None, *, nchunks=None):
    """Return canonical, identity-bound evidence for the narrow 1h exception.

    Hash pins are the authorization boundary. Text markers explain the reviewed
    control flow; they never authorize other source that happens to look similar.
    This reads source bytes only and never imports or executes the frozen runner.
    """
    positive(cpus, "cpus_per_task")
    positive(wall, "walltime_s")
    batch_size = max(cpus*2, 8)
    fit_budget = positive(config.get("fit_timeout_s", 300), "fit_timeout_s") * (2 if config.get("paired_standard") else 1)
    drain = ((batch_size+cpus-1)//cpus)*fit_budget
    if policy is None:
        safe_minimum = frozen["runtime_s"] + max(frozen["signal_lead_s"], drain) + 60
        if wall < safe_minimum:
            raise ValueError(f"Walltime {wall}s below runtime plus drain/checkpoint reserve {safe_minimum}s")
        return None
    if policy == DEV_CALIBRATION_1H_POLICY:
        return validate_dev_calibration_walltime(frozen, config, cpus, wall, nchunks)
    if policy == CORE_PAIRED_1H_POLICY:
        return validate_core_paired_walltime(frozen, config, cpus, wall, nchunks)
    if policy != SINGLE_ENGINE_1H_POLICY:
        raise ValueError("Unknown execution_walltime_policy")
    pin = SINGLE_ENGINE_1H_PINS.get(config.get("run_id"))
    if pin is None or frozen["sha256_manifest"] != pin[0]:
        raise ValueError("1h policy requires an explicitly allowlisted run and exact frozen snapshot")
    if nchunks != pin[2]:
        raise ValueError("1h policy requires the audited fixed chunk layout")
    config_bytes = Path(frozen["config_path"]).read_bytes()
    if hashlib.sha256(config_bytes).hexdigest() != pin[1] or json.loads(config_bytes) != config:
        raise ValueError("1h policy frozen configuration hash/content differs from the allowlist")
    if (wall != 3600 or cpus != 16 or config.get("engines") != ["repaired"]
            or config.get("paired_standard", False) is not False or config.get("fit_timeout_s") != 300
            or config.get("phase") != "confirmation" or frozen["runtime_s"] != 4000
            or frozen["signal_lead_s"] != 1500 or frozen["mapping"] != "chunk_major"):
        raise ValueError("1h policy requires 16 workers, one repaired engine, timeout300, runtime4000 and signal1500")
    snapshot = Path(frozen["script_path"]).parents[2]
    texts = {}
    for name, digest in SINGLE_ENGINE_1H_SOURCE_PINS.items():
        source = (snapshot/name).read_bytes()
        if hashlib.sha256(source).hexdigest() != digest:
            raise ValueError(f"1h policy frozen source hash differs: {name}")
        texts[name] = source.decode("utf-8")
    batch = texts["hpc/v2/run_array.sbatch"]
    runner = texts["scripts/v2/run_v2.py"]
    # The final logical command replaces the shell, so B:USR1 reaches Python.
    tail = ('exec "$HOME/ebmcal-env/bin/python" scripts/v2/run_v2.py '
            '--config "configs/v2/${config_name}.json" --cell-index "$cell" '
            '--chunk "$chunk" --nchunks "$nchunks" --workers "${SLURM_CPUS_PER_TASK}" '
            '--time-budget-s 4000 --output "$root/runs/$config_name/cell_${cell}/chunk_${chunk}"')
    logical_lines = batch.replace("\\\n", " ").strip().splitlines()
    if (re.sub(r"\s+", " ", logical_lines[-1]).strip() != tail
            or re.findall(r"^#SBATCH\s+--signal=(.+)$", batch, re.M) != ["B:USR1@1500"]):
        raise ValueError("1h policy requires the audited final exec and exact batch-only USR1 signal")
    markers = {
        "stop_handler": "def request_stop(signum, frame):\n    global STOP_REQUESTED\n    STOP_REQUESTED = True",
        "single_fit_alarm": "signal.alarm(config.get('fit_timeout_s', 300))",
        "unpaired_dispatch": "return fit_paired_job(job) if job[3] == '__paired__' else [fit_job(job)]",
        "no_pair_collapse": "if not config.get('paired_standard', False):\n        return jobs",
        "handler_installation": "for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGUSR1):\n            signal.signal(sig, request_stop)",
        "stop_predicate": "return STOP_REQUESTED or stop_file.exists() or time.monotonic()-started > args.time_budget_s",
        "stop_loop": "while not stopping():",
        "bounded_batch": "jobs = collapse_paired_jobs(jobs, config)[:max(args.workers*2, 8)]",
        "single_item_dispatch": "pool.imap_unordered(execute_job, jobs, chunksize=1)",
        "event_fsync": "stream.flush()\n                        os.fsync(stream.fileno())",
        "final_checkpoint": "        checkpoint(done)",
    }
    evidence = {}
    for name, marker in markers.items():
        if runner.count(marker) != 1:
            raise ValueError(f"1h policy frozen stop/checkpoint evidence differs: {name}")
        evidence[name] = {"line": runner[:runner.index(marker)].count("\n")+1,
                          "text_sha256": hashlib.sha256(marker.encode()).hexdigest()}
    if evidence["handler_installation"]["line"] >= evidence["stop_loop"]["line"]:
        raise ValueError("1h policy stop handler is not installed before the dispatch loop")
    return {"snapshot_sha256": pin[0], "config_sha256": pin[1], "nchunks": pin[2],
            "source_sha256": dict(SINGLE_ENGINE_1H_SOURCE_PINS),
            "batch_final_exec_python": True, "signal": "B:USR1@1500",
            "workers": 16, "single_fit_timeout_s": 300, "batch_max_jobs": batch_size,
            "nominal_batch_drain_s": drain, "signal_to_wall_s": 1500,
            "nominal_checkpoint_reserve_s": 1500-drain,
            "runtime_budget_s": 4000, "first_stop_timer": "slurm_signal",
            "earliest_signal_elapsed_s": 2040, "runner_evidence": evidence,
            "startup_requirement": "Handler installation must finish before the signal; startup and log replay precede it.",
            "drain_limit": "Nominal Python-alarm bound; 900s also covers planning/I/O/checkpoint overhead, not an unconditional native-code timeout guarantee."}


def validate_dev_calibration_walltime(frozen, config, cpus, wall, nchunks):
    """Separate exact old-run exception; no extension of confirmation pins.

    This development run has two engines in total, but the pinned old runner
    directly schedules one fit_job per engine, with no paired computational job.
    """
    pin = DEV_CALIBRATION_1H_PIN
    if config.get("run_id") != "dev_calibration_a" or frozen["sha256_manifest"] != pin[0]:
        raise ValueError("Development 1h policy requires the exact dev_calibration_a snapshot")
    config_bytes = Path(frozen["config_path"]).read_bytes()
    if hashlib.sha256(config_bytes).hexdigest() != pin[1] or json.loads(config_bytes) != config:
        raise ValueError("Development 1h policy frozen config hash/content differs")
    if (wall != 3600 or cpus != 16 or nchunks != 8 or frozen["mapping"] != "cell_major"
            or frozen["runtime_s"] != 11800 or frozen["signal_lead_s"] != 1500
            or config.get("phase") != "development" or config.get("base_seed") != 31000200
            or config.get("datasets") != 8 or config.get("engines") != ["original", "repaired"]
            or config.get("paired_standard", False) is not False or config.get("fit_timeout_s") != 300
            or config.get("bperm") != 599 or config.get("permutation_batch") != 8
            or config.get("cells") != [{"name": "IID_H0"}, {"name": "REF_H0"},
                                      {"name": "STAGE_H0", "overrides": {"stage_eta": 1.0}}]):
        raise ValueError("Development 1h policy requires the audited R8/cell_major/nchunks8/16worker/single-fit300/runtime11800 layout")
    snapshot = Path(frozen["script_path"]).parents[2]
    texts = {}
    for name, digest in DEV_CALIBRATION_1H_SOURCE_PINS.items():
        source = (snapshot/name).read_bytes()
        if hashlib.sha256(source).hexdigest() != digest:
            raise ValueError(f"Development 1h policy frozen source hash differs: {name}")
        texts[name] = source.decode("utf-8")
    batch, runner = texts["hpc/v2/run_array.sbatch"], texts["scripts/v2/run_v2.py"]
    tail = ('exec "$HOME/ebmcal-env/bin/python" scripts/v2/run_v2.py '
            '--config "configs/v2/${config_name}.json" --cell-index "$cell" '
            '--chunk "$chunk" --nchunks "$nchunks" --workers "${SLURM_CPUS_PER_TASK}" '
            '--time-budget-s 11800 --output "$root/runs/$config_name/cell_${cell}/chunk_${chunk}"')
    logical_lines = batch.replace("\\\n", " ").strip().splitlines()
    if (re.sub(r"\s+", " ", logical_lines[-1]).strip() != tail
            or re.findall(r"^#SBATCH\s+--signal=(.+)$", batch, re.M) != ["B:USR1@1500"]):
        raise ValueError("Development 1h policy requires final exec Python and exact batch-only USR1 signal")
    markers = {
        "stop_handler": "def request_stop(signum, frame):\n    global STOP_REQUESTED\n    STOP_REQUESTED = True",
        "single_fit_alarm": "signal.alarm(config.get('fit_timeout_s', 300))",
        "handler_installation": "for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGUSR1):\n            signal.signal(sig, request_stop)",
        "stop_predicate": "return STOP_REQUESTED or stop_file.exists() or time.monotonic()-started > args.time_budget_s",
        "stop_loop": "while not stopping():",
        "bounded_batch": "jobs = jobs[:max(args.workers*2, 8)]",
        "single_engine_item_dispatch": "pool.imap_unordered(fit_job, jobs, chunksize=1)",
        "event_fsync": "stream.flush()\n                    os.fsync(stream.fileno())",
        "final_checkpoint": "        checkpoint(done)",
    }
    evidence = {}
    for name, marker in markers.items():
        if runner.count(marker) != 1:
            raise ValueError(f"Development 1h policy frozen stop/checkpoint evidence differs: {name}")
        evidence[name] = {"line": runner[:runner.index(marker)].count("\n")+1,
                          "text_sha256": hashlib.sha256(marker.encode()).hexdigest()}
    if evidence["handler_installation"]["line"] >= evidence["stop_loop"]["line"]:
        raise ValueError("Development 1h stop handler is not installed before its loop")
    return {"snapshot_sha256": pin[0], "config_sha256": pin[1], "nchunks": 8,
            "mapping": "cell_major", "source_sha256": dict(DEV_CALIBRATION_1H_SOURCE_PINS),
            "batch_final_exec_python": True, "signal": "B:USR1@1500", "workers": 16,
            "engines_per_fit_job": 1, "configured_engines": ["original", "repaired"],
            "single_fit_timeout_s": 300, "batch_max_jobs": 32, "nominal_batch_drain_s": 600,
            "signal_to_wall_s": 1500, "nominal_checkpoint_reserve_s": 900,
            "runtime_budget_s": 11800, "first_stop_timer": "slurm_signal",
            "earliest_signal_elapsed_s": 2040, "runner_evidence": evidence,
            "startup_requirement": "Environment checks and existing event replay precede handler installation; finish before the signal.",
            "drain_limit": "Nominal Python-alarm bound; 900s covers planning/worker startup/I/O/checkpoint overhead, not an unconditional native-code timeout guarantee.",
            "scientific_scope": "Original development R8/B599 only; original/repaired remain separate single-engine jobs, with no source, seed or configuration changes."}


def validate_core_paired_walltime(frozen, config, cpus, wall, nchunks):
    """Exact core-only request; effective Slurm delivery remains unverified."""
    pin = CORE_PAIRED_1H_PIN
    if config.get("run_id") != "confirm_core_a" or frozen["sha256_manifest"] != pin[0]:
        raise ValueError("Core paired 1h policy requires the exact confirm_core_a snapshot")
    config_bytes = Path(frozen["config_path"]).read_bytes()
    if hashlib.sha256(config_bytes).hexdigest() != pin[1] or json.loads(config_bytes) != config:
        raise ValueError("Core paired 1h policy frozen config hash/content differs")
    if (wall != 3600 or cpus != 16 or nchunks != 500 or frozen["mapping"] != "chunk_major"
            or frozen["runtime_s"] != 4000 or frozen["signal_lead_s"] != 1500
            or config.get("phase") != "confirmation" or config.get("base_seed") != 42100000
            or config.get("datasets") != 1000 or config.get("engines") != ["original", "repaired"]
            or config.get("paired_standard") is not True or config.get("fit_timeout_s") != 300
            or config.get("fast_likelihood") is not True or config.get("bperm") != 599
            or config.get("permutation_batch") != 8 or config.get("full_p_first_n") != 100
            or config.get("alpha") != .05 or config.get("rule") != "le"
            or config.get("cells") != [{"name": "IID_H0"}, {"name": "REF_H0"}]
            or config.get("schemes") != {"unrestricted": {"stratify": "none", "require_max": True},
                                        "diagnosis": {"stratify": "diagnosis", "require_max": True}}):
        raise ValueError("Core paired 1h policy requires the exact R1000/B599/nchunks500/chunk_major/16worker/paired600 layout")
    snapshot = Path(frozen["script_path"]).parents[2]
    texts = {}
    for name, digest in CORE_PAIRED_1H_SOURCE_PINS.items():
        source = (snapshot/name).read_bytes()
        if hashlib.sha256(source).hexdigest() != digest:
            raise ValueError(f"Core paired 1h policy frozen source hash differs: {name}")
        texts[name] = source.decode("utf-8")
    batch, runner = texts["hpc/v2/run_array.sbatch"], texts["scripts/v2/run_v2.py"]
    tail = ('exec "$HOME/ebmcal-env/bin/python" scripts/v2/run_v2.py '
            '--config "configs/v2/${config_name}.json" --cell-index "$cell" '
            '--chunk "$chunk" --nchunks "$nchunks" --workers "${SLURM_CPUS_PER_TASK}" '
            '--time-budget-s 4000 --output "$root/runs/$config_name/cell_${cell}/chunk_${chunk}"')
    logical_lines = batch.replace("\\\n", " ").strip().splitlines()
    if (re.sub(r"\s+", " ", logical_lines[-1]).strip() != tail
            or re.findall(r"^#SBATCH\s+--signal=(.+)$", batch, re.M) != ["B:USR1@1500"]):
        raise ValueError("Core paired policy requires frozen final exec and unchanged default signal1500")
    markers = {
        "stop_handler": "def request_stop(signum, frame):\n    global STOP_REQUESTED\n    STOP_REQUESTED = True",
        "single_fit_alarm": "signal.alarm(config.get('fit_timeout_s', 300))",
        "paired_total_alarm": "total_budget = 2*config.get('fit_timeout_s', 300)\n    signal.alarm(total_budget)",
        "paired_dispatch": "return fit_paired_job(job) if job[3] == '__paired__' else [fit_job(job)]",
        "paired_completed_results": "completed = getattr(exc, 'completed_results', {})",
        "paired_restore_completed": "result = completed.get(engine)\n            if result is not None:",
        "environment_before_replay": "current_environment = verify_pyebm()",
        "startup_replay": "records = load_events(event_path)",
        "handler_installation": "for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGUSR1):\n            signal.signal(sig, request_stop)",
        "stop_predicate": "return STOP_REQUESTED or stop_file.exists() or time.monotonic()-started > args.time_budget_s",
        "stop_loop": "while not stopping():",
        "bounded_paired_batch": "jobs = collapse_paired_jobs(jobs, config)[:max(args.workers*2, 8)]",
        "one_computational_job_dispatch": "pool.imap_unordered(execute_job, jobs, chunksize=1)",
        "validate_saved_half_before_append": "replay_audits[candidate_key] = replay_consistent(records[candidate_key], candidate)",
        "preserve_saved_record": "'policy': 'keep_original_record'",
        "event_fsync": "stream.flush()\n                        os.fsync(stream.fileno())",
        "final_checkpoint": "        checkpoint(done)",
    }
    evidence = {}
    for name, marker in markers.items():
        if runner.count(marker) != 1:
            raise ValueError(f"Core paired frozen stop/checkpoint evidence differs: {name}")
        evidence[name] = {"line": runner[:runner.index(marker)].count("\n")+1,
                          "text_sha256": hashlib.sha256(marker.encode()).hexdigest()}
    if not (evidence["environment_before_replay"]["line"] < evidence["startup_replay"]["line"]
            < evidence["handler_installation"]["line"] < evidence["stop_loop"]["line"]
            < evidence["validate_saved_half_before_append"]["line"] < evidence["event_fsync"]["line"]
            < evidence["final_checkpoint"]["line"]):
        raise ValueError("Core paired startup/handler/replay/append/checkpoint ordering differs")
    paired = texts["scripts/v2/paired_engine_v2.py"]
    completed_marker = "abort.cause.completed_results = results"
    if paired.count(completed_marker) != 1:
        raise ValueError("Core paired engine must attach completed results on abort")
    return {"snapshot_sha256": pin[0], "config_sha256": pin[1], "nchunks": 500, "mapping": "chunk_major",
            "source_sha256": dict(CORE_PAIRED_1H_SOURCE_PINS), "batch_final_exec_python": True,
            "frozen_default_signal": "B:USR1@1500", "requested_signal": "B:USR1@2100",
            "requested_sbatch_option": CORE_PAIRED_1H_SIGNAL_OPTION,
            "signal_effective_verified": False, "signal_verification_status": "unverified_until_first_execution",
            "workers": 16, "configured_engines": ["original", "repaired"], "engines_per_computational_job": 2,
            "single_fit_timeout_s": 300, "paired_total_timeout_s": 600, "batch_max_jobs": 32,
            "timeout_scope": "300s is the configured single-fit baseline. The shared paired job has one total600s alarm, not two independent300s arm alarms.",
            "nominal_batch_waves": 2, "nominal_batch_drain_s": 1200, "requested_signal_to_wall_s": 2100,
            "nominal_checkpoint_reserve_s": 900, "runtime_budget_s": 4000,
            "first_stop_timer": "requested_slurm_signal_if_delivered_after_handler_installation",
            "nominal_requested_signal_elapsed_s": 1500, "earliest_requested_signal_elapsed_s": 1440,
            "signal_early_delivery_allowance_s": 60, "runner_evidence": evidence,
            "paired_engine_evidence": {"completed_results_on_abort": {"line": paired[:paired.index(completed_marker)].count("\n")+1,
                "text_sha256": hashlib.sha256(completed_marker.encode()).hexdigest()}},
            "startup_requirement": "Environment checks and full existing-event replay precede handler installation; startup must finish before elapsed1440s.",
            "drain_limit": "Nominal Python-alarm bound only. The900s reserve also covers startup/planning/IPC/I/O/cleanup; native-code delays can exceed it. No unconditional wall-clock guarantee.",
            "partial_pair_limit": "Completed results attached by a caught abort are retained; a worker result not returned and fsynced before SIGKILL is not guaranteed saved.",
            "verification_requirement": "Requested argv only; effective delivery and execution behavior remain unverified. Keep strict argv/request-parsing evidence separate from first-execution stopped/done flags, STOP absence, durable records and normal exit. No controller WarnSignal readback field is assumed available; an early done task does not validate the warning/checkpoint path."}


def array_index(cell, chunk, ncells, nchunks, mapping):
    return cell * nchunks + chunk if mapping == "cell_major" else chunk * ncells + cell


def queue_view(queue, project, now, max_age):
    if queue.get("schema_version") != "campaign_queue_v1" or queue.get("complete") is not True:
        raise ValueError("A complete campaign_queue_v1 queue snapshot is required")
    if queue.get("scope") != "all_user_jobs":
        raise ValueError("Queue snapshot must cover all user jobs, with explicit project labels")
    age = (now - parse_utc(queue["captured_utc"])).total_seconds()
    if not 0 <= age <= max_age:
        raise ValueError(f"Queue snapshot is stale or future-dated ({age:.1f} seconds)")
    jobs = queue.get("jobs")
    if not isinstance(jobs, list):
        raise ValueError("queue.jobs must be a list")
    reservations = Counter()
    submitted_members = Counter()
    active = defaultdict(list)
    ids = set()
    for job in jobs:
        job_id = job.get("job_id")
        if not isinstance(job_id, str) or not re.fullmatch(r"[0-9]+(?:_[0-9]+)?", job_id) or job_id in ids:
            raise ValueError("Queue job IDs must be unique numeric jobs/array parents")
        ids.add(job_id)
        if job.get("state") in TERMINAL_STATES:
            continue
        if job.get("state") not in ACTIVE_STATES:
            raise ValueError(f"Unknown queue state: {job.get('state')}")
        label = job.get("project")
        if not isinstance(label, str) or not label:
            raise ValueError("Every active job requires an explicit project label")
        cpus = positive(job.get("cpus_per_task"), "queue cpus_per_task")
        task_ids = job.get("array_task_ids")
        if task_ids is None:
            reservation = cpus
            members = 1
        else:
            if not isinstance(task_ids, list) or not task_ids or any(type(t) is not int or t < 0 for t in task_ids):
                raise ValueError("Array tasks must explicitly list all pending/running task indices")
            if len(set(task_ids)) != len(task_ids):
                raise ValueError("Duplicate queue array task ID")
            throttle = positive(job.get("max_concurrent"), "array max_concurrent")
            running = I.count(job.get("running_tasks"), "running_tasks")
            if running > len(task_ids):
                raise ValueError("running_tasks exceeds active array members")
            reservation = cpus * max(running, min(throttle, len(task_ids)))
            members = len(task_ids)
        reservations[label] += reservation
        submitted_members[label] += members
        if label == project and task_ids is not None and not job.get("run_id"):
            raise ValueError("Every project array requires a run_id for duplicate-target protection")
        if label == project and job.get("run_id"):
            if task_ids is None:
                raise ValueError("Registered simulation runs need array_task_ids, even for one task")
            for index in task_ids:
                active[(job["run_id"], index)].append(job)
    return {"reservations_by_project": dict(reservations),
            "active_job_members_by_project": dict(submitted_members),
            "all_user_active_job_members": sum(submitted_members.values()),
            "project_reserved_cpus": reservations[project],
            "all_user_reserved_cpus": sum(reservations.values()),
            "active": active, "jobs_by_id": {job['job_id']: job for job in jobs}, "age_s": age,
            "validated_utc": now.isoformat()}


def certified_terminal_members(entry, queue_jobs, now=None):
    """One release rule for target protection, CPU reservations and job slots.

    An uncertain submission remains wholly reserved even if it carries a job ID
    or accounting-looking fields. Only an accepted receipt can use strict sacct
    evidence, and currently active members always retain their reservation.
    """
    if entry.get("state") != "submitted":
        return set()
    visible = queue_jobs.get(entry.get("job_id"), {})
    active = (visible.get("array_task_ids") or []) if visible.get("state") in ACTIVE_STATES else []
    return A.releasable_members(entry, active_members=active, now=now or utc())


def ledger_view(ledger, queue_jobs, now=None):
    if ledger.get("schema_version") != "campaign_ledger_v1" or not isinstance(ledger.get("submissions"), list):
        raise ValueError("An explicit campaign_ledger_v1 ledger is required, even when empty")
    unresolved = defaultdict(list)
    unseen_reserved = 0
    plan_ids = set()
    for entry in ledger["submissions"]:
        state = entry.get("state")
        if state not in LEDGER_OPEN | LEDGER_CLOSED:
            raise ValueError("Unknown submission ledger state")
        if not isinstance(entry.get("run_id"), str) or not isinstance(entry.get("array_task_ids"), list):
            raise ValueError("Every ledger receipt must retain run and array task identities")
        if not I.is_hash(entry.get("snapshot_sha256")):
            raise ValueError("Every ledger receipt must retain the frozen snapshot digest")
        if not I.is_hash(entry.get("plan_id")) or entry["plan_id"] in plan_ids:
            raise ValueError("Ledger plan IDs must be unique SHA256 identities")
        plan_ids.add(entry["plan_id"])
        if state in LEDGER_CLOSED:
            if not entry.get("reconciliation_evidence"):
                raise ValueError("Closing an intent requires explicit accounting/no-submit evidence")
            continue
        cpus = positive(entry.get("cpus_per_task"), "ledger cpus_per_task")
        concurrency = positive(entry.get("max_concurrent"), "ledger max_concurrent")
        if not entry["array_task_ids"] or len(set(entry["array_task_ids"])) != len(entry["array_task_ids"]):
            raise ValueError("Open ledger intent needs unique task indices")
        visible = queue_jobs.get(entry.get("job_id"), {})
        if visible.get("state") in ACTIVE_STATES:
            if (visible.get("project") != ledger["project"] or visible.get("run_id") != entry["run_id"]
                    or visible.get("snapshot_sha256") != entry["snapshot_sha256"]
                    or not set(visible.get("array_task_ids") or []) <= set(entry["array_task_ids"])):
                raise ValueError("Visible job ID does not match its submission receipt identity")
        released = certified_terminal_members(entry, queue_jobs, now)
        remaining = set(entry["array_task_ids"]) - released
        if visible.get("state") not in ACTIVE_STATES:
            unseen_reserved += cpus * min(concurrency, len(remaining))
        # A submitted job absent from squeue remains unresolved until sacct or
        # equivalent evidence is recorded; absence alone never permits retry.
        # Strict off-queue terminal members can resume before their parent ends.
        for index in entry["array_task_ids"]:
            if type(index) is not int or index < 0:
                raise ValueError("Invalid ledger array index")
            if index in remaining:
                unresolved[(entry["run_id"], index)].append(entry)
    return unresolved, unseen_reserved


def user_job_budget(ledger, queue_state, cap=None):
    """Count submitted members, independently of array CPU throttles.

    Missing open-receipt members remain reserved unless accepted submission
    identity plus strict terminal sacct evidence certifies that member ended.
    Uncertain receipts and currently active members are never released.
    """
    if cap is not None:
        positive(cap, "max_user_submitted_jobs")
    now = parse_utc(queue_state["validated_utc"]) if queue_state.get("validated_utc") else utc()
    ledger_view(ledger, queue_state["jobs_by_id"], now)
    unseen = released_total = 0
    for entry in ledger["submissions"]:
        if entry["state"] not in LEDGER_OPEN:
            continue
        visible = queue_state["jobs_by_id"].get(entry.get("job_id"), {})
        members = set(visible.get("array_task_ids") or []) if visible.get("state") in ACTIVE_STATES else set()
        released = certified_terminal_members(entry, queue_state["jobs_by_id"], now)
        released_total += len(released)
        unseen += len(set(entry["array_task_ids"]) - members - released)
    active = queue_state["all_user_active_job_members"]
    return {"cap": cap, "queue_active_members": active,
            "queue_members_by_project": queue_state["active_job_members_by_project"],
            "unresolved_submission_reserved": unseen,
            "certified_terminal_members_released": released_total,
            "used_before_plan": active + unseen,
            "remaining_before_plan": None if cap is None else max(0, cap-active-unseen)}


def row_work(row, config, force_full):
    if row.get("status") != "ok":
        pending = row.get("status") in {"pending", "running"}
        return {"terminal": not pending, "satisfied": False, "failed_observed": int(not pending),
                "failed_permutations": 0, "states": ["observed_pending" if pending else "observed_failed"]}
    I.validate_observed(row)
    states = [I.scheme_state(row.get("schemes", {}).get(name, {
        "requested_nperm": config["bperm"], "nperm": 0, "attempted": 0,
        "complete": False, "gt": [0]*3, "eq": [0]*3,
        **({"maxgt": 0, "maxeq": 0} if definition.get("require_max") else {})
    }), definition, config, force_full) for name, definition in config.get("schemes", {}).items()]
    return {"terminal": all(s["terminal"] for s in states),
            "satisfied": all(s["obligation_satisfied"] for s in states),
            "failed_observed": 0, "failed_permutations": sum(s["failed"] for s in states),
            "states": [s["state"] for s in states]}


def chunk_work(directory, seeds, config, full_p):
    expected = {(seed, engine) for seed in seeds for engine in config["engines"]}
    rows = {}
    path = directory / "rows.jsonl"
    if path.exists():
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            key = (row["seed"], row["engine"])
            if key in rows or key not in expected:
                raise ValueError("Duplicate or unexpected chunk row; do not choose a winner")
            rows[key] = row
    states = [row_work(row, config, row["seed"] in full_p) for row in rows.values()]
    terminal = len(rows) == len(expected) and all(s["terminal"] for s in states)
    satisfied = len(rows) == len(expected) and all(s["satisfied"] for s in states)
    return {"missing_observed_rows": len(expected)-len(rows),
            "failed_observed": sum(s["failed_observed"] for s in states),
            "failed_permutations": sum(s["failed_permutations"] for s in states),
            "scheme_states": dict(Counter(v for s in states for v in s["states"])),
            "terminal": terminal, "satisfied": satisfied,
            "state": "complete" if satisfied else "terminal_with_failures" if terminal else "needs_resume"}


def make_plan(campaign, queue, ledger, now=None, weekly_used=None):
    now = now or utc()
    project = campaign.get("project", "ebmcal_v2")
    cpu_cap = positive(campaign.get("max_project_cpus", 256), "max_project_cpus")
    max_age = positive(campaign.get("queue_max_age_s", 180), "queue_max_age_s")
    runs = campaign.get("runs")
    if not isinstance(runs, list) or not runs:
        raise ValueError("campaign.runs must be a nonempty list")
    if len({run.get("run_id") for run in runs}) != len(runs):
        raise ValueError("Campaign run IDs must be unique")
    roots = [str(Path(run["output_root"]).resolve()) for run in runs]
    if len(set(roots)) != len(roots):
        raise ValueError("Registered runs must have distinct output roots")
    if ledger.get("project") != project:
        raise ValueError("Submission ledger must explicitly match the campaign project")
    q = queue_view(queue, project, now, max_age)
    intents, unknown_reserved = ledger_view(ledger, q["jobs_by_id"], now)
    slots = user_job_budget(ledger, q, campaign.get("max_user_submitted_jobs"))
    free_slots = slots["remaining_before_plan"]
    pause_usage = weekly_used is not None and weekly_used >= 90
    if weekly_used is not None and not 0 <= weekly_used <= 100:
        raise ValueError("weekly_used must be a percentage")
    report = {"schema_version": "campaign_plan_v1", "created_utc": now.isoformat(),
              "mode": "dry_run_only", "project": project,
              "planner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "inventory_sha256": hashlib.sha256(Path(I.__file__).read_bytes()).hexdigest(),
              "accounting_validator_sha256": hashlib.sha256(Path(A.__file__).read_bytes()).hexdigest(),
              "weekly_used_percent": weekly_used, "pause_new_submissions": pause_usage,
              "queue_captured_utc": queue["captured_utc"], "queue_age_s": q["age_s"],
              "cpu_budget": {"cap": cpu_cap, "project_reserved": q["project_reserved_cpus"],
                             "unresolved_submission_reserved": unknown_reserved,
                             "all_user_reserved": q["all_user_reserved_cpus"],
                             "other_projects": {k:v for k,v in q["reservations_by_project"].items() if k != project}},
              "user_job_budget": slots,
              "runs": [], "submission_plans": [], "issues": []}
    free = max(0, cpu_cap-q["project_reserved_cpus"]-unknown_reserved)
    if q["project_reserved_cpus"] + unknown_reserved > cpu_cap:
        report["issues"].append("Project reservations already exceed its cap; no new jobs planned. Other projects are not controlled.")
    if free_slots == 0:
        report["issues"].append("User submitted-job slots are exhausted; no new array members planned. All projects count toward this separate optional limit.")
    known_ids = {run["run_id"] for run in runs}
    for (run_id, index), jobs in q["active"].items():
        if run_id not in known_ids:
            report["issues"].append(f"Project run {run_id} is not registered here; its CPUs are still reserved.")
    for run in sorted(runs, key=lambda r: (r.get("priority", 100), r["run_id"])):
        result = {"run_id": run["run_id"], "chunks": [], "issues": []}
        report["runs"].append(result)
        try:
            name = run.get("config_name", run["run_id"])
            if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
                raise ValueError("config_name must be a simple filename stem")
            if name != run["run_id"]:
                raise ValueError("config_name must equal run_id; aliases cannot share an output directory")
            if not I.is_hash(run.get("snapshot_sha256")):
                raise ValueError("Register the exact SHA256SUMS digest before planning")
            frozen = snapshot_info(run["snapshot"], name, run["snapshot_sha256"])
            config = read_json(frozen["config_path"])
            cells, all_seeds, full_p = I.validate_config(config)
            if config["run_id"] != run["run_id"]:
                raise ValueError("Campaign run ID differs from frozen config")
            if config["phase"] not in {"development", "confirmation", "diagnostic"}:
                raise ValueError("Frozen phase is not supported by the runner")
            nchunks = positive(run["nchunks"], "nchunks")
            if nchunks > config["datasets"]:
                raise ValueError("Do not submit empty chunks; nchunks exceeds datasets")
            max_array_size = positive(campaign.get("max_array_size", 1000), "max_array_size")
            if len(cells)*nchunks > max_array_size:
                raise ValueError(f"Frozen mapping requires an index >= MaxArraySize ({max_array_size}); choose layout before freezing")
            if any(run_id == run["run_id"] and not 0 <= index < len(cells)*nchunks
                   for run_id, index in q["active"]):
                raise ValueError("Active queue contains an array index outside the registered run layout")
            if run.get("mapping", frozen["mapping"]) != frozen["mapping"]:
                raise ValueError("Registered mapping differs from frozen batch script")
            root = Path(run["output_root"]).resolve()
            remote_snapshot = str(Path(run.get("remote_snapshot", run["snapshot"])))
            expected_output = f"{frozen['remote_root']}/runs/{name}"
            if Path(run.get("remote_output_root", run["output_root"])).resolve() != Path(expected_output).resolve():
                raise ValueError("Remote output root differs from frozen batch script destination")
            cpus = positive(run.get("cpus_per_task", 16), "cpus_per_task")
            wall = positive(run.get("walltime_s", 14400), "walltime_s")
            policy = run.get("execution_walltime_policy")
            execution_evidence = validate_execution_walltime(frozen, config, cpus, wall, policy, nchunks=nchunks)
            inventory = I.inventory(frozen["config_path"], root,
                                    source_root=Path(run["snapshot"])/"scripts/v2", nchunks_override=nchunks)
            result.update(snapshot_sha256=frozen["sha256_manifest"], mapping=frozen["mapping"],
                          inventory_ready=inventory["ready_for_final_analysis"],
                          inventory_totals=inventory["totals"], output_root=str(root))
            errors = [issue for issue in inventory["issues"] if issue["severity"] == "error"]
            if errors or not inventory["stable_rows_snapshot"]:
                result["issues"].extend(errors or [{"code": "unstable_rows", "message": "Repeat after checkpoint settles"}])
                result["state"] = "blocked_integrity"
                continue
            stop_file = Path(config["stop_file"]) if config.get("stop_file") else None
            candidates = []
            for chunk in range(nchunks):
                seeds = sorted(all_seeds)[chunk::nchunks]
                for cell in range(len(cells)):
                    index = array_index(cell, chunk, len(cells), nchunks, frozen["mapping"])
                    directory = root/f"cell_{cell}"/f"chunk_{chunk}"
                    item = {"cell_index": cell, "cell": cells[cell], "chunk": chunk,
                            "array_index": index, "seeds": seeds, "output": str(directory)}
                    result["chunks"].append(item)
                    work = chunk_work(directory, seeds, config, full_p)
                    item.update(work)
                    active = q["active"].get((run["run_id"], index), [])
                    pending = intents.get((run["run_id"], index), [])
                    if any(job.get("snapshot_sha256") != frozen["sha256_manifest"] for job in active):
                        raise ValueError("Active job lacks matching frozen snapshot identity")
                    if any(entry["snapshot_sha256"] != frozen["sha256_manifest"] for entry in pending):
                        raise ValueError("Open submission intent belongs to a different snapshot")
                    if len(active) > 1:
                        item.update(state="duplicate_active_jobs", active_job_ids=[j["job_id"] for j in active])
                        result["issues"].append(f"Multiple jobs target array index {index}; do not dispatch another.")
                    elif active:
                        item.update(state="active", active_job_ids=[active[0]["job_id"]])
                    elif pending:
                        item.update(state="submission_reconciliation_required", intent_count=len(pending))
                    elif (stop_file is not None and stop_file.exists()) or (stop_file is None and (directory.parent/"STOP").exists()):
                        item["state"] = "stop_file_present"
                    elif work["state"] == "needs_resume":
                        item["state"] = "eligible_resume"
                        candidates.append(index)
            result["state_counts"] = dict(Counter(c["state"] for c in result["chunks"]))
            result["state"] = "complete" if all(c["state"] == "complete" for c in result["chunks"]) else "pending_or_blocked"
            if candidates and not pause_usage and free >= cpus and free_slots != 0:
                parallel = min(free//cpus, positive(run.get("max_parallel_tasks", 16), "max_parallel_tasks"), len(candidates))
                # Limit command length and job-array cardinality; array index limits
                # remain a site setting the dispatcher must verify separately.
                task_limit = positive(campaign.get("max_tasks_per_submission", 500), "max_tasks_per_submission")
                if free_slots is not None:
                    task_limit = min(task_limit, free_slots)
                selected = sorted(candidates)[:task_limit]
                parallel = min(parallel, len(selected))
                identity = {"run_id": run["run_id"], "snapshot_sha256": frozen["sha256_manifest"],
                            "array_task_ids": selected, "nchunks": nchunks, "mapping": frozen["mapping"],
                            "cpus_per_task": cpus, "max_concurrent": parallel, "walltime_s": wall,
                            "dispatch_generation": 1+sum(e["run_id"] == run["run_id"] for e in ledger["submissions"])}
                # Optional extension: ordinary historical identity bytes stay unchanged.
                if execution_evidence is not None:
                    identity.update(execution_walltime_policy=policy, execution_evidence=execution_evidence)
                plan_id = I.digest(identity)
                array = ",".join(map(str, selected))+f"%{parallel}"
                hh, rem = divmod(wall, 3600)
                mm, ss = divmod(rem, 60)
                logs = f"{frozen['remote_root']}/logs/{name}"
                argv = ["sbatch", "--parsable", f"--array={array}", f"--cpus-per-task={cpus}",
                        f"--time={hh:02d}:{mm:02d}:{ss:02d}", f"--comment={project}:{plan_id}",
                        f"--output={logs}_%A_%a.out", f"--error={logs}_%A_%a.err",
                        f"{remote_snapshot}/hpc/v2/run_array.sbatch", remote_snapshot, name, str(nchunks)]
                if policy == CORE_PAIRED_1H_POLICY:
                    argv.insert(8, CORE_PAIRED_1H_SIGNAL_OPTION)
                report["submission_plans"].append(identity | {"plan_id": plan_id, "argv": argv,
                    "shell_preview": shlex.join(argv), "requires_fresh_locked_dispatch": True,
                    "ledger_intent": identity | {"plan_id": plan_id, "state": "intent", "created_utc": now.isoformat()}})
                free -= cpus*parallel
                if free_slots is not None:
                    free_slots -= len(selected)
                for item in result["chunks"]:
                    if item["array_index"] in selected and item["state"] == "eligible_resume":
                        item["state"] = "planned_dry_run"
                result["state_counts"] = dict(Counter(c["state"] for c in result["chunks"]))
        except (ValueError, KeyError, TypeError, OSError) as exc:
            result["state"] = "blocked_validation"
            result["issues"].append(str(exc))
    report["cpu_budget"]["new_reserved"] = sum(p["cpus_per_task"]*p["max_concurrent"] for p in report["submission_plans"])
    report["cpu_budget"]["remaining_after_plan"] = free
    report["user_job_budget"].update(new_reserved=sum(len(p["array_task_ids"]) for p in report["submission_plans"]),
                                     remaining_after_plan=free_slots)
    report["all_registered_runs_complete"] = all(r["state"] == "complete" for r in report["runs"])
    report["notes"] = ["No sbatch/scontrol/scancel was executed, and no ledger was mutated.",
        "Project CPU reservations include pending arrays at their throttle; unrelated projects are reported only.",
        "The optional user submitted-job cap counts every active array member and ordinary job across all projects, plus unseen open-receipt members without strict terminal sacct evidence. It does not use CPU throttle counts.",
        "Only accepted receipt members with exact terminal sacct evidence and absent from active queue targets can be resumed before their array parent ends; unknown receipts are never released.",
        "Slot-limited batches retain the same nchunks, cell mapping and scientific seeds; open receipts still block duplicate targets.",
        "A plan is not a submission receipt. Dispatch under one project lock after fresh queue/ledger reconciliation.",
        "Terminal failed fits remain in planned denominators; this tool never deletes or retries failed permutation IDs.",
        "Inventory does not read large fit logs. Runner owns exact event replay and rejects duplicate/interior-corrupt logs.",
        "Passing a weekly usage percentage >=90 prevents new plans; this tool does not query usage or cancel running jobs."]
    return report


def markdown(plan):
    lines = ["# Frozen simulation campaign plan", "", f"Generated: {plan['created_utc']}", "",
             "Mode: **dry run only**. No jobs were submitted or changed.", "",
             f"Project CPU cap: {plan['cpu_budget']['cap']}; reserved: {plan['cpu_budget']['project_reserved']}; "
             f"proposed additional: {plan['cpu_budget']['new_reserved']}.", "",
             "| Run | State | Chunk states |", "|---|---|---|"]
    for run in plan["runs"]:
        states = dict(Counter(c["state"] for c in run["chunks"]))
        lines.append(f"| {run['run_id']} | {run['state']} | {json.dumps(states)} |")
    slots = plan.get("user_job_budget", {})
    if slots.get("cap") is not None:
        lines.extend(["", f"User submitted-job cap: {slots['cap']}; active across all projects: "
                      f"{slots['queue_active_members']}; unseen receipt members: {slots['unresolved_submission_reserved']}; "
                      f"proposed members: {slots['new_reserved']}; remaining: {slots['remaining_after_plan']}."])
    for run in plan["runs"]:
        for issue in run["issues"]:
            lines.extend(["", f"- {run['run_id']}: {json.dumps(issue) if isinstance(issue, dict) else issue}"])
    for submission in plan["submission_plans"]:
        lines.extend(["", f"## Proposed {submission['run_id']} ({submission['plan_id'][:12]})", "",
                      "The dispatcher must revalidate and write the intent before using this preview.", "",
                      "```sh", submission["shell_preview"], "```"])
    lines.extend(["", *[f"- {note}" for note in plan["notes"]], ""])
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--queue", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--weekly-used-percent", type=float)
    args = parser.parse_args(argv)
    try:
        plan = make_plan(read_json(args.campaign), read_json(args.queue), read_json(args.ledger),
                         weekly_used=args.weekly_used_percent)
    except (ValueError, TypeError, KeyError, OSError) as exc:
        plan = {"schema_version": "campaign_plan_v1", "created_utc": utc().isoformat(),
                "mode": "dry_run_only", "submission_plans": [], "error": str(exc)}
        write_json(args.output, plan)
        if args.markdown:
            args.markdown.parent.mkdir(parents=True, exist_ok=True)
            args.markdown.write_text(f"# Campaign planning blocked\n\n{plan['created_utc']}\n\n{exc}\n\nNo submission plan was generated.\n")
        return 2
    write_json(args.output, plan)
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(markdown(plan))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
