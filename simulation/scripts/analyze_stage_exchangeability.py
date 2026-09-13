"""Post-run DX x true-stage overlap diagnostics; no fitting or new experiment.

The CLI requires explicit existing chunk evidence and explicit seeds. It verifies
all selected observed records BEFORE importing the immutable design module or
reconstructing any data. Never point it at unrun confirmation seeds. Output is
JSON/Markdown only; it neither edits logs nor evaluates permutation test results.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import math
from numbers import Integral
import os
from pathlib import Path
import re
import sys
import tempfile

GROUPS = ("e2", "e33", "e4")
DIAGNOSES = ("CN", "MCI", "AD")
PAIRS = ((0, 1), (0, 2), (1, 2))
HASH = re.compile(r"[0-9a-f]{64}\Z")
RUNTIME_NAMES = frozenset({"run_v2.py", "design_v2.py", "engine_v2.py", "dependency_v2.py",
                           "oracle_engine_v2.py", "oracle_estimated_prior_v2.py",
                           "paired_engine_v2.py", "fast_likelihood_v2.py"})
TERMINAL_OBSERVED_STATUSES = frozenset({"ok", "error", "invalid"})
INTERPRETATION = (
    "These are realized overlap and distinct-label-assignment diagnostics, not an "
    "empirical proof of conditional exchangeability or population positivity. "
    "The oracle requires common measurement/missingness input laws conditional on "
    "Diagnosis and true discrete stage. Single-group blocks remain fixed: they "
    "can make an otherwise valid conditional randomization uninformative. Low "
    "rejection with a degenerate or small randomization space cannot establish "
    "nontrivial stage-oracle calibration. A large label space likewise need not "
    "produce many distinct ordering statistics. No rows are excluded or replaced "
    "because of poor overlap."
)


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _log_multinomial(counts):
    positive = [n for n in counts if n]
    if len(positive) <= 1:
        return 0.0
    return max(0.0, math.lgamma(sum(positive) + 1) -
               math.fsum(math.lgamma(n + 1) for n in positive))


def _block(dx, stage, counts):
    n = sum(counts)
    represented = sum(x > 0 for x in counts)
    movable = n if represented >= 2 else 0
    expected = n - sum(x * x for x in counts) / n if n else 0.0
    return {"diagnosis": dx, "stage": stage, "n": n,
            "group_counts": dict(zip(GROUPS, counts)),
            "represented_groups": represented,
            "single_group": represented == 1,
            "movable_sample_count": movable,
            "expected_label_changes_one_uniform_permutation": expected,
            "log_distinct_label_assignments": _log_multinomial(counts)}


def _summarize(blocks):
    n = sum(b["n"] for b in blocks)
    occupied = [b for b in blocks if b["n"]]
    single = [b for b in occupied if b["single_group"]]
    movable = sum(b["movable_sample_count"] for b in blocks)
    log_space = math.fsum(b["log_distinct_label_assignments"] for b in blocks)
    pairs = {}
    for a, b in PAIRS:
        ga, gb = GROUPS[a], GROUPS[b]
        shared = [x for x in blocks if x["group_counts"][ga] and x["group_counts"][gb]]
        pair_n = sum(x["group_counts"][ga] + x["group_counts"][gb] for x in blocks)
        shared_n = sum(x["group_counts"][ga] + x["group_counts"][gb] for x in shared)
        pairs[f"{ga}_vs_{gb}"] = {
            "blocks_containing_both_groups": len(shared),
            "pair_samples": pair_n,
            "pair_samples_in_shared_blocks": shared_n,
            "pair_sample_overlap_fraction": shared_n / pair_n if pair_n else None,
            "log_pair_only_label_assignments": math.fsum(
                _log_multinomial([x["group_counts"][ga], x["group_counts"][gb]]) for x in shared),
            "note": "Pair-only space is an overlap diagnostic; confirm_stage_a permutes all three labels."
        }
    by_group = {}
    for group in GROUPS:
        total = sum(b["group_counts"][group] for b in blocks)
        free = sum(b["group_counts"][group] for b in blocks if b["represented_groups"] >= 2)
        by_group[group] = {"n": total, "movable_samples": free,
                           "movable_fraction": free / total if total else None}
    return {
        "n": n, "possible_blocks": len(blocks), "occupied_blocks": len(occupied),
        "empty_blocks": len(blocks) - len(occupied),
        "single_group_blocks": len(single),
        "single_group_block_fraction_among_occupied": len(single) / len(occupied) if occupied else None,
        "samples_in_single_group_blocks": sum(b["n"] for b in single),
        "single_group_sample_fraction": sum(b["n"] for b in single) / n if n else None,
        "all_three_group_blocks": sum(b["represented_groups"] == 3 for b in blocks),
        "movable_samples": movable, "movable_sample_fraction": movable / n if n else None,
        "expected_label_changes_one_uniform_permutation": math.fsum(
            b["expected_label_changes_one_uniform_permutation"] for b in blocks),
        "log_distinct_label_assignments": log_space,
        "log10_distinct_label_assignments": log_space / math.log(10),
        "randomization_space_degenerate": log_space == 0,
        "by_group": by_group, "pair_overlap": pairs,
        "blocks": blocks
    }


def block_metrics(diagnoses, stages, labels, *, max_stage):
    """Summarize all DX x k blocks, including empty ones, and the DX-only comparator.

    ``movable`` means a row can change APOE under at least one legal label
    assignment, not that a random permutation always changes it. ln(Omega) is
    sum_b ln[n_b! / product_g n_bg!], the number of distinct label vectors, not
    n_b! index permutations and not the number of distinct fitted statistics.
    """
    diagnoses, stages, labels = list(diagnoses), list(stages), list(labels)
    max_stage = _integer(max_stage, "max_stage")
    if not diagnoses or not len(diagnoses) == len(stages) == len(labels):
        raise ValueError("Nonempty diagnosis/stage/label arrays must have equal lengths")
    counts = Counter()
    for dx, stage, group in zip(diagnoses, stages, labels):
        stage, group = _integer(stage, "stage"), _integer(group, "APOE")
        if dx not in DIAGNOSES or stage > max_stage or group >= len(GROUPS):
            raise ValueError("Unknown diagnosis, stage outside 0..I, or APOE outside0..2")
        counts[dx, stage, group] += 1
    oracle = [_block(dx, k, [counts[dx, k, g] for g in range(3)])
              for dx in DIAGNOSES for k in range(max_stage + 1)]
    dx_only = [_block(dx, None, [sum(counts[dx, k, g] for k in range(max_stage + 1))
                                for g in range(3)]) for dx in DIAGNOSES]
    return {"group_order": list(GROUPS), "max_stage": max_stage,
            "oracle_dx_stage": _summarize(oracle), "diagnosis": _summarize(dx_only),
            "interpretation": INTERPRETATION}


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _read_json(path):
    return json.loads(Path(path).read_text())


def _verify_snapshot(workspace, run):
    snapshot = (workspace / run["snapshot"]).resolve()
    if not snapshot.is_relative_to(workspace):
        raise ValueError("Snapshot path escapes workspace")
    checks = snapshot / "SHA256SUMS"
    if _sha(checks) != run["snapshot_sha256"]:
        raise ValueError("Snapshot manifest hash differs from frozen protocol")
    declared = {}
    for line in checks.read_text().splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if not match:
            raise ValueError("Malformed snapshot SHA256SUMS")
        digest, name = match.groups()
        path = (snapshot / name).resolve()
        if name in declared or not path.is_relative_to(snapshot) or path == snapshot:
            raise ValueError("Unsafe or duplicate snapshot path")
        if _sha(path) != digest:
            raise ValueError(f"Frozen source bytes changed: {name}")
        declared[name] = digest
    required = {run["config"], "scripts/v2/design_v2.py", "scripts/v2/run_v2.py"}
    if not required <= declared.keys() or declared[run["config"]] != run["config_sha256"]:
        raise ValueError("Frozen config or generator/runner is not correctly checksummed")
    return snapshot, declared, _read_json(snapshot / run["config"])


def _observed_evidence(path, seeds):
    """Read only the initial file prefix; ignore an uncommitted tail without editing it."""
    found, digest = {}, hashlib.sha256()
    size = Path(path).stat().st_size
    tail = None
    with Path(path).open("rb") as stream:
        remaining = size
        while remaining:
            line = stream.readline(remaining)
            if not line:
                raise ValueError("Fit log was truncated during evidence read")
            remaining -= len(line)
            digest.update(line)
            if not line.endswith(b"\n"):
                tail = hashlib.sha256(line).hexdigest()
                continue
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except (ValueError, UnicodeDecodeError) as exc:
                raise ValueError("Corrupt committed fit record; no reconstruction attempted") from exc
            if not isinstance(record, dict):
                raise ValueError("Fit record must be an object")
            if record.get("kind") != "observed" or record.get("seed") not in seeds:
                continue
            if record.get("engine") != "repaired":
                continue
            seed = _integer(record["seed"], "observed seed")
            if record.get("scheme") is not None or record.get("perm_id") is not None:
                raise ValueError("Observed evidence has permutation fields")
            if record.get("status") not in TERMINAL_OBSERVED_STATUSES:
                raise ValueError("Observed evidence must have a terminal fit status")
            if seed in found:
                raise ValueError("Duplicate selected observed record")
            found[seed] = record
    missing = set(seeds) - found.keys()
    if missing:
        raise ValueError(f"No committed repaired observed evidence for seeds {sorted(missing)}; generation forbidden")
    return found, {"fit_records": str(Path(path).resolve()), "prefix_bytes": size,
                   "prefix_sha256": digest.hexdigest(), "ignored_uncommitted_tail_sha256": tail}


def _load_frozen_design(path, expected_versions, expected_sha256=None):
    # This is intentionally after evidence validation. No engine or runner import.
    import numpy as np
    import pandas as pd
    actual = {"numpy_version": np.__version__, "pandas_version": pd.__version__}
    if actual != expected_versions:
        raise ValueError(f"Reconstruct in the original NumPy/Pandas environment: {expected_versions}; local {actual}")
    # SourceFileLoader may READ an unchecked pre-existing .pyc even when
    # dont_write_bytecode=True. Compile the verified bytes directly instead.
    source = Path(path).read_bytes()
    digest = hashlib.sha256(source).hexdigest()
    if expected_sha256 is not None and digest != expected_sha256:
        raise ValueError("Frozen generator changed before source compilation")
    name = "_stage_diagnostic_frozen_design_" + digest
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    previous_module = sys.modules.get(name)
    sys.modules[name] = module  # dataclasses requires the defining module while loading.
    old = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        exec(compile(source, str(path), "exec"), module.__dict__)
    finally:
        sys.dont_write_bytecode = old
        if previous_module is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = previous_module
    return module


def analyze_chunk(protocol_path, chunk_dir, seeds, *, workspace=None):
    """Verify existing evidence, then reconstruct ONLY explicitly requested observed seeds.

    Missing observed evidence rejects the entire selection before any generation.
    An observed failure with no saved generation hash is retained as unavailable,
    without reconstruction. A saved-generation/hash mismatch is a hard error.
    """
    workspace = Path(workspace or Path(__file__).resolve().parents[2]).resolve()
    protocol_path, chunk_dir = Path(protocol_path), Path(chunk_dir)
    seeds = [_integer(x, "requested seed") for x in seeds]
    if not seeds or len(seeds) != len(set(seeds)):
        raise ValueError("Select a nonempty, unique explicit seed list")
    protocol = _read_json(protocol_path)
    runs = [r for r in protocol["runs"] if r["run_id"] == "confirm_stage_a"]
    if len(runs) != 1:
        raise ValueError("Exactly one frozen confirm_stage_a protocol entry is required")
    run = runs[0]
    snapshot, declared, config = _verify_snapshot(workspace, run)
    definition = config.get("schemes", {}).get("oracle_dx_stage", {})
    if (config["run_id"] != "confirm_stage_a" or config["engines"] != ["repaired"]
            or definition.get("stratify") != "oracle_dx_stage"
            or definition.get("groups", [0, 1, 2]) != [0, 1, 2]):
        raise ValueError("Unsupported frozen stage scheme/engine; do not guess strata")
    manifest = _read_json(chunk_dir / "manifest.json")
    if manifest["config"] != config:
        raise ValueError("Chunk config differs from frozen stage config")
    cell_index = _integer(manifest["cell_index"], "cell index")
    chunk = _integer(manifest["chunk"], "chunk")
    nchunks = _integer(manifest["nchunks"], "nchunks", 1)
    if cell_index >= len(config["cells"]) or chunk >= nchunks:
        raise ValueError("Invalid chunk/cell layout")
    planned = list(range(config["base_seed"], config["base_seed"] + config["datasets"]))[chunk::nchunks]
    if manifest["seeds"] != planned or not set(seeds) <= set(planned):
        raise ValueError("Selected seeds are not in the manifest's frozen chunk layout")
    sources = manifest["source_sha256"]
    expected_sources = {name: declared["scripts/v2/" + name] for name in RUNTIME_NAMES
                        if "scripts/v2/" + name in declared}
    if sources != expected_sources:
        raise ValueError("Chunk runtime source closure differs from immutable snapshot")
    identity = {"config": config, "cell_index": cell_index, "chunk": chunk,
                "nchunks": nchunks, "source_sha256": sources}
    if hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest() != manifest["signature"]:
        raise ValueError("Chunk manifest signature mismatch")
    observed, evidence = _observed_evidence(chunk_dir / "fit_records.jsonl", set(seeds))
    # All selections now have committed observed records. Validate all available
    # generation evidence before the first call to simulate.
    available = {}
    for seed in seeds:
        record, truth = observed[seed], observed[seed].get("truth", {})
        if not isinstance(truth, dict) or not isinstance(truth.get("manifest", {}), dict):
            raise ValueError("Saved generation truth/manifest must be objects")
        tm = truth.get("manifest", {})
        if tm.get("data_sha256") is None:
            if record["status"] == "ok":
                raise ValueError("Successful observed fit lacks saved generation data hash")
            continue
        if not isinstance(tm["data_sha256"], str) or not HASH.fullmatch(tm["data_sha256"]):
            raise ValueError("Invalid saved data hash")
        if (truth.get("seed") != seed or tm.get("seed") != seed or
                tm.get("source_sha256") != declared["scripts/v2/design_v2.py"]):
            raise ValueError("Saved truth seed/generator provenance mismatch")
        if truth.get("group_order") != list(GROUPS):
            raise ValueError("Unexpected saved truth group ordering")
        for key in ("config", "biomarker_names", "latent_stage", "latent_row_ptid"):
            if key not in truth:
                raise ValueError(f"Saved truth lacks {key}")
        config_hash = hashlib.sha256(json.dumps(truth["config"], sort_keys=True,
            separators=(",", ":"), allow_nan=False).encode()).hexdigest()
        if tm.get("config_sha256") != config_hash:
            raise ValueError("Saved generation config hash mismatch")
        environment = manifest.get("environment", {}).get("versions", {})
        if any(not isinstance(environment.get(package), str) or
               tm.get(package + "_version") != environment[package] for package in ("numpy", "pandas")):
            raise ValueError("Saved generation versions differ from chunk environment")
        names, stages, ptids = truth["biomarker_names"], truth["latent_stage"], truth["latent_row_ptid"]
        if (not isinstance(names, list) or not names or any(not isinstance(x, str) for x in names)
                or len(set(names)) != len(names) or not isinstance(stages, list)
                or not isinstance(ptids, list) or not stages or len(stages) != len(ptids)
                or any(not isinstance(x, str) for x in ptids) or len(set(ptids)) != len(ptids)):
            raise ValueError("Invalid saved biomarker/latent-stage/PTID layout")
        if any(_integer(k, "saved latent stage") > len(names) for k in stages):
            raise ValueError("Saved latent stage outside 0..I")
        available[seed] = truth
    rows = []
    cell = config["cells"][cell_index]
    prepared = {}
    if available:
        first = next(iter(available.values()))
        versions = {k: first["manifest"][k] for k in ("numpy_version", "pandas_version")}
        design = _load_frozen_design(snapshot / "scripts/v2/design_v2.py", versions,
                                     declared["scripts/v2/design_v2.py"])
        cfg = design.resolve(cell["name"], **cell.get("overrides", {}))
        for seed, saved in available.items():
            if cfg.to_dict() != saved["config"]:
                raise ValueError("Resolved frozen design does not match saved truth config")
            prepared[seed] = (design, cfg)
    for seed in seeds:
        record = observed[seed]
        if seed not in available:
            rows.append({"seed": seed, "observed_status": record.get("status"),
                         "diagnostic_status": "unavailable_saved_generation_hash_missing",
                         "reconstructed": False})
            continue
        saved = available[seed]
        design, cfg = prepared[seed]
        df, truth = design.simulate(cfg, seed)
        data_hash = hashlib.sha256(df.to_csv(index=False, float_format="%.17g").encode()).hexdigest()
        if data_hash != saved["manifest"]["data_sha256"] or truth["manifest"]["data_sha256"] != data_hash:
            raise ValueError("Reconstructed data SHA differs from original observed data; no diagnostic accepted")
        if (truth["latent_stage"] != saved["latent_stage"] or
                df.PTID.tolist() != saved["latent_row_ptid"] or
                truth.get("latent_row_ptid") != saved["latent_row_ptid"] or
                truth["manifest"].get("source_sha256") != declared["scripts/v2/design_v2.py"] or
                truth["manifest"].get("config_sha256") != saved["manifest"]["config_sha256"] or
                truth["biomarker_names"] != saved["biomarker_names"]):
            raise ValueError("Latent-stage/row alignment differs from saved truth")
        rows.append({"seed": seed, "observed_status": record.get("status"),
                     "diagnostic_status": "verified", "reconstructed": True,
                     "data_sha256": data_hash,
                     "metrics": block_metrics(df.Diagnosis, truth["latent_stage"], df.APOE,
                                              max_stage=len(truth["biomarker_names"]))})
    return {"schema_version": "stage-exchangeability-diagnostics-v1", "run_id": config["run_id"],
            "cell": cell.get("id", cell["name"]), "requested_seeds": seeds,
            "planned_chunk_datasets": len(planned), "observed_evidence_required": True,
            "verified_datasets": sum(r["diagnostic_status"] == "verified" for r in rows),
            "observed_status_counts": dict(Counter(r["observed_status"] for r in rows)),
            "diagnostic_status_counts": dict(Counter(r["diagnostic_status"] for r in rows)),
            "selection_note": "Explicit existing-record subset only; not all planned seeds or an unbiased completion sample.",
            "protocol_path": str(protocol_path.resolve()), "protocol_sha256": _sha(protocol_path),
            "snapshot": str(snapshot), "snapshot_sha256": run["snapshot_sha256"],
            "chunk_manifest_sha256": _sha(chunk_dir / "manifest.json"),
            "generator_sha256": declared["scripts/v2/design_v2.py"],
            "oracle_scheme_definition": definition,
            "fitting_performed": False, "new_experiment": False, "evidence": evidence,
            "interpretation": INTERPRETATION, "datasets": rows}


def render_markdown(result):
    lines = ["# Stage-oracle block diagnostics", "", INTERPRETATION, "",
             f"Run `{result['run_id']}`, cell `{result['cell']}`: "
             f"{result['verified_datasets']}/{len(result['requested_seeds'])} explicitly requested datasets verified.",
             "This is an existing-record subset; execution failures remain in the main analysis denominator.", "",
             "| Seed | Original fit | Diagnostic | Occupied blocks | Single-group blocks | Movable samples | ln(label space) |",
             "|---|---|---|---:|---:|---:|---:|"]
    for row in result["datasets"]:
        if row["diagnostic_status"] != "verified":
            lines.append(f"| {row['seed']} | {row['observed_status']} | unavailable hash | — | — | — | — |")
        else:
            m = row["metrics"]["oracle_dx_stage"]
            lines.append(f"| {row['seed']} | {row['observed_status']} | verified | {m['occupied_blocks']} | "
                         f"{m['single_group_blocks']} | {m['movable_samples']}/{m['n']} | "
                         f"{m['log_distinct_label_assignments']:.6f} |")
    lines += ["", "JSON includes every DX × true discrete-stage block (empty blocks included), "
              "group counts, pair overlap, DX-only comparisons, exact denominators and provenance.", ""]
    return "\n".join(lines)


def validate_output_paths(protocol_path, chunk_dir, outputs, *, workspace=None):
    """Reject output aliases to frozen inputs or original chunk evidence first."""
    workspace = Path(workspace or Path(__file__).resolve().parents[2]).resolve()
    protocol_path = Path(protocol_path).resolve()
    protocol = _read_json(protocol_path)
    protected_files = {protocol_path}
    protected_roots = {Path(chunk_dir).resolve()}
    for run in protocol["runs"]:
        protected_files.add((workspace / run["config"]).resolve())
        protected_roots.add((workspace / run["snapshot"]).resolve())
    targets = [Path(path).resolve() for path in outputs if path is not None]
    if len(targets) != len(set(targets)):
        raise ValueError("JSON and Markdown outputs must be different paths")
    for target in targets:
        if target in protected_files or any(target.is_relative_to(root) for root in protected_roots):
            raise ValueError("Diagnostic output would overwrite frozen inputs or original chunk evidence")
        if target.exists() and target.stat().st_nlink > 1:
            raise ValueError("Diagnostic output has multiple hard links; refusing possible input overwrite")


def write_outputs(protocol_path, chunk_dir, outputs, *, workspace=None):
    """Recheck destinations and atomically replace reports without writing through aliases."""
    outputs = [(Path(path), content) for path, content in outputs if path is not None]
    destinations = [path for path, _ in outputs]
    validate_output_paths(protocol_path, chunk_dir, destinations, workspace=workspace)
    for path, content in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                             prefix=f".{path.name}.", suffix=".tmp", delete=False) as stream:
                temporary = Path(stream.name)
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
            validate_output_paths(protocol_path, chunk_dir, destinations, workspace=workspace)
            # Even if an alias appeared after the last check, replacement does
            # not truncate the inode shared with an existing input hard link.
            os.replace(temporary, path)
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--chunk-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, action="append", required=True,
                        help="Explicit seed already confirmed to have a committed HPC observed record; repeatable")
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path)
    args = parser.parse_args()
    validate_output_paths(args.protocol, args.chunk_dir, [args.output_json, args.output_md], workspace=args.workspace)
    result = analyze_chunk(args.protocol, args.chunk_dir, args.seed, workspace=args.workspace)
    write_outputs(args.protocol, args.chunk_dir,
                  [(args.output_json, json.dumps(result, indent=2, allow_nan=False) + "\n"),
                   (args.output_md, render_markdown(result))], workspace=args.workspace)


if __name__ == "__main__":
    main()
