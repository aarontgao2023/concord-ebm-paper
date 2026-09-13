"""Standard-library helpers for read-only, planned-identity development audits."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "hpc/v2"))
import inventory_run as I
import analyze_v2 as A

GROUPS = ("e2", "e33", "e4")


def frozen_config(run):
    return ROOT / "runs/v2/snapshots" / run / "configs/v2" / (run + ".json")


def cell_id(spec):
    return spec.get("id", spec["name"])


def identity(cell, record):
    return cell, record["seed"], record["engine"], record.get("scheme"), record.get("perm_id")


def load_run(config_path, root, wanted=None, check_inventory=True):
    """Stream fit logs, reject duplicate identities, and retain input digests.

    wanted restricts a *reference* run to requested identities; irrelevant chunks
    are skipped by signed manifest seed allocation. It does not waive missingness.
    A torn line is reported and never repaired by this read-only analysis.
    """
    config_path, root = Path(config_path), Path(root)
    cfg = json.loads(config_path.read_text())
    records, inputs, manifests, issues, duplicates = {}, [], {}, [], set()
    seeds = set(range(cfg["base_seed"], cfg["base_seed"] + cfg["datasets"]))
    if wanted is not None:
        seeds &= {k[1] for k in wanted}
    for path in sorted(root.glob("cell_*/chunk_*/manifest.json")):
        try:
            manifest = json.loads(path.read_text())
            index = int(path.parent.parent.name.removeprefix("cell_"))
            cell = cell_id(cfg["cells"][index])
            if not seeds.intersection(manifest.get("seeds", [])):
                continue
            signed = {name: manifest[name] for name in
                      ("config", "cell_index", "chunk", "nchunks", "source_sha256")}
            if manifest["signature"] != I.digest(signed) or manifest["config"] != cfg:
                raise ValueError("Manifest signature or frozen config mismatch")
            manifests[str(path)] = manifest
            inputs.append({"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
            log = path.parent / "fit_records.jsonl"
            digest, n, selected = hashlib.sha256(), 0, 0
            before = log.stat()
            with log.open("rb") as stream:
                for n, line in enumerate(stream, 1):
                    digest.update(line)
                    if not line.strip():
                        continue
                    r = json.loads(line)
                    k = identity(cell, r)
                    if wanted is not None and k not in wanted:
                        continue
                    selected += 1
                    if k in records or k in duplicates:
                        duplicates.add(k)
                        records.pop(k, None)
                        issues.append({"code": "duplicate_fit_identity", "key": list(k), "path": str(log), "line": n})
                        continue
                    records[k] = r
            after = log.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                issues.append({"code": "input_changed_during_read", "path": str(log)})
            inputs.append({"path": str(log.resolve()), "sha256": digest.hexdigest(),
                           "n_lines": n, "n_selected": selected})
        except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
            issues.append({"code": "invalid_or_missing_input", "path": str(path),
                           "message": f"{type(exc).__name__}: {exc}"})
    inventory = I.inventory(config_path, root) if check_inventory else None
    return {"config": cfg, "records": records, "inputs": inputs, "issues": issues,
            "manifests": manifests, "inventory": inventory,
            "config_path": str(config_path.resolve()),
            "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest()}


def expected_observed(cfg, engines=None):
    return {(cell_id(cell), seed, engine, None, None)
            for cell in cfg["cells"]
            for seed in range(cfg["base_seed"], cfg["base_seed"] + cfg["datasets"])
            for engine in (cfg["engines"] if engines is None else engines)}


def data_hash(record):
    return record.get("truth", {}).get("manifest", {}).get("data_sha256")


def compare_fit(a, b):
    if a is None or b is None:
        return {"both_present": False, "successful_exact": False,
                "missing": "both" if a is b else "new" if a is None else "reference"}
    checks = {name: a.get(name) == b.get(name) for name in ("status", "orderings", "taus")}
    success = a.get("status") == b.get("status") == "ok"
    return {"both_present": True, "both_successful": success, "checks": checks,
            "successful_exact": success and all(checks.values()),
            "new_status": a.get("status"), "reference_status": b.get("status")}


def kendall(a, b):
    if len(a) != len(b) or sorted(a) != sorted(b) or len(set(a)) != len(a) or len(a) < 2:
        raise ValueError("Ordering must be a unique permutation with at least two events")
    rank = {v: i for i, v in enumerate(b)}
    return sum(rank[a[i]] > rank[a[j]] for i in range(len(a)) for j in range(i + 1, len(a))) / (len(a) * (len(a) - 1) / 2)


def numbers(values):
    values = list(values)
    return {"n": len(values), "mean": statistics.mean(values) if values else None,
            "median": statistics.median(values) if values else None,
            "min": min(values, default=None), "max": max(values, default=None)}


def paired_effect(a, b, context, repeats=10000):
    differences = [y - x for x, y in zip(a, b)]
    return {"first": numbers(a), "second": numbers(b),
            "difference_second_minus_first": A.continuous_summary(differences, repeats, 20260908, context),
            "n_lower": sum(v < -1e-12 for v in differences),
            "n_equal": sum(abs(v) <= 1e-12 for v in differences),
            "n_higher": sum(v > 1e-12 for v in differences)}


def manifest_profiles(bundle):
    fields = ("versions", "source_sha256", "wheel_sha256")
    return {I.digest({k: m.get("environment", {}).get(k) for k in fields})
            for m in bundle["manifests"].values()}


def fmt(value):
    return "NA" if value is None else f"{value:.6g}"


def save(out, stem, payload, lines):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / (stem + ".json")).write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    (out / (stem + ".md")).write_text("\n".join(lines) + "\n")
    print(json.dumps({"output": str(out.resolve()), "ready": payload.get("ready", False),
                      "n_issues": len(payload.get("issues", []))}))
