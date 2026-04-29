#!/usr/bin/env python3
"""Run hyperelliptic_finder_c shards in parallel and merge their JSON output."""

from __future__ import annotations

import argparse
import json
import select
import subprocess
import sys
import time
from pathlib import Path


DEFAULT_WORKERS = 6
SEPARATOR = "=" * 72


def output_paths(path: Path) -> tuple[Path, Path]:
    if path.suffix == ".json":
        return path.with_suffix(".txt"), path
    if path.suffix == ".txt":
        return path, path.with_suffix(".json")
    return path, path.with_suffix(path.suffix + ".json")


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(rows, handle, indent=2)
        handle.write("\n")


def ensure_binary(root: Path, binary: Path, arch: str, no_build: bool) -> None:
    if binary.exists():
        return
    if no_build:
        raise FileNotFoundError(f"{binary} does not exist; run make first")
    command = ["make"]
    if arch:
        command.append(f"ARCH={arch}")
    subprocess.run(command, cwd=root, check=True)


def saved_presentation_count(path: Path) -> int:
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return 0
    try:
        return int(data.get("total_presentations_found", 0))
    except (TypeError, ValueError):
        return 0


def format_polynomial(coeffs: list[int]) -> str:
    terms = []
    for exponent, coeff in enumerate(coeffs):
        if coeff == 0:
            continue
        if exponent == 0:
            terms.append(str(coeff))
        elif exponent == 1:
            terms.append("x" if coeff == 1 else f"{coeff}x")
        else:
            terms.append(f"x^{exponent}" if coeff == 1 else f"{coeff}x^{exponent}")
    return " + ".join(terms) if terms else "0"


def key_tuple(values: object) -> tuple[int, ...]:
    if not isinstance(values, list):
        return ()
    return tuple(int(v) for v in values)


def merge_worker_jsons(
    worker_jsons: list[Path],
    output_txt: Path,
    output_json: Path,
    dedupe: str,
) -> dict[str, object]:
    merged_curves: list[dict[str, object]] = []
    seen: set[tuple[int, ...]] = set()
    first_data: dict[str, object] | None = None
    worker_statuses = []

    for path in worker_jsons:
        if not path.exists():
            worker_statuses.append({"path": str(path), "status": "missing"})
            continue
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        if first_data is None:
            first_data = data
        worker_statuses.append(
            {
                "path": str(path),
                "status": data.get("search_status"),
                "complete_list": data.get("complete_list"),
                "total_presentations_found": data.get("total_presentations_found"),
            }
        )
        for curve in data.get("curves", []):
            if dedupe == "reduction-key":
                key = key_tuple(curve.get("reduction_key"))
            else:
                key = key_tuple(curve.get("f_coeffs"))
            if key in seen:
                continue
            seen.add(key)
            merged_curves.append(dict(curve))

    for index, curve in enumerate(merged_curves, start=1):
        curve["index"] = index

    class_keys = {key_tuple(curve.get("reduction_key")) for curve in merged_curves if curve.get("reduction_key") is not None}
    if not class_keys:
        class_keys = {key_tuple(curve.get("canonical_f_coeffs")) for curve in merged_curves}

    p = int(first_data.get("p", 0)) if first_data else 0
    g = int(first_data.get("g", 0)) if first_data else 0
    reduction = str(first_data.get("presentation_reduction", "")) if first_data else ""
    complete = bool(worker_statuses) and all(row.get("complete_list") is True for row in worker_statuses)
    status = "complete" if complete else "partial"

    output = {
        "p": p,
        "g": g,
        "field": f"F_{p}" if p else "",
        "presentation_reduction": reduction,
        "runner": "parallel_c",
        "dedupe": dedupe,
        "worker_count": len(worker_jsons),
        "worker_statuses": worker_statuses,
        "search_status": status,
        "complete_list": complete,
        "total_presentations_found": len(merged_curves),
        "reduction_class_count": len(class_keys),
        "isomorphism_class_count": len(class_keys) if reduction in ("pgl2", "pgl2save") else None,
        "curves": merged_curves,
    }

    with open(output_json, "w", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)
        handle.write("\n")

    with open(output_txt, "w", encoding="utf-8") as handle:
        handle.write("Hyperelliptic curves with trinomial L-polynomial\n")
        handle.write("=" * 56 + "\n")
        handle.write(f"p = {p}\n")
        handle.write(f"g = {g}\n")
        handle.write(f"base field: F_{p}\n")
        handle.write(f"presentation reduction: {reduction}\n")
        handle.write(f"runner: parallel_c\n")
        handle.write(f"dedupe: {dedupe}\n")
        handle.write(f"search status: {status}\n")
        handle.write(f"complete list: {complete}\n")
        handle.write(f"total presentations found: {len(merged_curves)}\n")
        handle.write(f"reduction class count: {len(class_keys)}\n\n")
        for curve in merged_curves:
            coeffs = [int(v) for v in curve.get("f_coeffs", [])]
            handle.write(f"[{curve['index']}]\n")
            handle.write(f"f(x) = {format_polynomial(coeffs)}\n")
            handle.write(f"middle_coefficient_a_{g} = {curve.get('middle_coefficient')}\n")
            handle.write("\n")

    return output


def main() -> int:
    parser = argparse.ArgumentParser(description="Run hyperelliptic_finder_c in parallel shards and merge output.")
    parser.add_argument("p", type=int, help="odd prime p")
    parser.add_argument("g", type=int, help="genus g")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS, help="number of worker processes; default: 6")
    parser.add_argument("--binary", default="hyperelliptic_finder_c", help="C finder binary path")
    parser.add_argument("--arch", default="-arch x86_64", help="ARCH value passed to make if the binary is missing")
    parser.add_argument("--no-build", action="store_true", help="do not run make if the C binary is missing")
    parser.add_argument("--outdir", default="parallel_results", help="directory for worker and merged output")
    parser.add_argument("--output", default="", help="merged output path; defaults inside --outdir")
    parser.add_argument("--reduction", choices=("pgl2", "affine", "pgl2save", "affinesave"), default="pgl2save")
    parser.add_argument("--max", type=int, default=0, help="max curves per worker; 0 means no per-worker max")
    parser.add_argument("--case-timeout", "--timeout", dest="timeout", type=int, default=1800)
    parser.add_argument("--dedupe", choices=("exact", "reduction-key"), default="exact")
    parser.add_argument("--quiet", action="store_true", help="pass --quiet to workers")
    parser.add_argument("--verbose", action="store_true", help="pass --verbose to workers")
    parser.add_argument("--monic-only", action="store_true")
    parser.add_argument("--no-hasse-witt-prefilter", action="store_true")
    parser.add_argument("--resume", action="store_true", help="reuse existing worker JSON files")
    args = parser.parse_args()

    if args.workers < 1:
        print("Error: --workers must be at least 1", file=sys.stderr)
        return 1

    root = Path(__file__).resolve().parent
    binary = Path(args.binary)
    if not binary.is_absolute():
        binary = root / binary
    ensure_binary(root, binary, args.arch, args.no_build)

    outdir = (root / args.outdir).resolve()
    worker_dir = outdir / f"p{args.p}_g{args.g}_{args.reduction}_workers"
    worker_dir.mkdir(parents=True, exist_ok=True)
    output_base = Path(args.output) if args.output else outdir / f"p{args.p}_g{args.g}_{args.reduction}_parallel.txt"
    if not output_base.is_absolute():
        output_base = root / output_base
    output_txt, output_json = output_paths(output_base)
    output_txt.parent.mkdir(parents=True, exist_ok=True)
    summary_path = worker_dir / "parallel_summary.json"

    processes: list[dict[str, object]] = []
    for worker_index in range(args.workers):
        worker_base = worker_dir / f"worker_{worker_index}.txt"
        worker_json = worker_base.with_suffix(".json")
        if args.resume and worker_json.exists():
            processes.append({"worker_index": worker_index, "status": "skipped", "json": worker_json})
            continue

        command = [
            str(binary),
            str(args.p),
            str(args.g),
            "--workers",
            str(args.workers),
            "--worker-index",
            str(worker_index),
            "--reduction",
            args.reduction,
            "--output",
            str(worker_base),
        ]
        if args.max > 0:
            command.extend(["--max", str(args.max)])
        if args.quiet:
            command.append("--quiet")
        if args.verbose:
            command.append("--verbose")
        if args.monic_only:
            command.append("--monic-only")
        if args.no_hasse_witt_prefilter:
            command.append("--no-hasse-witt-prefilter")

        process = subprocess.Popen(
            command,
            cwd=root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
        )
        processes.append(
            {
                "worker_index": worker_index,
                "process": process,
                "command": command,
                "json": worker_json,
                "stdout": [],
                "start": time.monotonic(),
                "timed_out": False,
            }
        )

    print(SEPARATOR, flush=True)
    print(f"Running p={args.p}, g={args.g} with {args.workers} C workers", flush=True)
    print(f"worker output: {worker_dir}", flush=True)
    print(f"merged output: {output_json}", flush=True)

    timeout = args.timeout if args.timeout > 0 else None
    unfinished = {i for i, row in enumerate(processes) if "process" in row}
    try:
        while unfinished:
            for row_index in list(unfinished):
                row = processes[row_index]
                process: subprocess.Popen[str] = row["process"]  # type: ignore[assignment]
                stdout = process.stdout
                if stdout is not None:
                    ready, _, _ = select.select([stdout], [], [], 0)
                    if ready:
                        line = stdout.readline()
                        if line:
                            row["stdout"].append(line)  # type: ignore[index]
                            print(f"[w{row['worker_index']}] {line}", end="" if line.endswith("\n") else "\n", flush=True)

                if (
                    timeout is not None
                    and not row["timed_out"]
                    and time.monotonic() - float(row["start"]) > timeout
                    and saved_presentation_count(row["json"]) > 0  # type: ignore[arg-type]
                ):
                    row["timed_out"] = True
                    process.terminate()

                if process.poll() is not None:
                    if stdout is not None:
                        for line in stdout:
                            row["stdout"].append(line)  # type: ignore[index]
                            print(f"[w{row['worker_index']}] {line}", end="" if line.endswith("\n") else "\n", flush=True)
                    row["returncode"] = process.returncode
                    row["status"] = "timeout" if row["timed_out"] else "ok" if process.returncode == 0 else "failed"
                    unfinished.remove(row_index)
            time.sleep(0.1)
    except KeyboardInterrupt:
        for row_index in unfinished:
            process = processes[row_index]["process"]
            process.terminate()  # type: ignore[union-attr]
        raise

    worker_rows = []
    worker_jsons = []
    for row in processes:
        worker_jsons.append(row["json"])  # type: ignore[arg-type]
        clean = {k: v for k, v in row.items() if k not in {"process"}}
        if isinstance(clean.get("json"), Path):
            clean["json"] = str(clean["json"])
        worker_rows.append(clean)
    write_summary(summary_path, worker_rows)

    merged = merge_worker_jsons(worker_jsons, output_txt, output_json, args.dedupe)
    print(SEPARATOR, flush=True)
    print(f"Merged {merged['total_presentations_found']} presentations", flush=True)
    print(f"Reduction classes: {merged['reduction_class_count']}", flush=True)
    print(f"Wrote {output_txt} and {output_json}", flush=True)
    print(f"Wrote {summary_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
