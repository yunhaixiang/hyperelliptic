#!/usr/bin/env python3
"""Batch-run hyperelliptic_finder_c for p=3 and genera 7 through 20."""

from __future__ import annotations

import argparse
import json
import select
import subprocess
import sys
import time
from pathlib import Path


PRIMES = (3,)
DEFAULT_MIN_GENUS = 7
DEFAULT_MAX_GENUS = 20
SEPARATOR = "=" * 72


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data: str) -> int:
        for stream in self.streams:
            stream.write(data)
            stream.flush()
        return len(data)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()


def write_summary(path: Path, summary: list[dict[str, object]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")


def print_case_header(p: int, g: int, reduction: str, output_json: Path) -> None:
    print("\n" + SEPARATOR, flush=True)
    print(f"CASE p={p}, g={g}, reduction={reduction}, runner=C", flush=True)
    print(f"output: {output_json}", flush=True)


def print_case_footer(case: dict[str, object]) -> None:
    stderr = str(case.get("stderr") or "")
    if stderr:
        print("stderr:", flush=True)
        print(stderr, end="" if stderr.endswith("\n") else "\n", flush=True)


def saved_presentation_count(output_json: Path) -> int:
    try:
        with open(output_json, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return 0
    try:
        return int(data.get("total_presentations_found", 0))
    except (TypeError, ValueError):
        return 0


def ensure_binary(root: Path, binary: Path, arch: str, no_build: bool) -> None:
    if binary.exists():
        return
    if no_build:
        raise FileNotFoundError(f"{binary} does not exist; run make first")

    command = ["make"]
    if arch:
        command.append(f"ARCH={arch}")
    subprocess.run(command, cwd=root, check=True)


def run_case(
    binary: Path,
    outdir: Path,
    p: int,
    g: int,
    reduction: str,
    max_curves: int,
    timeout: int | None,
    quiet: bool,
    verbose: bool,
    monic_only: bool,
    no_hasse_witt_prefilter: bool,
) -> dict[str, object]:
    mode_suffix = "_monic" if monic_only else ""
    output_base = outdir / f"p{p}_g{g}_{reduction}_c{mode_suffix}.txt"
    command = [
        str(binary),
        str(p),
        str(g),
        "--reduction",
        reduction,
        "--output",
        str(output_base),
    ]
    if max_curves > 0:
        command.extend(["--max", str(max_curves)])
    if quiet:
        command.append("--quiet")
    if verbose:
        command.append("--verbose")
    if monic_only:
        command.append("--monic-only")
    if no_hasse_witt_prefilter:
        command.append("--no-hasse-witt-prefilter")

    case = {
        "p": p,
        "g": g,
        "runner": "c",
        "reduction": reduction,
        "output_txt": str(output_base),
        "output_json": str(output_base.with_suffix(".json")),
        "command": command,
    }
    output_lines: list[str] = []
    start = time.monotonic()
    process = subprocess.Popen(
        command,
        cwd=binary.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
    )
    assert process.stdout is not None

    timed_out = False
    while True:
        if timeout is not None and time.monotonic() - start > timeout and saved_presentation_count(output_base.with_suffix(".json")) > 0:
            timed_out = True
            process.terminate()
            break

        ready, _, _ = select.select([process.stdout], [], [], 0.1)
        if ready:
            line = process.stdout.readline()
            if line:
                output_lines.append(line)
                print(line, end="" if line.endswith("\n") else "\n", flush=True)

        if process.poll() is not None:
            for line in process.stdout:
                output_lines.append(line)
                print(line, end="" if line.endswith("\n") else "\n", flush=True)
            break

    if timed_out:
        try:
            remaining, _ = process.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            remaining, _ = process.communicate()
        if remaining:
            output_lines.append(remaining)
            print(remaining, end="" if remaining.endswith("\n") else "\n", flush=True)
        case.update({"returncode": None, "stdout": "".join(output_lines), "stderr": "", "status": "timeout"})
    else:
        returncode = process.returncode
        case.update(
            {
                "returncode": returncode,
                "stdout": "".join(output_lines),
                "stderr": "",
                "status": "ok" if returncode == 0 else "failed",
            }
        )
    return case


def main() -> int:
    parser = argparse.ArgumentParser(description="Batch-run C trinomial hyperelliptic searches.")
    parser.add_argument("--outdir", default="batch_results_c", help="directory for C output files")
    parser.add_argument("--binary", default="hyperelliptic_finder_c", help="C finder binary path")
    parser.add_argument("--arch", default="-arch x86_64", help="ARCH value passed to make if the binary is missing")
    parser.add_argument("--no-build", action="store_true", help="do not run make if the C binary is missing")
    parser.add_argument("--reduction", choices=("pgl2", "affine", "pgl2save", "affinesave"), default="pgl2save")
    parser.add_argument("--max", type=int, default=0, help="max curves per case; 0 means complete search")
    parser.add_argument(
        "--timeout",
        "--case-timeout",
        dest="timeout",
        type=int,
        default=1800,
        help=(
            "seconds per case before stopping that run after at least one presentation has been saved; "
            "default: 1800; 0 means no timeout"
        ),
    )
    parser.add_argument("--min-genus", type=int, default=DEFAULT_MIN_GENUS)
    parser.add_argument("--max-genus", type=int, default=DEFAULT_MAX_GENUS)
    parser.add_argument("--resume", action="store_true", help="skip cases whose JSON output already exists")
    output_group = parser.add_mutually_exclusive_group()
    output_group.add_argument("--quiet", action="store_true", help="pass --quiet to hyperelliptic_finder_c")
    output_group.add_argument("--verbose", action="store_true", help="pass --verbose to hyperelliptic_finder_c")
    parser.add_argument("--log", action="store_true", help="write batch_c.log.txt")
    parser.add_argument("--monic-only", action="store_true", help="enumerate only monic presentations")
    parser.add_argument(
        "--no-hasse-witt-prefilter",
        action="store_true",
        help="disable the default Hasse-Witt prefilter in hyperelliptic_finder_c",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    binary = Path(args.binary)
    if not binary.is_absolute():
        binary = root / binary
    ensure_binary(root, binary, args.arch, args.no_build)

    outdir = (root / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "batch_summary_c.json"
    original_stdout = sys.stdout
    log_handle = None
    if args.log:
        log_path = outdir / "batch_c.log.txt"
        log_handle = open(log_path, "w", encoding="utf-8", buffering=1)
        sys.stdout = Tee(sys.stdout, log_handle)

    timeout = args.timeout if args.timeout > 0 else None
    summary = []
    try:
        if log_handle is not None:
            print(f"Writing C batch log to {log_path}.", flush=True)
        for p in PRIMES:
            for g in range(args.min_genus, args.max_genus + 1):
                mode_suffix = "_monic" if args.monic_only else ""
                output_json = outdir / f"p{p}_g{g}_{args.reduction}_c{mode_suffix}.json"
                if args.resume and output_json.exists():
                    summary.append({"p": p, "g": g, "runner": "c", "status": "skipped", "output_json": str(output_json)})
                    write_summary(summary_path, summary)
                    print("\n" + SEPARATOR, flush=True)
                    print(f"SKIPPED p={p}, g={g}: existing {output_json}", flush=True)
                    continue
                print_case_header(p, g, args.reduction, output_json)
                case = run_case(
                    binary,
                    outdir,
                    p,
                    g,
                    args.reduction,
                    args.max,
                    timeout,
                    args.quiet,
                    args.verbose,
                    args.monic_only,
                    args.no_hasse_witt_prefilter,
                )
                summary.append(case)
                write_summary(summary_path, summary)
                print_case_footer(case)
    except KeyboardInterrupt:
        summary.append({"runner": "c", "status": "interrupted"})
        write_summary(summary_path, summary)
        print(f"\nInterrupted; wrote partial summary to {summary_path}.")
        sys.stdout = original_stdout
        if log_handle is not None:
            log_handle.close()
        return 130

    print(f"Wrote {summary_path}.")
    sys.stdout = original_stdout
    if log_handle is not None:
        log_handle.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
