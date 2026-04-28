#!/usr/bin/env python3
"""Batch-run hyperelliptic_finder.py for primes <= 11 and genera <= 20."""

from __future__ import annotations

import argparse
import json
import select
import subprocess
import sys
import time
from pathlib import Path


PRIMES = (3, 5, 7)
DEFAULT_MAX_GENUS = 5
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
    print(f"CASE p={p}, g={g}, reduction={reduction}", flush=True)
    print(f"output: {output_json}", flush=True)


def print_case_footer(case: dict[str, object]) -> None:
    stderr = str(case.get("stderr") or "")
    if stderr:
        print("stderr:", flush=True)
        print(stderr, end="" if stderr.endswith("\n") else "\n", flush=True)


def run_case(
    script: Path,
    outdir: Path,
    p: int,
    g: int,
    reduction: str,
    max_curves: int,
    timeout: int | None,
    quiet: bool,
) -> dict[str, object]:
    output_base = outdir / f"p{p}_g{g}_{reduction}.txt"
    command = [
        sys.executable,
        str(script),
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

    case = {
        "p": p,
        "g": g,
        "reduction": reduction,
        "output_txt": str(output_base),
        "output_json": str(output_base.with_suffix(".json")),
        "command": command,
    }
    output_lines: list[str] = []
    start = time.monotonic()
    process = subprocess.Popen(
        command,
        cwd=script.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
    )
    assert process.stdout is not None

    timed_out = False
    while True:
        if timeout is not None and time.monotonic() - start > timeout:
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
        case.update(
            {
                "returncode": None,
                "stdout": "".join(output_lines),
                "stderr": "",
                "status": "timeout",
            }
        )
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
    parser = argparse.ArgumentParser(description="Batch-run trinomial hyperelliptic searches.")
    parser.add_argument("--outdir", default="batch_results", help="directory for output files")
    parser.add_argument("--reduction", choices=("pgl2", "affine"), default="pgl2")
    parser.add_argument("--max", type=int, default=0, help="max curves per case; 0 means complete search")
    parser.add_argument("--timeout", type=int, default=0, help="seconds per case; 0 means no timeout")
    parser.add_argument("--min-genus", type=int, default=1)
    parser.add_argument("--max-genus", type=int, default=DEFAULT_MAX_GENUS)
    parser.add_argument("--resume", action="store_true", help="skip cases whose JSON output already exists")
    parser.add_argument("--quiet", action="store_true", help="pass --quiet to hyperelliptic_finder.py")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    script = root / "hyperelliptic_finder.py"
    outdir = (root / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "batch_summary.json"
    log_path = outdir / "batch.log.txt"
    log_handle = open(log_path, "w", encoding="utf-8", buffering=1)
    original_stdout = sys.stdout
    sys.stdout = Tee(sys.stdout, log_handle)

    timeout = args.timeout if args.timeout > 0 else None
    summary = []
    try:
        print(f"Writing batch log to {log_path}.", flush=True)
        for p in PRIMES:
            for g in range(args.min_genus, args.max_genus + 1):
                output_json = outdir / f"p{p}_g{g}_{args.reduction}.json"
                if args.resume and output_json.exists():
                    summary.append({"p": p, "g": g, "status": "skipped", "output_json": str(output_json)})
                    write_summary(summary_path, summary)
                    print("\n" + SEPARATOR, flush=True)
                    print(f"SKIPPED p={p}, g={g}: existing {output_json}", flush=True)
                    continue
                print_case_header(p, g, args.reduction, output_json)
                case = run_case(script, outdir, p, g, args.reduction, args.max, timeout, args.quiet)
                summary.append(case)
                write_summary(summary_path, summary)
                print_case_footer(case)
    except KeyboardInterrupt:
        summary.append({"status": "interrupted"})
        write_summary(summary_path, summary)
        print(f"\nInterrupted; wrote partial summary to {summary_path}.")
        sys.stdout = original_stdout
        log_handle.close()
        return 130

    print(f"Wrote {summary_path}.")
    sys.stdout = original_stdout
    log_handle.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
