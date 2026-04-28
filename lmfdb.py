#!/usr/bin/env python3
"""Compare local trinomial L-polynomial data with available LMFDB records.

LMFDB's finite-field data is organized by abelian-variety isogeny class, not by
our curve presentations. This script checks whether each local trinomial
L-polynomial corresponds to an available LMFDB isogeny class and, when present,
whether LMFDB records hyperelliptic Jacobians for that class.
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


LMFDB_BASE = "https://www.lmfdb.org"


def encode_lmfdb_integer(value: int) -> str:
    """Encode an integer as used in LMFDB abelian-variety finite-field labels."""
    if value == 0:
        return "a"
    sign = "a" if value < 0 else ""
    n = abs(value)
    digits = []
    while n:
        digits.append(chr(ord("a") + (n % 26)))
        n //= 26
    return sign + "".join(reversed(digits))


def lmfdb_label(p: int, g: int, middle_coefficient: int) -> str:
    coeffs = [0] * g
    coeffs[g - 1] = middle_coefficient
    iso = "_".join(encode_lmfdb_integer(coeff) for coeff in coeffs)
    return f"{g}.{p}.{iso}"


def load_json_url(url: str, timeout: int) -> dict[str, Any] | None:
    request = urllib.request.Request(url, headers={"User-Agent": "hyperelliptic-finder/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            data = response.read().decode("utf-8")
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
        return None

    try:
        return json.loads(data)
    except json.JSONDecodeError:
        return None


def fetch_lmfdb_record(label: str, timeout: int) -> tuple[str, dict[str, Any] | None, str | None]:
    quoted = urllib.parse.quote(label)
    urls = [
        f"{LMFDB_BASE}/api/av_fq_isog/?label={quoted}&_format=json",
        f"{LMFDB_BASE}/Variety/Abelian/Fq/data/{quoted}?_format=json",
    ]
    for url in urls:
        payload = load_json_url(url, timeout)
        if not payload:
            continue
        if "data" in payload and isinstance(payload["data"], list) and payload["data"]:
            data = payload["data"][0]
            if isinstance(data, list) and data and isinstance(data[0], dict):
                return "found", data[0], url
            if isinstance(data, dict):
                return "found", data, url
        if "av_fq_isog" in payload:
            record = payload["av_fq_isog"]
            if isinstance(record, list):
                record = record[0] if record else None
            if isinstance(record, dict):
                return "found", record, url
        if payload.get("label") == label:
            return "found", payload, url
    return "missing_or_unavailable", None, None


def local_l_polynomial_key(curve: dict[str, Any], p: int, g: int) -> tuple[int, int, int]:
    return (p, g, int(curve["middle_coefficient"]))


def is_result_file(data: dict[str, Any]) -> bool:
    return "p" in data and "g" in data and isinstance(data.get("curves"), list)


def compare_file(path: Path, timeout: int, delay: float) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        local = json.load(handle)

    if not is_result_file(local):
        return {
            "local_file": str(path),
            "status": "skipped",
            "reason": "not a hyperelliptic_finder result JSON file",
        }

    p = int(local["p"])
    g = int(local["g"])
    unique = {}
    for curve in local.get("curves", []):
        unique[local_l_polynomial_key(curve, p, g)] = curve

    comparisons = []
    for (_p, _g, middle), curve in sorted(unique.items()):
        label = lmfdb_label(p, g, middle)
        status, record, url = fetch_lmfdb_record(label, timeout)
        if delay:
            time.sleep(delay)

        hyp_count = None
        curve_count = None
        jacobian_count = None
        if record:
            hyp_count = record.get("hyp_count")
            curve_count = record.get("curve_count")
            jacobian_count = record.get("jacobian_count")

        local_curves = [
            item
            for item in local.get("curves", [])
            if int(item["middle_coefficient"]) == middle
        ]
        canonical_indices = sorted(
            {
                int(item.get("canonical_presentation_index", item["index"]))
                for item in local_curves
            }
        )
        comparisons.append(
            {
                "middle_coefficient": middle,
                "local_presentation_count": len(local_curves),
                "local_curve_indices": [item["index"] for item in local_curves],
                "local_canonical_presentation_indices": canonical_indices,
                "local_reused_l_polynomial_count": sum(
                    1 for item in local_curves if bool(item.get("reused_l_polynomial", False))
                ),
                "lmfdb_label": label,
                "lmfdb_status": status,
                "lmfdb_url": url,
                "lmfdb_hyp_count": hyp_count,
                "lmfdb_curve_count": curve_count,
                "lmfdb_jacobian_count": jacobian_count,
                "matches_hyperelliptic_available": bool(hyp_count and int(hyp_count) > 0),
            }
        )

    available = [row for row in comparisons if row["lmfdb_status"] == "found"]
    hyperelliptic_matches = [row for row in comparisons if row["matches_hyperelliptic_available"]]
    return {
        "local_file": str(path),
        "status": "compared",
        "p": p,
        "g": g,
        "presentation_reduction": local.get("presentation_reduction"),
        "search_status": local.get("search_status"),
        "complete_list": local.get("complete_list"),
        "local_presentations": len(local.get("curves", [])),
        "local_unique_l_polynomials": len(comparisons),
        "lmfdb_records_found": len(available),
        "lmfdb_hyperelliptic_matches": len(hyperelliptic_matches),
        "accuracy_among_available": (
            len(hyperelliptic_matches) / len(available) if available else None
        ),
        "note": (
            "LMFDB comparison is by abelian-variety isogeny class. Saved equivalent "
            "presentations are grouped by middle coefficient before querying LMFDB."
        ),
        "comparisons": comparisons,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare local results with available LMFDB data.")
    parser.add_argument("files", nargs="+", help="local result JSON files, e.g. p3_g2_results.json")
    parser.add_argument("--out", default="lmfdb_accuracy.json", help="output JSON report")
    parser.add_argument("--timeout", type=int, default=20, help="HTTP timeout in seconds")
    parser.add_argument("--delay", type=float, default=0.2, help="delay between LMFDB requests")
    args = parser.parse_args()

    reports = [compare_file(Path(path), args.timeout, args.delay) for path in args.files]
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(reports, handle, indent=2)
        handle.write("\n")

    for report in reports:
        if report.get("status") == "skipped":
            print(f"{report['local_file']}: skipped ({report['reason']})")
            continue
        total = report["local_unique_l_polynomials"]
        found = report["lmfdb_records_found"]
        matched = report["lmfdb_hyperelliptic_matches"]
        presentations = report["local_presentations"]
        accuracy = report["accuracy_among_available"]
        accuracy_text = "n/a" if accuracy is None else f"{accuracy:.3f}"
        print(
            f"{report['local_file']}: presentations={presentations}, unique={total}, lmfdb_found={found}, "
            f"hyperelliptic_matches={matched}, accuracy_among_available={accuracy_text}"
        )
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
