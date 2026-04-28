#!/usr/bin/env python3
"""
Find hyperelliptic curves y^2 = f(x) over F_p whose L-polynomial is trinomial:

    L(t) = 1 + a_g t^g + p^g t^{2g},

with a_g nonzero.
"""

from __future__ import annotations

import argparse
import atexit
import json
import math
import os
import signal
import sys
from dataclasses import asdict, dataclass
from functools import lru_cache
from itertools import product
from typing import Iterable, Sequence


DEFAULT_OUTPUT = "hyperelliptic_results.txt"
results: list["CurveResult"] = []
run_context: dict[str, object] = {}
log_handle = None


@dataclass
class CurveResult:
    index: int
    f_coeffs: list[int]
    f_polynomial: str
    middle_coefficient: int


@dataclass
class SearchStats:
    considered: int = 0
    skipped_by_reduction: int = 0
    checked: int = 0
    rejected_by_early_l_coefficient: int = 0
    rejected_by_zero_middle_coefficient: int = 0
    saved: int = 0


def output_paths(path: str) -> tuple[str, str]:
    if path.endswith(".json"):
        return path[:-5] + ".txt", path
    if path.endswith(".txt"):
        return path, path[:-4] + ".json"
    return path, path + ".json"


def log_path_for_output(path: str) -> str:
    if path.endswith(".txt"):
        return path[:-4] + ".log.txt"
    return path + ".log.txt"


def open_log(path: str) -> None:
    global log_handle
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    log_handle = open(path, "w", encoding="utf-8", buffering=1)


def close_log() -> None:
    global log_handle
    if log_handle is not None:
        log_handle.close()
        log_handle = None


def emit(message: str = "") -> None:
    print(message, flush=True)
    if log_handle is not None:
        log_handle.write(message + "\n")


def trim(poly: Sequence[int]) -> list[int]:
    out = list(poly)
    while out and out[-1] == 0:
        out.pop()
    return out


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, math.isqrt(n) + 1, 2):
        if n % i == 0:
            return False
    return True


def mod_inv(a: int, p: int) -> int:
    a %= p
    if a == 0:
        raise ZeroDivisionError("0 has no inverse")
    return pow(a, -1, p)


def poly_add(a: Sequence[int], b: Sequence[int], p: int) -> list[int]:
    length = max(len(a), len(b))
    return trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(length)])


def poly_sub(a: Sequence[int], b: Sequence[int], p: int) -> list[int]:
    length = max(len(a), len(b))
    return trim([((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p for i in range(length)])


def poly_mul(a: Sequence[int], b: Sequence[int], p: int) -> list[int]:
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        for j, cb in enumerate(b):
            out[i + j] = (out[i + j] + ca * cb) % p
    return trim(out)


def poly_divmod(a: Sequence[int], b: Sequence[int], p: int) -> tuple[list[int], list[int]]:
    dividend = trim([c % p for c in a])
    divisor = trim([c % p for c in b])
    if not divisor:
        raise ZeroDivisionError("polynomial division by zero")
    if len(dividend) < len(divisor):
        return [], dividend

    quotient = [0] * (len(dividend) - len(divisor) + 1)
    inv_lead = mod_inv(divisor[-1], p)
    while len(dividend) >= len(divisor) and dividend:
        shift = len(dividend) - len(divisor)
        coeff = dividend[-1] * inv_lead % p
        quotient[shift] = coeff
        for i, c in enumerate(divisor):
            dividend[shift + i] = (dividend[shift + i] - coeff * c) % p
        dividend = trim(dividend)
    return trim(quotient), dividend


def poly_mod(a: Sequence[int], modulus: Sequence[int], p: int) -> list[int]:
    return poly_divmod(a, modulus, p)[1]


def poly_gcd(a: Sequence[int], b: Sequence[int], p: int) -> list[int]:
    x = trim([c % p for c in a])
    y = trim([c % p for c in b])
    while y:
        x, y = y, poly_mod(x, y, p)
    if not x:
        return []
    return [(c * mod_inv(x[-1], p)) % p for c in x]


def poly_derivative(poly: Sequence[int], p: int) -> list[int]:
    return trim([(i * poly[i]) % p for i in range(1, len(poly))])


def poly_squarefree(poly: Sequence[int], p: int) -> bool:
    deriv = poly_derivative(poly, p)
    if not deriv:
        return False
    return poly_gcd(poly, deriv, p) == [1]


def poly_pow_mod(base: Sequence[int], exponent: int, modulus: Sequence[int], p: int) -> list[int]:
    result = [1]
    power = poly_mod(base, modulus, p)
    while exponent:
        if exponent & 1:
            result = poly_mod(poly_mul(result, power, p), modulus, p)
        power = poly_mod(poly_mul(power, power, p), modulus, p)
        exponent >>= 1
    return result


def is_irreducible(poly: Sequence[int], p: int) -> bool:
    f = trim(poly)
    degree = len(f) - 1
    if degree <= 0 or f[-1] % p == 0:
        return False
    if degree == 1:
        return True

    x = [0, 1]
    power = x[:]
    for _ in range(1, degree // 2 + 1):
        power = poly_pow_mod(power, p, f, p)
        if poly_gcd(poly_sub(power, x, p), f, p) != [1]:
            return False
    return poly_sub(poly_pow_mod(x, p**degree, f, p), x, p) == []


@lru_cache(maxsize=None)
def find_irreducible_polynomial(p: int, degree: int) -> tuple[int, ...]:
    if degree == 1:
        return (0, 1)

    known = {
        (2, 2): (1, 1, 1),
        (2, 3): (1, 1, 0, 1),
        (2, 4): (1, 1, 0, 0, 1),
        (3, 2): (2, 0, 1),
        (5, 2): (2, 0, 1),
    }
    if (p, degree) in known and is_irreducible(known[(p, degree)], p):
        return known[(p, degree)]

    for coeffs in product(range(p), repeat=degree):
        candidate = tuple(coeffs) + (1,)
        if is_irreducible(candidate, p):
            return candidate
    raise ValueError(f"could not find an irreducible polynomial of degree {degree} over F_{p}")


class FiniteField:
    def __init__(self, p: int, degree: int):
        self.p = p
        self.degree = degree
        self.order = p**degree
        self.modulus = list(find_irreducible_polynomial(p, degree))

    def coeffs(self, value: int) -> list[int]:
        out = []
        x = value
        for _ in range(self.degree):
            out.append(x % self.p)
            x //= self.p
        return out

    def element(self, coeffs: Sequence[int]) -> int:
        value = 0
        multiplier = 1
        for c in coeffs[: self.degree]:
            value += (c % self.p) * multiplier
            multiplier *= self.p
        return value

    def add(self, a: int, b: int) -> int:
        return self.element(poly_add(self.coeffs(a), self.coeffs(b), self.p))

    def mul(self, a: int, b: int) -> int:
        reduced = poly_mod(poly_mul(self.coeffs(a), self.coeffs(b), self.p), self.modulus, self.p)
        return self.element(reduced)

    def pow(self, a: int, exponent: int) -> int:
        result = 1
        power = a
        while exponent:
            if exponent & 1:
                result = self.mul(result, power)
            power = self.mul(power, power)
            exponent >>= 1
        return result

    def __str__(self) -> str:
        if self.degree == 1:
            return f"F_{self.p}"
        return f"F_{{{self.p}^{self.degree}}} = F_{self.p}[u]/({format_polynomial(self.modulus, 'u')})"


@lru_cache(maxsize=None)
def get_field(p: int, degree: int) -> FiniteField:
    return FiniteField(p, degree)


@lru_cache(maxsize=None)
def square_solution_counts(p: int, degree: int) -> tuple[int, ...]:
    field = get_field(p, degree)
    counts = [0] * field.order
    for y in range(field.order):
        counts[field.mul(y, y)] += 1
    return tuple(counts)


def eval_polynomial_over_field(f_coeffs: Sequence[int], x: int, field: FiniteField) -> int:
    result = 0
    for coeff in reversed(f_coeffs):
        result = field.add(field.mul(result, x), coeff % field.p)
    return result


def count_points(p: int, k: int, f_coeffs: Sequence[int]) -> int:
    field = get_field(p, k)
    counts = square_solution_counts(p, k)
    finite = 0
    for x in range(field.order):
        finite += counts[eval_polynomial_over_field(f_coeffs, x, field)]

    degree = len(trim(f_coeffs)) - 1
    leading_coeff = f_coeffs[degree] % p
    points_at_infinity = 1 if degree % 2 == 1 else counts[leading_coeff]
    return finite + points_at_infinity


def next_l_coefficient(coeffs: Sequence[int], power_sums: Sequence[int], k: int) -> int:
    numerator = power_sums[k - 1]
    for i in range(1, k):
        numerator += coeffs[i] * power_sums[k - i - 1]
    if numerator % k != 0:
        raise ArithmeticError(f"Newton identity produced nonintegral coefficient for k={k}")
    return -numerator // k


def trinomial_middle_coefficient(p: int, g: int, f_coeffs: Sequence[int]) -> tuple[str, int | None]:
    coeffs = [1]
    power_sums = []
    for k in range(1, g + 1):
        power_sums.append(p**k + 1 - count_points(p, k, f_coeffs))
        coeff = next_l_coefficient(coeffs, power_sums, k)
        if k < g and coeff != 0:
            return "early_l_coefficient", None
        coeffs.append(coeff)

    if coeffs[g] == 0:
        return "zero_middle_coefficient", None
    return "valid", coeffs[g]


def poly_compose_linear(poly: Sequence[int], scale: int, shift: int, p: int) -> list[int]:
    result: list[int] = []
    linear = [shift % p, scale % p]
    for coeff in reversed(poly):
        result = poly_add(poly_mul(result, linear, p), [coeff % p], p)
    return result


def has_square_root(value: int, p: int) -> bool:
    value %= p
    if value == 0 or p == 2:
        return True
    return pow(value, (p - 1) // 2, p) == 1


def affine_normalized_transform(poly: Sequence[int], scale: int, shift: int, p: int) -> tuple[int, ...] | None:
    transformed = trim(poly_compose_linear(poly, scale, shift, p))
    if len(transformed) != len(poly):
        return None

    leading = transformed[-1] % p
    if not has_square_root(leading, p):
        return None

    inv_leading = mod_inv(leading, p)
    return tuple((coeff * inv_leading) % p for coeff in transformed)


@lru_cache(maxsize=None)
def pgl2_matrices(p: int) -> tuple[tuple[int, int, int, int], ...]:
    matrices = set()
    for a, b, c, d in product(range(p), repeat=4):
        determinant = (a * d - b * c) % p
        if determinant == 0:
            continue
        entries = [a, b, c, d]
        first_nonzero = next(entry for entry in entries if entry != 0)
        inv = mod_inv(first_nonzero, p)
        matrices.add(tuple((entry * inv) % p for entry in entries))
    return tuple(sorted(matrices))


def binary_form_pgl2_transform(
    poly: Sequence[int],
    genus: int,
    matrix: tuple[int, int, int, int],
    p: int,
) -> tuple[int, ...] | None:
    a, b, c, d = matrix
    binary_degree = 2 * genus + 2
    if len(poly) > binary_degree + 1:
        return None

    first_linear = [b % p, a % p]
    second_linear = [d % p, c % p]
    first_powers = [[1]]
    second_powers = [[1]]
    for _ in range(binary_degree):
        first_powers.append(poly_mul(first_powers[-1], first_linear, p))
        second_powers.append(poly_mul(second_powers[-1], second_linear, p))

    transformed: list[int] = []
    for i, coeff in enumerate(poly):
        if coeff % p == 0:
            continue
        term = poly_mul(first_powers[i], second_powers[binary_degree - i], p)
        term = [(coeff * term_coeff) % p for term_coeff in term]
        transformed = poly_add(transformed, term, p)

    transformed = trim(transformed)
    degree = len(transformed) - 1
    if degree not in (2 * genus + 1, 2 * genus + 2):
        return None

    leading = transformed[-1] % p
    if not has_square_root(leading, p):
        return None

    inv_leading = mod_inv(leading, p)
    return tuple((coeff * inv_leading) % p for coeff in transformed)


def canonical_key(poly: Sequence[int], p: int, genus: int, reduction: str) -> tuple[int, ...]:
    transforms = []
    if reduction == "affine":
        for scale in range(1, p):
            for shift in range(p):
                transformed = affine_normalized_transform(poly, scale, shift, p)
                if transformed is not None:
                    transforms.append(transformed)
    elif reduction == "pgl2":
        for matrix in pgl2_matrices(p):
            transformed = binary_form_pgl2_transform(poly, genus, matrix, p)
            if transformed is not None:
                transforms.append(transformed)
    else:
        raise ValueError(f"unknown reduction mode: {reduction}")
    return min(transforms)


def generate_polynomials(p: int, g: int) -> Iterable[list[int]]:
    for degree in (2 * g + 1, 2 * g + 2):
        for lower_coeffs in product(range(p), repeat=degree):
            poly = list(lower_coeffs) + [1]
            if poly_squarefree(poly, p):
                yield poly


@lru_cache(maxsize=None)
def presentation_number(target: tuple[int, ...], p: int, g: int) -> int | None:
    for number, poly in enumerate(generate_polynomials(p, g), 1):
        if tuple(poly) == target:
            return number
    return None


def format_polynomial(coeffs: Sequence[int], var: str = "x") -> str:
    terms = []
    for exponent, coeff in enumerate(coeffs):
        if coeff == 0:
            continue
        if exponent == 0:
            terms.append(str(coeff))
        elif exponent == 1:
            terms.append(var if coeff == 1 else f"{coeff}{var}")
        else:
            terms.append(f"{var}^{exponent}" if coeff == 1 else f"{coeff}{var}^{exponent}")
    return " + ".join(terms) if terms else "0"


def write_results() -> None:
    if not run_context:
        return

    output_path = str(run_context["output_file"])
    json_path = str(run_context["json_file"])
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(json_path)), exist_ok=True)

    p = int(run_context["p"])
    g = int(run_context["g"])

    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("Hyperelliptic curves with trinomial L-polynomial\n")
        handle.write("=" * 56 + "\n")
        handle.write(f"p = {p}\n")
        handle.write(f"g = {g}\n")
        handle.write(f"base field: F_{p}\n")
        handle.write(f"presentation reduction: {run_context['reduction']}\n")
        handle.write(f"log file: {run_context.get('log_file', '')}\n")
        handle.write(f"search status: {run_context.get('search_status', 'incomplete')}\n")
        handle.write(f"complete list: {run_context.get('complete_list', False)}\n")
        handle.write(f"total presentations found: {len(results)}\n\n")
        handle.write("Only the middle L-polynomial coefficient a_g is listed.\n")
        handle.write("Each entry is a presentation y^2 = f(x).\n\n")
        for result in results:
            handle.write(f"[{result.index}]\n")
            handle.write(f"f(x) = {result.f_polynomial}\n")
            handle.write(f"f_coeffs = {result.f_coeffs}\n")
            handle.write(f"middle_coefficient_a_{g} = {result.middle_coefficient}\n\n")

    data = {
        "p": p,
        "g": g,
        "field": f"F_{p}",
        "presentation_reduction": str(run_context["reduction"]),
        "log_file": str(run_context.get("log_file", "")),
        "search_status": str(run_context.get("search_status", "incomplete")),
        "complete_list": bool(run_context.get("complete_list", False)),
        "total_presentations_found": len(results),
        "curves": [asdict(result) for result in results],
    }
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)
        handle.write("\n")


def find_curves(p: int, g: int, max_curves: int, reduction: str, verbose: bool) -> list[CurveResult]:
    stats = SearchStats()
    canonical_owner: dict[tuple[int, ...], int] = {}
    for f_coeffs in generate_polynomials(p, g):
        stats.considered += 1
        if verbose:
            emit(f"[{stats.considered}]")
            emit(f"Considering f(x) = {format_polynomial(f_coeffs)}")
            emit(f"Checking {reduction} repeats.")

        key = canonical_key(f_coeffs, p, g, reduction)
        owner = canonical_owner.get(key)
        if tuple(f_coeffs) != key:
            stats.skipped_by_reduction += 1
            if verbose:
                if owner is None:
                    owner = presentation_number(key, p, g)
                    canonical_owner[key] = owner if owner is not None else stats.considered
                if owner is None:
                    emit(f"Skipped: {reduction}-equivalent repeat.")
                else:
                    emit(f"Skipped: {reduction}-equivalent repeat of [{owner}].")
            continue
        canonical_owner.setdefault(key, stats.considered)

        stats.checked += 1
        if verbose:
            emit("Canonical; checking L-polynomial coefficients.")
        status, middle_coefficient = trinomial_middle_coefficient(p, g, f_coeffs)
        if status == "early_l_coefficient":
            stats.rejected_by_early_l_coefficient += 1
            if verbose:
                emit("Early rejected: pre-middle coefficient is nonzero.")
            continue
        if status == "zero_middle_coefficient":
            stats.rejected_by_zero_middle_coefficient += 1
            if verbose:
                emit("Rejected: middle coefficient is zero.")
            continue
        if middle_coefficient is None:
            raise AssertionError("valid trinomial candidate is missing its middle coefficient")

        result = CurveResult(
            index=len(results) + 1,
            f_coeffs=f_coeffs,
            f_polynomial=format_polynomial(f_coeffs),
            middle_coefficient=middle_coefficient,
        )
        results.append(result)
        stats.saved += 1
        write_results()
        if verbose:
            emit(f"Saved presentation {result.index}; a_{g} = {result.middle_coefficient}.")
        else:
            emit(f"Saved presentation {result.index}; a_{g} = {result.middle_coefficient}.")

        if max_curves > 0 and len(results) >= max_curves:
            break

    emit("")
    emit("SUMMARY")
    emit(f"Considered {stats.considered} squarefree monic presentations.")
    emit(f"Checked {stats.checked} canonical presentations after {reduction} reduction.")
    emit(f"Skipped {stats.skipped_by_reduction} {reduction}-equivalent presentations.")
    emit(f"Early-rejected {stats.rejected_by_early_l_coefficient} presentations.")
    emit(f"Rejected {stats.rejected_by_zero_middle_coefficient} presentations with zero middle coefficient.")
    return results


def handle_interrupt(_sig, _frame) -> None:
    run_context["search_status"] = "interrupted"
    run_context["complete_list"] = False
    write_results()
    output_file = run_context.get("output_file", DEFAULT_OUTPUT)
    emit("")
    emit(f"Interrupted. Saved {len(results)} presentations to {output_file}.")
    raise SystemExit(130)


def save_on_exit() -> None:
    write_results()


def main() -> int:
    parser = argparse.ArgumentParser(description="Find hyperelliptic curves over F_p with trinomial L-polynomials")
    parser.add_argument("p", type=int, help="prime characteristic p")
    parser.add_argument("g", type=int, help="genus g")
    parser.add_argument("--max", type=int, default=0, help="maximum number of presentations to save; 0 means all")
    parser.add_argument("--output", type=str, default=DEFAULT_OUTPUT, help="text or JSON output path")
    parser.add_argument(
        "--reduction",
        choices=("pgl2", "affine"),
        default="pgl2",
        help="presentation reduction to use before point counting; default: pgl2",
    )
    parser.add_argument("--quiet", action="store_true", help="only print saved matches and final summary")
    args = parser.parse_args()

    if not is_prime(args.p):
        print(f"Error: {args.p} is not prime", file=sys.stderr)
        return 1
    if args.p == 2:
        print("Error: characteristic 2 is not supported for models y^2 = f(x); use an odd prime p", file=sys.stderr)
        return 1
    if args.g < 1:
        print("Error: g must be at least 1", file=sys.stderr)
        return 1
    if args.max < 0:
        print("Error: --max must be nonnegative", file=sys.stderr)
        return 1

    text_path, json_path = output_paths(args.output)
    log_path = log_path_for_output(text_path)
    open_log(log_path)
    run_context.update(
        {
            "p": args.p,
            "g": args.g,
            "output_file": text_path,
            "json_file": json_path,
            "log_file": log_path,
            "reduction": args.reduction,
            "search_status": "incomplete",
            "complete_list": False,
        }
    )
    atexit.register(save_on_exit)
    atexit.register(close_log)
    signal.signal(signal.SIGINT, handle_interrupt)
    signal.signal(signal.SIGTERM, handle_interrupt)

    write_results()
    emit(f"Searching over F_{args.p}.")
    emit(f"Using {args.reduction} presentation reduction.")
    emit(f"Saving detailed results to {text_path} and {json_path}.")
    emit(f"Writing progress log to {log_path}.")
    try:
        find_curves(args.p, args.g, args.max, args.reduction, verbose=not args.quiet)
        run_context["search_status"] = "complete" if args.max == 0 else "max_reached"
        run_context["complete_list"] = args.max == 0
        return_code = 0
    except Exception:
        run_context["search_status"] = "error"
        run_context["complete_list"] = False
        raise
    finally:
        write_results()

    emit(f"Done. Saved {len(results)} presentations to {text_path} and {json_path}.")
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
