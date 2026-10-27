#!/usr/bin/env python3
"""Convert exact physical-series JSON files to Mathematica expression lists.

The default input layout is the one produced by remote_bundle/all22_series:

    INPUT/FAMILY/FAMILY_cK_wW.json

For every selected file the output is:

    OUTPUT/FAMILY/phys_cK_wW_ordN.m

Each Mathematica file is a list in original master order, matching the format
of B1_series_reference.  Algebraic-sector masters are written with half-integer
powers in the positive-root convention encoded by the JSON prefactor metadata.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import os
from pathlib import Path
import sys
import time


FAMILIES = (
    "A", "B1", "B2", "C", "D", "E1", "E2", "F1", "F2", "G1", "G2",
    "G3", "G4", "H", "I1", "I2", "I3", "L3", "M3", "N3", "P3", "R2",
)
CORNERS = ("c1", "c2", "c3")
WEIGHTS = tuple(range(7))
VARIABLES = {"c1": ("x", "y"), "c2": ("y", "z"), "c3": ("x", "z")}
TOKEN_FIELDS = (
    "power_u", "power_v", "log_u", "log_v", "pi", "zeta3", "zeta5",
)
PREFACTORS = {"1", "sqrt_uv"}


def source_path(root: Path, family: str, corner: str, weight: int) -> Path:
    name = f"{family}_{corner}_w{weight}.json"
    nested = root / family / name
    return nested if nested.is_file() else root / name


def power(symbol: str, exponent: int) -> str | None:
    if exponent == 0:
        return None
    return symbol if exponent == 1 else f"{symbol}^{exponent}"


def half_power(symbol: str, integer_part: int) -> str:
    numerator = 2 * integer_part + 1
    return f"Sqrt[{symbol}]" if numerator == 1 else f"{symbol}^({numerator}/2)"


def token_factors(token: tuple[int, ...], variables: tuple[str, str], half: bool) -> list[str]:
    u, v, log_u, log_v, pi, zeta3, zeta5 = token
    result: list[str] = []
    if half:
        result.extend((half_power(variables[0], u), half_power(variables[1], v)))
    else:
        for symbol, exponent in zip(variables, (u, v)):
            if factor := power(symbol, exponent):
                result.append(factor)
    for symbol, exponent in ((f"Log[{variables[0]}]", log_u),
                             (f"Log[{variables[1]}]", log_v),
                             ("Pi", pi), ("Zeta[3]", zeta3), ("Zeta[5]", zeta5)):
        if factor := power(symbol, exponent):
            result.append(factor)
    return result


def term_text(numerator: int, denominator: int, factors: list[str]) -> tuple[bool, str]:
    negative = numerator < 0
    numerator = abs(numerator)
    pieces = ([] if numerator == 1 and factors else [str(numerator)]) + factors
    body = "*".join(pieces)
    if denominator != 1:
        body = f"({body})/{denominator}" if "*" in body else f"{body}/{denominator}"
    return negative, body


class MathematicaWriter:
    def __init__(self, stream, width: int):
        self.stream = stream
        self.width = width
        self.column = 0

    def raw(self, text: str) -> None:
        self.stream.write(text)
        if "\n" in text:
            self.column = len(text.rsplit("\n", 1)[1])
        else:
            self.column += len(text)

    def term(self, separator: str, text: str) -> None:
        piece = separator + text
        if self.column > 1 and self.column + len(piece) > self.width:
            self.raw("\n  " + piece.lstrip())
        else:
            self.raw(piece)


def validate(data: dict, path: Path) -> tuple[str, str, int, int, tuple[str, str], list[tuple[int, ...]], list[str]]:
    try:
        family = data["family"]
        corner = data["corner"]
        weight = data["weight"]
        dimension = data["dimension"]
        variables = tuple(data["variables"])
        bounds = data["truncation"]["max_powers"]
        tokens = [tuple(token) for token in data["tokens"]]
        rows = data["rows"]
    except (KeyError, TypeError) as error:
        raise ValueError(f"{path}: incomplete series metadata") from error

    version = data.get("schema_version")
    if version not in (1, 3):
        raise ValueError(f"{path}: unsupported schema version {version!r}")
    if (family not in FAMILIES or corner not in CORNERS or weight not in WEIGHTS
            or variables != VARIABLES[corner] or data.get("crossing") != ""
            or data.get("coordinates") != "physical"
            or data.get("coefficient_encoding") != "rational"
            or tuple(data.get("token_fields", ())) != TOKEN_FIELDS
            or data.get("index_base") != 0
            or data.get("truncation", {}).get("kind") != "rectangular"
            or not isinstance(dimension, int) or dimension <= 0 or len(rows) != dimension
            or not isinstance(bounds, list) or len(bounds) != 2
            or any(type(n) is not int or n < 0 for n in bounds)):
        raise ValueError(f"{path}: invalid physical-series metadata")

    prefactors = data.get("master_prefactors", ["1"] * dimension)
    if ((version == 1 and "master_prefactors" in data)
            or (version == 3 and "master_prefactors" not in data)
            or len(prefactors) != dimension):
        raise ValueError(f"{path}: invalid master-prefactor metadata")
    if any(prefactor not in PREFACTORS for prefactor in prefactors):
        raise ValueError(f"{path}: unsupported master prefactor")

    if len(tokens) != len(set(tokens)):
        raise ValueError(f"{path}: duplicate tokens")
    for token in tokens:
        if (len(token) != 7 or any(type(n) is not int or n < 0 for n in token)
                or token[0] > bounds[0] or token[1] > bounds[1]
                or token[2] + token[3] + token[4] + 3 * token[5] + 5 * token[6] > weight):
            raise ValueError(f"{path}: invalid token {token!r}")

    return family, corner, weight, dimension, variables, tokens, prefactors


def convert_file(source: Path, destination: Path, width: int, expected: tuple[str, str, int, int]) -> tuple[int, int]:
    with source.open(encoding="utf-8") as stream:
        data = json.load(stream)
    family, corner, weight, dimension, variables, tokens, prefactors = validate(data, source)
    expected_family, expected_corner, expected_weight, expected_order = expected
    if (family, corner, weight) != (expected_family, expected_corner, expected_weight):
        raise ValueError(f"{source}: filename/job metadata mismatch")
    if data["truncation"]["max_powers"] != [expected_order, expected_order]:
        raise ValueError(f"{source}: expected coordinate order {expected_order}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    terms = 0
    try:
        with temporary.open("w", encoding="utf-8") as stream:
            writer = MathematicaWriter(stream, width)
            writer.raw("{")
            for master, row in enumerate(data["rows"]):
                if master:
                    writer.raw(",\n ")
                if not isinstance(row, list):
                    raise ValueError(f"{source}: master {master + 1} is not a term list")
                if not row:
                    writer.raw("0")
                    continue
                seen: set[int] = set()
                half = prefactors[master] == "sqrt_uv"
                for position, encoded in enumerate(row):
                    if not isinstance(encoded, list) or len(encoded) != 3:
                        raise ValueError(f"{source}: malformed term in master {master + 1}")
                    token_id, numerator_text, denominator_text = encoded
                    if (type(token_id) is not int or not 0 <= token_id < len(tokens)
                            or token_id in seen or not isinstance(numerator_text, str)
                            or not isinstance(denominator_text, str)):
                        raise ValueError(f"{source}: invalid term in master {master + 1}")
                    try:
                        coefficient = Fraction(int(numerator_text), int(denominator_text))
                    except (ValueError, ZeroDivisionError) as error:
                        raise ValueError(f"{source}: invalid rational in master {master + 1}") from error
                    if (coefficient == 0 or str(coefficient.numerator) != numerator_text
                            or str(coefficient.denominator) != denominator_text):
                        raise ValueError(f"{source}: noncanonical rational in master {master + 1}")
                    negative, text = term_text(
                        coefficient.numerator,
                        coefficient.denominator,
                        token_factors(tokens[token_id], variables, half),
                    )
                    separator = "-" if position == 0 and negative else "" if position == 0 else " - " if negative else " + "
                    writer.term(separator, text)
                    seen.add(token_id)
                    terms += 1
            writer.raw("}\n")
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return dimension, terms


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="results directory containing family subdirectories")
    parser.add_argument("output", type=Path, help="destination root for Mathematica family directories")
    parser.add_argument("--families", nargs="+", choices=FAMILIES, default=list(FAMILIES))
    parser.add_argument("--corners", nargs="+", choices=CORNERS, default=list(CORNERS))
    parser.add_argument("--weights", nargs="+", type=int, choices=WEIGHTS, default=list(WEIGHTS))
    parser.add_argument("--line-width", type=int, default=100)
    parser.add_argument("--order", type=int, default=20, help="expected equal coordinate order (default: 20)")
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()
    if args.line_width < 40:
        parser.error("--line-width must be at least 40")
    if args.order < 0:
        parser.error("--order must be nonnegative")
    if args.input.resolve() == args.output.resolve():
        parser.error("input and output directories must differ")
    return args


def main() -> int:
    args = parse_args()
    jobs = [(family, corner, weight, source_path(args.input, family, corner, weight))
            for family in args.families for corner in args.corners for weight in args.weights]
    missing = [path for _, _, _, path in jobs if not path.is_file()]
    if missing:
        print(f"error: {len(missing)} expected input file(s) are missing:", file=sys.stderr)
        for path in missing[:20]:
            print(f"  {path}", file=sys.stderr)
        if len(missing) > 20:
            print(f"  ... and {len(missing) - 20} more", file=sys.stderr)
        return 2

    started = time.monotonic()
    total_terms = 0
    for index, (family, corner, weight, source) in enumerate(jobs, 1):
        destination = args.output / family / f"phys_{corner}_w{weight}_ord{args.order}.m"
        if args.skip_existing and destination.is_file():
            print(f"[{index:03d}/{len(jobs):03d}] skip {destination}", flush=True)
            continue
        item_started = time.monotonic()
        masters, terms = convert_file(
            source, destination, args.line_width, (family, corner, weight, args.order)
        )
        total_terms += terms
        print(f"[{index:03d}/{len(jobs):03d}] {family} {corner} w{weight}: "
              f"{masters} masters, {terms} terms, {time.monotonic() - item_started:.1f}s",
              flush=True)
    print(f"Converted {len(jobs)} files with {total_terms} terms in "
          f"{time.monotonic() - started:.1f}s into {args.output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
