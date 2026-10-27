#!/usr/bin/env python3
"""Convert exact trident physical-series JSON files to Mathematica expression lists.

Input layout (series_trident_all_weights):

    INPUT/trident_cK_wW.json

For every selected file the output is:

    OUTPUT/phys_cK_wW_ordN.m

Each Mathematica file is a list of the trident functions in original order, in
the same format produced by convert_series_to_mathematica.py.  Functions with a
sqrt_uv prefactor are written with half-integer powers.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import os
from pathlib import Path
import sys
import time

from convert_series_to_mathematica import (
    CORNERS, PREFACTORS, TOKEN_FIELDS, VARIABLES, WEIGHTS,
    MathematicaWriter, term_text, token_factors,
)


def validate(data: dict, path: Path) -> tuple[str, int, int, tuple[str, str], list[tuple[int, ...]], list[str]]:
    try:
        corner = data["corner"]
        weight = data["weight"]
        dimension = data["dimension"]
        variables = tuple(data["variables"])
        bounds = data["truncation"]["max_powers"]
        tokens = [tuple(token) for token in data["tokens"]]
        rows = data["rows"]
        prefactors = data["function_prefactors"]
    except (KeyError, TypeError) as error:
        raise ValueError(f"{path}: incomplete series metadata") from error

    if (data.get("schema_version") != 1 or data.get("functions") != "published_trident"
            or corner not in CORNERS or weight not in WEIGHTS
            or variables != VARIABLES[corner]
            or data.get("coordinates") != "physical"
            or data.get("coefficient_encoding") != "rational"
            or tuple(data.get("token_fields", ())) != TOKEN_FIELDS
            or data.get("truncation", {}).get("kind") != "rectangular"
            or not isinstance(dimension, int) or dimension <= 0 or len(rows) != dimension
            or not isinstance(bounds, list) or len(bounds) != 2
            or any(type(n) is not int or n < 0 for n in bounds)):
        raise ValueError(f"{path}: invalid trident-series metadata")

    if len(prefactors) != dimension or any(p not in PREFACTORS for p in prefactors):
        raise ValueError(f"{path}: invalid function-prefactor metadata")

    if len(tokens) != len(set(tokens)):
        raise ValueError(f"{path}: duplicate tokens")
    for token in tokens:
        if (len(token) != 7 or any(type(n) is not int or n < 0 for n in token)
                or token[0] > bounds[0] or token[1] > bounds[1]
                or token[2] + token[3] + token[4] + 3 * token[5] + 5 * token[6] > weight):
            raise ValueError(f"{path}: invalid token {token!r}")

    return corner, weight, dimension, variables, tokens, prefactors


def convert_file(source: Path, destination: Path, width: int, expected: tuple[str, int, int]) -> tuple[int, int]:
    with source.open(encoding="utf-8") as stream:
        data = json.load(stream)
    corner, weight, dimension, variables, tokens, prefactors = validate(data, source)
    expected_corner, expected_weight, expected_order = expected
    if (corner, weight) != (expected_corner, expected_weight):
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
            for function, row in enumerate(data["rows"]):
                if function:
                    writer.raw(",\n ")
                if not isinstance(row, list):
                    raise ValueError(f"{source}: function {function + 1} is not a term list")
                if not row:
                    writer.raw("0")
                    continue
                seen: set[int] = set()
                half = prefactors[function] == "sqrt_uv"
                for position, encoded in enumerate(row):
                    if not isinstance(encoded, list) or len(encoded) != 3:
                        raise ValueError(f"{source}: malformed term in function {function + 1}")
                    token_id, numerator_text, denominator_text = encoded
                    if (type(token_id) is not int or not 0 <= token_id < len(tokens)
                            or token_id in seen or not isinstance(numerator_text, str)
                            or not isinstance(denominator_text, str)):
                        raise ValueError(f"{source}: invalid term in function {function + 1}")
                    try:
                        coefficient = Fraction(int(numerator_text), int(denominator_text))
                    except (ValueError, ZeroDivisionError) as error:
                        raise ValueError(f"{source}: invalid rational in function {function + 1}") from error
                    if (coefficient == 0 or str(coefficient.numerator) != numerator_text
                            or str(coefficient.denominator) != denominator_text):
                        raise ValueError(f"{source}: noncanonical rational in function {function + 1}")
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
    parser.add_argument("input", type=Path, help="directory containing trident_cK_wW.json files")
    parser.add_argument("output", type=Path, help="destination directory for the Mathematica files")
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
    jobs = [(corner, weight, args.input / f"trident_{corner}_w{weight}.json")
            for corner in args.corners for weight in args.weights]
    missing = [path for _, _, path in jobs if not path.is_file()]
    if missing:
        print(f"error: {len(missing)} expected input file(s) are missing:", file=sys.stderr)
        for path in missing:
            print(f"  {path}", file=sys.stderr)
        return 2

    started = time.monotonic()
    total_terms = 0
    for index, (corner, weight, source) in enumerate(jobs, 1):
        destination = args.output / f"phys_{corner}_w{weight}_ord{args.order}.m"
        if args.skip_existing and destination.is_file():
            print(f"[{index:02d}/{len(jobs):02d}] skip {destination}", flush=True)
            continue
        item_started = time.monotonic()
        functions, terms = convert_file(
            source, destination, args.line_width, (corner, weight, args.order)
        )
        total_terms += terms
        print(f"[{index:02d}/{len(jobs):02d}] trident {corner} w{weight}: "
              f"{functions} functions, {terms} terms, {time.monotonic() - item_started:.1f}s",
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
