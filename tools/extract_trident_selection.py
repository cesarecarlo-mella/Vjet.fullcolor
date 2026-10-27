#!/usr/bin/env python3
"""Extract a selection of Trident[index, weight] series into one small file per corner.

Reads the expression lists written by convert_trident_to_mathematica.py:

    INPUT/phys_cK_wW_ordN.m

and a Mathematica list of Trident[index, weight] heads (index 1-based).  For each
corner the output is

    OUTPUT/trident_selected_cK_ordN.m

containing an Association <|Trident[i, w] -> series, ...|> in selection order.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import os
from pathlib import Path
import re
import sys
import time

CORNERS = ("c1", "c2", "c3")
DIMENSION = 1318
SELECTION = re.compile(r"Trident\[\s*(\d+)\s*,\s*(\d+)\s*\]")


def read_selection(path: Path) -> list[tuple[int, int]]:
    selection = [(int(i), int(w)) for i, w in SELECTION.findall(path.read_text(encoding="utf-8"))]
    if not selection:
        raise ValueError(f"{path}: no Trident[index, weight] entries found")
    if len(selection) != len(set(selection)):
        raise ValueError(f"{path}: duplicate entries")
    for index, weight in selection:
        if not 1 <= index <= DIMENSION or not 0 <= weight <= 6:
            raise ValueError(f"{path}: entry Trident[{index}, {weight}] out of range")
    return selection


def read_functions(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8").strip()
    if not (text.startswith("{") and text.endswith("}")):
        raise ValueError(f"{path}: not a Mathematica list")
    # Series terms contain no commas, so top-level entries are split by ",\n ".
    functions = text[1:-1].split(",\n ")
    if len(functions) != DIMENSION:
        raise ValueError(f"{path}: expected {DIMENSION} functions, found {len(functions)}")
    return functions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selection", type=Path, help="file with a list of Trident[index, weight]")
    parser.add_argument("input", type=Path, help="directory with phys_cK_wW_ordN.m files")
    parser.add_argument("output", type=Path, help="destination directory")
    parser.add_argument("--corners", nargs="+", choices=CORNERS, default=list(CORNERS))
    parser.add_argument("--order", type=int, default=20)
    args = parser.parse_args()

    selection = read_selection(args.selection)
    by_weight: dict[int, list[int]] = defaultdict(list)
    for index, weight in selection:
        by_weight[weight].append(index)

    args.output.mkdir(parents=True, exist_ok=True)
    for corner in args.corners:
        started = time.monotonic()
        series: dict[tuple[int, int], str] = {}
        for weight, indices in sorted(by_weight.items()):
            functions = read_functions(args.input / f"phys_{corner}_w{weight}_ord{args.order}.m")
            for index in indices:
                series[index, weight] = functions[index - 1]
            del functions
        destination = args.output / f"trident_selected_{corner}_ord{args.order}.m"
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        try:
            with temporary.open("w", encoding="utf-8") as stream:
                stream.write("<|")
                for position, key in enumerate(selection):
                    if position:
                        stream.write(",\n ")
                    stream.write(f"Trident[{key[0]}, {key[1]}] -> {series[key]}")
                stream.write("|>\n")
            os.replace(temporary, destination)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
        zeros = sum(series[key] == "0" for key in selection)
        print(f"{corner}: {len(selection)} entries ({zeros} identically zero) -> {destination}, "
              f"{time.monotonic() - started:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
