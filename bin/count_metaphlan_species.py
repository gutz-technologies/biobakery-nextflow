#!/usr/bin/env python3
"""Count detected species per sample in a merged MetaPhlAn profile."""

import argparse
import csv
import math
from pathlib import Path
from typing import List, Tuple


def count_species(source: Path, index: str) -> List[Tuple[str, int]]:
    """Read a merged profile and count its positive terminal species rows.

    Parameters
    ----------
    source : Path
        Merged MetaPhlAn TSV, optionally preceded by metadata comments.
    index : str
        Database index suffix to remove from the merged sample names.

    Returns
    -------
    list of tuple
        Sample IDs and counts in the order of the table's columns.

    Raises
    ------
    ValueError
        If the table has invalid headers, dimensions, or abundances.
    """
    suffix = "_" + index
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for header in reader:
            if not header:
                continue
            if header[0] in ("clade_name", "#clade_name"):
                break
            if not header[0].startswith("#"):
                raise ValueError("expected a clade_name header")
        else:
            raise ValueError("missing clade_name header")

        samples = [name[:-len(suffix)] if name.endswith(suffix) else name
                   for name in header[1:]]
        if not samples or not all(name.strip() for name in samples):
            raise ValueError("expected nonempty sample names")
        if len(set(samples)) != len(samples):
            raise ValueError("duplicate sample names after removing the index suffix")
        counts = [0] * len(samples)

        for row in reader:
            if not row:
                continue
            if len(row) != len(header) or not row[0]:
                raise ValueError(f"line {reader.line_num}: incorrect table row")
            terminal_rank = row[0].rsplit("|", 1)[-1]
            is_species = terminal_rank.startswith("s__") and len(terminal_rank) > 3
            for i, text in enumerate(row[1:]):
                try:
                    abundance = float(text)
                except ValueError as exc:
                    raise ValueError(f"line {reader.line_num}: non-numeric abundance") from exc
                if not math.isfinite(abundance) or abundance < 0:
                    raise ValueError(f"line {reader.line_num}: invalid abundance")
                if is_species and abundance > 0:
                    counts[i] += 1
    return list(zip(samples, counts))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--index", required=True)
    args = parser.parse_args()
    try:
        counts = count_species(args.input, args.index)
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("# samples", "total features"))
            writer.writerows(counts)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
