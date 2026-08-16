#!/usr/bin/env python3
"""Convert featureCounts output to a clean genes x samples integer matrix."""

import argparse
import csv
from pathlib import Path


def read_samples(path):
    with open(path, newline="") as fh:
        return [row["sample"] for row in csv.DictReader(fh, delimiter="\t")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    samples = read_samples(args.samples)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    with open(args.input) as src:
        lines = (line for line in src if not line.startswith("#"))
        reader = csv.reader(lines, delimiter="\t")
        header = next(reader)
        if len(header) < 7:
            raise SystemExit("Unexpected featureCounts output: fewer than 7 columns")
        count_columns = header[6:]
        if len(count_columns) != len(samples):
            raise SystemExit(
                f"featureCounts has {len(count_columns)} sample columns but samples.tsv has {len(samples)}"
            )

        with open(args.output, "w", newline="") as dst:
            writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id", *samples])
            for row in reader:
                writer.writerow([row[0], *row[6:]])


if __name__ == "__main__":
    main()
