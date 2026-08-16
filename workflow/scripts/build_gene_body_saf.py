#!/usr/bin/env python3
"""Create a featureCounts SAF file containing one interval per annotated gene.

This intentionally treats the complete genomic span of each gene as a feature,
so intronic as well as exonic reads can contribute to nuclear-RNA abundance.
The conventional exon-count matrix is generated separately by the workflow.
"""

import argparse
import gzip
import re
from pathlib import Path


def op(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def parse_attrs(s: str) -> dict[str, str]:
    out = {}
    for field in s.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        m = re.match(r'([^ ]+)\s+"([^"]*)"', field)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--saf", required=True)
    ap.add_argument("--map", required=True, help="Gene-ID to symbol/biotype map")
    args = ap.parse_args()

    Path(args.saf).parent.mkdir(parents=True, exist_ok=True)
    genes = []
    with op(args.gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "gene":
                continue
            attrs = parse_attrs(f[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            genes.append(
                (
                    gene_id,
                    f[0],
                    int(f[3]),
                    int(f[4]),
                    f[6],
                    attrs.get("gene_name", ""),
                    attrs.get("gene_type", attrs.get("gene_biotype", "")),
                )
            )

    if not genes:
        raise SystemExit("No gene records found in the GTF; cannot build SAF.")

    with open(args.saf, "w") as out:
        out.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for gid, chrom, start, end, strand, _, _ in genes:
            out.write(f"{gid}\t{chrom}\t{start}\t{end}\t{strand}\n")

    with open(args.map, "w") as out:
        out.write("gene_id\tgene_name\tgene_type\n")
        for gid, _, _, _, _, name, biotype in genes:
            out.write(f"{gid}\t{name}\t{biotype}\n")

    print(f"Wrote {len(genes):,} gene bodies to {args.saf}")


if __name__ == "__main__":
    main()
