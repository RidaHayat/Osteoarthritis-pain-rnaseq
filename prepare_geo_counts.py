#!/usr/bin/env python3
"""Extract the 10 mRNA samples from GEO's GSE99662 gene-count matrix.

Download GSE99662_gene_counts_matrix.txt.gz from the GEO Series supplementary
files and place it in data/ before running this script.
"""
import csv
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "data" / "GSE99662_gene_counts_matrix.txt.gz"
COUNTS_OUT = ROOT / "data" / "counts_geo.tsv"
METADATA_OUT = ROOT / "data" / "metadata.csv"

# GEO matrix columns -> stable sample IDs used by the study metadata -> phenotype.
SAMPLES = [
    ("P_12_H", "GENS00000733", "Pain:High"),
    ("P_13_L", "GENS00000734", "Pain:Low"),
    ("P_14_L", "GENS00000735", "Pain:Low"),
    ("P_15_L", "GENS00000736", "Pain:Low"),
    ("P_16_H", "GENS00000737", "Pain:High"),
    ("P_17_H", "GENS00000738", "Pain:High"),
    ("P_18_H", "GENS00000739", "Pain:High"),
    ("P_19_H", "GENS00000740", "Pain:High"),
    ("P_21_L", "GENS00000741", "Pain:Low"),
    ("P_22_L", "GENS00000742", "Pain:Low"),
]

if not SOURCE.exists():
    raise SystemExit(f"Missing {SOURCE}. Download the GEO supplementary matrix first.")

with gzip.open(SOURCE, "rt", newline="") as source:
    reader = csv.DictReader(source, delimiter="\t")
    required = {"gene_id", "gene_name", *(column for column, _, _ in SAMPLES)}
    missing = required.difference(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"GEO matrix is missing columns: {', '.join(sorted(missing))}")
    with COUNTS_OUT.open("w", newline="") as target:
        writer = csv.writer(target, delimiter="\t")
        writer.writerow(["ID", *(sample_id for _, sample_id, _ in SAMPLES)])
        for row in reader:
            gene_id = row["gene_id"]
            gene_name = row["gene_name"]
            identifier = f"{gene_id}_{gene_name}" if gene_name and gene_name != "-" else gene_id
            writer.writerow([identifier, *(row[column] for column, _, _ in SAMPLES)])

with METADATA_OUT.open("w", newline="") as target:
    writer = csv.writer(target)
    writer.writerow(["Sample.ID", "Phenotype"])
    writer.writerows((sample_id, phenotype) for _, sample_id, phenotype in SAMPLES)

print(f"Wrote {COUNTS_OUT}")
print(f"Wrote {METADATA_OUT} (sample IDs and phenotype only)")
