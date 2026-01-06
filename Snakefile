# Snakefile for ProteoBench benchmark download and conversion

import glob
from pathlib import Path

BENCHMARK = "dia_aif"
OUTPUT_DIR = "benchmark_data"


def get_input_files():
    """Find all downloaded input_file.txt files."""
    return glob.glob(f"{OUTPUT_DIR}/{BENCHMARK}/*/*/input_file.txt")


def get_h5ad_targets(wildcards=None):
    """Get expected h5ad output paths from existing input files."""
    return [f.replace("input_file.txt", "output.h5ad") for f in get_input_files()]


rule all:
    input:
        get_h5ad_targets


rule download:
    output:
        directory(f"{OUTPUT_DIR}")
    shell:
        "download_data --benchmark {BENCHMARK} --output {OUTPUT_DIR}"


rule convert:
    input:
        "{path}/input_file.txt"
    output:
        "{path}/output.h5ad"
    shell:
        "prot2ad convert {input} {output}"


rule convert_all:
    input:
        get_h5ad_targets


rule help:
    """Show available rules."""
    run:
        print("""
Snakefile for ProteoBench benchmark download and conversion

Rules:
  snakemake help -c1         Show this help message
  snakemake download -c1     Download benchmark data (dia_aif)
  snakemake convert_all -c1  Convert all downloaded files to h5ad
  snakemake -c1              Default: convert all (requires download first)

Configuration:
  BENCHMARK = "dia_aif"
  OUTPUT_DIR = "benchmark_data"
""")
