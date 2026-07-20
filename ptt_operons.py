#!/usr/bin/env python3
"""
=============================================================================
Operon Prediction from PTT Genome Annotation Files
=============================================================================

Description
-----------
Predicts operons (longest contiguous multi-gene transcriptional units)
from one or more NCBI PTT (Protein Table) files.

Rule : Adjacent co-directional genes with an intergenic distance of
       less than 50 bp are grouped into the same operon.

       gap = next_gene.start - prev_gene.stop - 1  <  50 bp

HOW TO USE
----------
Run with no arguments to process every *.ptt.gz file in ../data/:

    python ptt_operons.py

Or pass specific file(s) and/or an output path explicitly:

    python ptt_operons.py path/to/genome1.ptt.gz path/to/genome2.ptt.gz -o results.txt

PTT file format (tab-delimited)
--------------------------------
    Line 1  : Genome description
    Line 2  : Number of proteins
    Line 3  : Column headers
    Lines 4+: Location | Strand | Length | PID | Gene | Synonym | Code | COG | Product
=============================================================================
"""

import argparse
import glob
import gzip
import os
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(REPO_ROOT, "data")
DEFAULT_OUTPUT = os.path.join(REPO_ROOT, "results", "ptt_operon_results.txt")


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

class Gene:
    """Represents a single annotated gene from a PTT file."""

    def __init__(self, contig, start, stop, strand, gene_id, synonym="", product=""):
        self.contig  = contig
        self.start   = int(start)
        self.stop    = int(stop)
        self.strand  = strand
        self.gene_id = gene_id
        self.synonym = synonym
        self.product = product

    def __repr__(self):
        return f"Gene({self.gene_id}, {self.start}-{self.stop} {self.strand})"


class Operon:
    """Represents a predicted operon — a group of co-directional adjacent genes."""

    def __init__(self, genes):
        self.genes  = genes
        self.contig = genes[0].contig
        self.strand = genes[0].strand
        self.start  = genes[0].start
        self.stop   = genes[-1].stop

    @property
    def size(self):
        return len(self.genes)

    @property
    def gene_ids(self):
        return [g.gene_id for g in self.genes]

    def __repr__(self):
        return f"Operon({self.start}-{self.stop} {self.strand} [{self.size} genes])"


# ---------------------------------------------------------------------------
# Silent warning log
# ---------------------------------------------------------------------------

_warnings = []

def _warn(msg):
    _warnings.append(msg)


# ---------------------------------------------------------------------------
# File opener
# ---------------------------------------------------------------------------

def _open_file(filepath):
    """Open a plain text or .gz compressed file in UTF-8 text mode."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: '{filepath}'")
    if not os.path.isfile(filepath):
        raise OSError(f"'{filepath}' is not a file.")
    try:
        if filepath.endswith(".gz"):
            return gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
        return open(filepath, "r", encoding="utf-8", errors="replace")
    except OSError as e:
        raise OSError(f"Cannot open '{filepath}': {e}") from e


# ---------------------------------------------------------------------------
# PTT parser
# ---------------------------------------------------------------------------

def parse_ptt(filepath):
    """
    Parse an NCBI PTT file.
    Returns (list of Gene sorted by start position, genome_name string).
    Wrap-around genes on circular chromosomes (start > stop) are included.
    """
    genes       = []
    genome_name = ""
    skipped     = 0

    try:
        fh = _open_file(filepath)
    except (FileNotFoundError, OSError) as e:
        _warn(str(e))
        return [], ""

    try:
        for line_num, line in enumerate(fh):
            line = line.rstrip("\n")

            if line_num == 0:
                genome_name = line.split(" - ")[0].strip() or os.path.basename(filepath)
                continue
            if line_num in (1, 2):
                continue

            parts = line.split("\t")
            if len(parts) < 6:
                skipped += 1
                continue

            location = parts[0]
            strand   = parts[1]
            gene_id  = parts[4]
            synonym  = parts[5]
            product  = parts[8] if len(parts) > 8 else ""

            if strand not in ("+", "-"):
                _warn(f"Line {line_num+1}: unexpected strand '{strand}' — skipping.")
                skipped += 1
                continue

            try:
                loc_parts = location.split("..")
                if len(loc_parts) != 2:
                    raise ValueError
                start, stop = int(loc_parts[0]), int(loc_parts[1])
            except ValueError:
                _warn(f"Line {line_num+1}: cannot parse location '{location}' — skipping.")
                skipped += 1
                continue

            if start <= 0 or stop <= 0:
                _warn(f"Line {line_num+1}: non-positive coordinates ({start},{stop}) — skipping.")
                skipped += 1
                continue

            primary_id = gene_id if gene_id not in ("-", "", None) else synonym
            if not primary_id or primary_id in ("-", ""):
                primary_id = f"gene_{start}"

            genes.append(Gene(
                contig  = genome_name,
                start   = start,
                stop    = stop,
                strand  = strand,
                gene_id = primary_id,
                synonym = synonym,
                product = product,
            ))

    except (UnicodeDecodeError, EOFError) as e:
        _warn(f"Read error in '{filepath}': {e}. Partial data used.")
    finally:
        fh.close()

    if skipped:
        _warn(f"Skipped {skipped} malformed line(s) in '{os.path.basename(filepath)}'.")
    if not genes:
        _warn(f"No valid genes found in '{os.path.basename(filepath)}'.")

    genes.sort(key=lambda g: g.start)
    return genes, genome_name


# ---------------------------------------------------------------------------
# Operon prediction
# ---------------------------------------------------------------------------

def predict_operons(genes, max_gap=49):
    """
    Predict operons from a sorted list of Gene objects.

    Two adjacent genes join the same operon when:
      1. Same strand (co-directional)
      2. Intergenic gap <= max_gap  (gap = next.start - prev.stop - 1)

    Only groups of >= 2 genes are returned as operons.
    """
    if not genes:
        return []

    operons       = []
    current_group = [genes[0]]

    for i in range(1, len(genes)):
        prev = current_group[-1]
        curr = genes[i]
        gap  = curr.start - prev.stop - 1

        if curr.strand == prev.strand and gap <= max_gap:
            current_group.append(curr)
        else:
            if len(current_group) >= 2:
                operons.append(Operon(current_group))
            current_group = [curr]

    if len(current_group) >= 2:
        operons.append(Operon(current_group))

    return operons


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

SEPARATOR  = "=" * 80
SUB_SEP    = "-" * 80
COL_HEADER = f"  {'#':<5} {'Start':>10} {'Stop':>10} {'Strand':^8} {'NGenes':>6}  Gene IDs"


def _size_distribution(operons):
    dist = defaultdict(int)
    for op in operons:
        dist[op.size] += 1
    return dict(sorted(dist.items()))


def compute_summary_stats(operons, total_genes):
    if not operons:
        return {}

    sizes        = [op.size for op in operons]
    genes_in_ops = sum(sizes)
    singletons   = max(0, total_genes - genes_in_ops)
    n            = len(sizes)
    sorted_sizes = sorted(sizes)
    mean_size    = sum(sorted_sizes) / n
    median_size  = (sorted_sizes[n // 2] if n % 2 == 1
                    else (sorted_sizes[n // 2 - 1] + sorted_sizes[n // 2]) / 2.0)

    return {
        "total_operons"        : len(operons),
        "total_genes"          : total_genes,
        "genes_in_operons"     : genes_in_ops,
        "singleton_genes"      : singletons,
        "pct_genes_in_operons" : round(100.0 * genes_in_ops / total_genes, 2) if total_genes else 0,
        "min_operon_size"      : min(sizes),
        "max_operon_size"      : max(sizes),
        "mean_operon_size"     : round(mean_size, 2),
        "median_operon_size"   : median_size,
        "fwd_strand_operons"   : sum(1 for op in operons if op.strand == "+"),
        "rev_strand_operons"   : sum(1 for op in operons if op.strand == "-"),
        "largest_operon"       : max(operons, key=lambda op: op.size),
        "longest_span_op"      : max(operons, key=lambda op: op.stop - op.start),
        "longest_span_bp"      : max(op.stop - op.start + 1 for op in operons),
        "size_distribution"    : _size_distribution(operons),
    }


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_genome_results(operons, genome_name, total_genes, fh):
    """Write one genome's operon predictions and summary stats to file."""
    fh.write(f"\n{SEPARATOR}\n")
    fh.write(f"  GENOME : {genome_name}\n")
    fh.write(f"{SEPARATOR}\n\n")

    if not operons:
        fh.write("  No operons predicted for this genome.\n\n")
        return

    stats = compute_summary_stats(operons, total_genes)

    # Summary statistics block
    fh.write("  --- Summary Statistics ---\n")
    fh.write(f"  Total genes (annotated)       : {stats['total_genes']}\n")
    fh.write(f"  Total operons predicted       : {stats['total_operons']}\n")
    fh.write(f"  Genes in operons              : {stats['genes_in_operons']}"
             f"  ({stats['pct_genes_in_operons']}% of all genes)\n")
    fh.write(f"  Singleton genes (no operon)   : {stats['singleton_genes']}"
             f"  ({round(100 - stats['pct_genes_in_operons'], 2)}% of all genes)\n\n")
    fh.write(f"  Operon size (number of genes):\n")
    fh.write(f"    Minimum   : {stats['min_operon_size']}\n")
    fh.write(f"    Maximum   : {stats['max_operon_size']}\n")
    fh.write(f"    Mean      : {stats['mean_operon_size']}\n")
    fh.write(f"    Median    : {stats['median_operon_size']}\n\n")
    fh.write(f"  Strand breakdown:\n")
    fh.write(f"    Forward (+) operons : {stats['fwd_strand_operons']}\n")
    fh.write(f"    Reverse (-) operons : {stats['rev_strand_operons']}\n\n")

    lop = stats['largest_operon']
    fh.write(f"  Largest operon ({lop.size} genes):\n")
    fh.write(f"    Location : {lop.start} - {lop.stop}  [{lop.strand}]\n")
    fh.write(f"    Genes    : {', '.join(lop.gene_ids)}\n\n")

    sop = stats['longest_span_op']
    fh.write(f"  Longest genomic span operon ({stats['longest_span_bp']} bp):\n")
    fh.write(f"    Location : {sop.start} - {sop.stop}  [{sop.strand}]\n")
    fh.write(f"    Genes    : {', '.join(sop.gene_ids)}\n\n")

    fh.write(f"  Operon size distribution:\n")
    fh.write(f"    {'Size':>6}  {'Count':>6}  Bar\n")
    fh.write(f"    {'------':>6}  {'------':>6}  {'-'*40}\n")
    for size, count in stats['size_distribution'].items():
        fh.write(f"    {size:>6}  {count:>6}  {'|' * min(count, 40)}\n")
    fh.write("\n")

    # Detailed operon list
    fh.write(f"  Detailed operon list:\n")
    fh.write(f"{COL_HEADER}\n")
    fh.write(f"  {SUB_SEP}\n")
    for idx, op in enumerate(operons, 1):
        fh.write(f"  {idx:<5} {op.start:>10} {op.stop:>10} {op.strand:^8} "
                 f"{op.size:>6}  {', '.join(op.gene_ids)}\n")


def write_combined_summary(summary_rows, fh=None):
    """
    Print combined summary table to console and optionally to output file.
    Handles both successfully processed genomes and skipped files.
    """
    lines = [
        "",
        SEPARATOR,
        f"  COMBINED SUMMARY  —  {len(summary_rows)} genome(s) processed",
        SEPARATOR,
        f"  {'Genome':<45} {'Genes':>7}  {'Operons':>8}  {'In operons':>10}  {'Coverage':>9}",
        f"  {'-'*45} {'-'*7}  {'-'*8}  {'-'*10}  {'-'*9}",
    ]

    for genome_name, n_genes, operons in summary_rows:
        in_op    = sum(op.size for op in operons)
        if n_genes > 0:
            coverage = round(100.0 * in_op / n_genes, 1)
            status   = f"{coverage:>8}%"
        else:
            status   = "  SKIPPED"
        lines.append(
            f"  {genome_name[:45]:<45} {n_genes:>7}  {len(operons):>8}  {in_op:>10}  {status}"
        )

    lines += [SEPARATOR, ""]

    for line in lines:
        print(line)
        if fh:
            fh.write(line + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    default_inputs = sorted(glob.glob(os.path.join(DATA_DIR, "*.ptt.gz")))

    parser = argparse.ArgumentParser(
        description="Predict operons from one or more NCBI PTT genome annotation files."
    )
    parser.add_argument(
        "input_files", nargs="*", default=default_inputs,
        help=f"PTT file(s) (.ptt or .ptt.gz). Defaults to every *.ptt.gz in {DATA_DIR}/",
    )
    parser.add_argument(
        "-o", "--output", default=DEFAULT_OUTPUT,
        help=f"Output file path (default: {DEFAULT_OUTPUT})",
    )
    args = parser.parse_args()
    input_files = args.input_files
    output_file = args.output

    # --- Validate ---
    if not input_files:
        print("[ERROR] No input files specified. Pass PTT path(s) as arguments, "
              f"or place *.ptt.gz files in {DATA_DIR}/.")
        sys.exit(1)

    # --- Create output directory if it doesn't exist ---
    out_dir = os.path.dirname(output_file)
    if out_dir:
        try:
            os.makedirs(out_dir, exist_ok=True)
        except OSError as e:
            print(f"[ERROR] Cannot create output directory '{out_dir}': {e}")
            sys.exit(1)

    # --- Open output file ---
    try:
        out = open(output_file, "w", encoding="utf-8")
    except OSError as e:
        print(f"[ERROR] Cannot create output file '{output_file}': {e}")
        sys.exit(1)

    n_total      = len(input_files)
    summary_rows = []

    with out:
        # File header
        out.write("OPERON PREDICTION RESULTS  —  PTT Genome(s)\n")
        out.write("Method : Adjacent co-directional genes with intergenic distance < 50 bp\n")
        out.write("Rule   : gap = next_gene.start - prev_gene.stop - 1 <= 49\n")
        out.write(f"Files  : {n_total} PTT file(s)\n")
        out.write(SEPARATOR + "\n")

        # Process each file
        for idx, input_file in enumerate(input_files, 1):
            fname = os.path.basename(input_file)
            print(f"[{idx}/{n_total}] Parsing : {fname} ...", flush=True)

            _warnings.clear()

            genes, genome_name = parse_ptt(input_file)

            if not genes:
                msg = f"SKIPPED — no genes could be parsed from '{fname}'"
                print(f"        -> {msg}")
                out.write(f"\n{SEPARATOR}\n  GENOME : {fname}\n{SEPARATOR}\n  {msg}\n")
                summary_rows.append((fname, 0, []))
                continue

            operons  = predict_operons(genes, max_gap=49)
            in_ops   = sum(op.size for op in operons)
            coverage = round(100.0 * in_ops / len(genes), 1)

            print(f"        -> {len(genes):,} genes | {len(operons):,} operons | {coverage}% in operons")

            write_genome_results(operons, genome_name, len(genes), out)

            if _warnings:
                out.write(f"  Warnings for this file:\n")
                for w in _warnings:
                    out.write(f"    ! {w}\n")
                out.write("\n")

            summary_rows.append((genome_name, len(genes), operons))

        # Combined summary at end of file
        write_combined_summary(summary_rows, fh=out)

    print(f"\nResults written -> {output_file}")


if __name__ == "__main__":
    main()