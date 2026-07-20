#!/usr/bin/env python3
"""
=============================================================================
Operon Prediction from Metagenome GFF Annotation Files
=============================================================================

Usage
-----
    python metagenome_operons.py [gff_file] [-o output_file]

    With no arguments, defaults to ../data/2088090036.gff.

Example
-------
    python metagenome_operons.py ../data/2088090036.gff

Description
-----------
Predicts operons (longest contiguous multi-gene transcriptional units)
from a GFF3 metagenome annotation file (equivalent to a PTT file).

Rule : Adjacent co-directional genes on the same contig with an intergenic
       distance of less than 50 bp are grouped into the same operon.

       gap = next_gene.start - prev_gene.stop - 1  <  50 bp

GFF3 file format (tab-delimited)
---------------------------------
    Col 1  Contig / sequence ID
    Col 2  Source
    Col 3  Feature type  (only CDS rows are used)
    Col 4  Start position (1-based)
    Col 5  Stop position  (1-based, inclusive)
    Col 6  Score
    Col 7  Strand         "+" or "-"
    Col 8  Frame
    Col 9  Attributes     semicolon-separated key=value pairs

Output
------
Results are written to a text file named:
    <input_filename>_operon_results.txt
in the same directory as the input file.
=============================================================================
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
DEFAULT_INPUT = os.path.join(REPO_ROOT, "data", "2088090036.gff")
RESULTS_DIR = os.path.join(REPO_ROOT, "results")


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

class Gene:
    """Represents a single CDS feature from a GFF file."""

    def __init__(self, contig, start, stop, strand, gene_id, product=""):
        self.contig  = contig
        self.start   = int(start)
        self.stop    = int(stop)
        self.strand  = strand       # "+" or "-"
        self.gene_id = gene_id      # locus tag or ID
        self.product = product      # functional annotation

    def __repr__(self):
        return f"Gene({self.gene_id}, {self.contig}:{self.start}-{self.stop} {self.strand})"


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
        return f"Operon({self.contig}:{self.start}-{self.stop} {self.strand} [{self.size} genes])"


# ---------------------------------------------------------------------------
# Silent warning log
# ---------------------------------------------------------------------------

_warnings = []

def _warn(msg):
    """Log a warning silently — written to output file, not printed to screen."""
    _warnings.append(msg)


# ---------------------------------------------------------------------------
# File opener
# ---------------------------------------------------------------------------

def _open_file(filepath):
    """Open a plain text or .gz compressed file in UTF-8 text mode."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"File not found: '{filepath}'\n"
            f"  Please check the filename and folder path."
        )
    if not os.path.isfile(filepath):
        raise OSError(f"'{filepath}' is not a file.")
    try:
        if filepath.endswith(".gz"):
            return gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
        return open(filepath, "r", encoding="utf-8", errors="replace")
    except OSError as e:
        raise OSError(f"Cannot open '{filepath}': {e}") from e


# ---------------------------------------------------------------------------
# GFF parser
# ---------------------------------------------------------------------------

def parse_gff(filepath):
    """
    Parse a GFF3 annotation file and group CDS features by contig.

    Only rows where column 3 == 'CDS' are processed.
    Gene IDs are extracted from the attributes column (col 9),
    preferring 'locus_tag' over 'ID'.

    Returns
    -------
    dict : contig_name -> list of Gene, each list sorted by start position
    """
    contigs   = defaultdict(list)
    skipped   = 0
    total_cds = 0

    try:
        fh = _open_file(filepath)
    except (FileNotFoundError, OSError) as e:
        _warn(str(e))
        return {}

    try:
        for line_num, line in enumerate(fh, 1):
            line = line.rstrip("\n")

            # Skip comment / directive lines and blank lines
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            # Must have at least 7 columns and be a CDS feature
            if len(parts) < 7 or parts[2] != "CDS":
                continue

            total_cds += 1
            contig = parts[0].strip()
            strand = parts[6].strip()
            attrs  = parts[8].strip() if len(parts) > 8 else ""

            # Validate strand
            if strand not in ("+", "-"):
                _warn(f"Line {line_num}: unexpected strand '{strand}' — skipping.")
                skipped += 1
                continue

            # Parse start/stop coordinates
            try:
                start = int(parts[3])
                stop  = int(parts[4])
            except ValueError:
                _warn(f"Line {line_num}: cannot parse coordinates — skipping.")
                skipped += 1
                continue

            # Reject non-positive coordinates
            if start <= 0 or stop <= 0:
                _warn(f"Line {line_num}: non-positive coordinates ({start}, {stop}) — skipping.")
                skipped += 1
                continue

            # Reject empty contig name
            if not contig:
                _warn(f"Line {line_num}: empty contig name — skipping.")
                skipped += 1
                continue

            # Extract gene identifier from attributes column
            # Priority: locus_tag > ID > fallback to contig_start
            gene_id = ""
            for field in attrs.split(";"):
                field = field.strip()
                if field.startswith("locus_tag="):
                    gene_id = field.split("=", 1)[1].strip()
                    break
                if field.startswith("ID=") and not gene_id:
                    gene_id = field.split("=", 1)[1].strip()
            if not gene_id:
                gene_id = f"{contig}_{start}"

            contigs[contig].append(Gene(
                contig  = contig,
                start   = start,
                stop    = stop,
                strand  = strand,
                gene_id = gene_id,
                product = attrs,
            ))

    except (UnicodeDecodeError, EOFError) as e:
        _warn(f"Read error: {e}. Partial data may be used.")
    finally:
        fh.close()

    if total_cds == 0:
        _warn("No CDS features found. Check that column 3 contains 'CDS'.")
    if skipped:
        _warn(f"Skipped {skipped} malformed line(s).")
    if not contigs:
        _warn("No valid genes were parsed.")

    # Sort genes within each contig by start position
    for contig in contigs:
        contigs[contig].sort(key=lambda g: g.start)

    return contigs


# ---------------------------------------------------------------------------
# Operon prediction
# ---------------------------------------------------------------------------

def predict_operons(genes, max_gap=49):
    """
    Predict operons from a position-sorted list of Gene objects.

    Two adjacent genes are placed in the same operon when:
      1. They are on the same strand  (co-directional).
      2. Their intergenic distance <= max_gap  (default 49 = "< 50 bp").

    Intergenic distance:
        gap = next_gene.start - prev_gene.stop - 1

    Only groups of >= 2 genes are returned as operons.

    Parameters
    ----------
    genes   : list of Gene, sorted by start position on the same contig
    max_gap : int, maximum intergenic distance (default 49)

    Returns
    -------
    list of Operon
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
    """Compute summary statistics for a list of predicted operons."""
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
# Output writer
# ---------------------------------------------------------------------------

def write_results(contig_operons, dataset_name, total_genes, fh):
    """Write full operon predictions and summary stats to an open file handle."""
    all_ops = [op for ops in contig_operons.values() for op in ops]

    fh.write(f"{SEPARATOR}\n")
    fh.write(f"  DATASET : {dataset_name}\n")
    fh.write(f"{SEPARATOR}\n\n")
    fh.write(f"  Contigs containing >= 1 operon : {len(contig_operons)}\n\n")

    if not all_ops:
        fh.write("  No operons predicted for this dataset.\n")
        return

    stats = compute_summary_stats(all_ops, total_genes)

    # --- Global summary statistics ---
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
    fh.write(f"    Location : {lop.contig} : {lop.start} - {lop.stop}  [{lop.strand}]\n")
    fh.write(f"    Genes    : {', '.join(lop.gene_ids)}\n\n")

    sop = stats['longest_span_op']
    fh.write(f"  Longest genomic span operon ({stats['longest_span_bp']} bp):\n")
    fh.write(f"    Location : {sop.contig} : {sop.start} - {sop.stop}  [{sop.strand}]\n")
    fh.write(f"    Genes    : {', '.join(sop.gene_ids)}\n\n")

    fh.write(f"  Operon size distribution:\n")
    fh.write(f"    {'Size':>6}  {'Count':>6}  Bar\n")
    fh.write(f"    {'------':>6}  {'------':>6}  {'-'*40}\n")
    for size, count in stats['size_distribution'].items():
        fh.write(f"    {size:>6}  {count:>6}  {'|' * min(count, 40)}\n")
    fh.write("\n")

    # --- Per-contig detail ---
    fh.write(f"\n  Per-contig detail:\n")
    for contig, ops in sorted(contig_operons.items()):
        fh.write(f"\n  Contig : {contig}   ({len(ops)} operon(s))\n")
        fh.write(f"{COL_HEADER}\n")
        fh.write(f"  {SUB_SEP}\n")
        for idx, op in enumerate(ops, 1):
            fh.write(f"  {idx:<5} {op.start:>10} {op.stop:>10} {op.strand:^8} "
                     f"{op.size:>6}  {', '.join(op.gene_ids)}\n")

    # --- Warnings ---
    if _warnings:
        fh.write(f"\n{SEPARATOR}\n  WARNINGS ({len(_warnings)} total)\n{SEPARATOR}\n")
        for w in _warnings:
            fh.write(f"  ! {w}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Predict operons from a GFF3 metagenome annotation file."
    )
    parser.add_argument(
        "input_file", nargs="?", default=DEFAULT_INPUT,
        help=f"GFF3 file (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output file path (default: results/<input_name>_operon_results.txt)",
    )
    args = parser.parse_args()
    input_file = args.input_file

    # --- Determine output file name ---
    os.makedirs(RESULTS_DIR, exist_ok=True)
    if args.output:
        output_file = args.output
    else:
        base = os.path.basename(input_file)
        out_name = os.path.splitext(base)[0]  # strip .gff extension
        output_file = os.path.join(RESULTS_DIR, out_name + "_operon_results.txt")

    # --- Parse GFF ---
    print(f"Parsing GFF file : {input_file}")
    contigs = parse_gff(input_file)

    if not contigs:
        print("No contigs could be parsed. Please check the input file.")
        sys.exit(1)

    total_genes = sum(len(g) for g in contigs.values())
    print(f"Contigs parsed   : {len(contigs)}")
    print(f"Genes parsed     : {total_genes}")

    # --- Predict operons per contig ---
    contig_operons = {}
    for contig, genes in contigs.items():
        ops = predict_operons(genes, max_gap=49)
        if ops:
            contig_operons[contig] = ops

    all_ops = [op for ops in contig_operons.values() for op in ops]
    in_ops  = sum(op.size for op in all_ops)
    print(f"Contigs with operons : {len(contig_operons)}")
    print(f"Operons predicted    : {len(all_ops)}")
    print(f"Genes in operons     : {in_ops}  ({round(100*in_ops/total_genes, 1)}%)")

    # --- Write output ---
    dataset_name = f"Metagenome: {os.path.basename(input_file)}"
    try:
        with open(output_file, "w", encoding="utf-8") as out:
            out.write("METAGENOME OPERON PREDICTION RESULTS\n")
            out.write("Method : Adjacent co-directional genes with intergenic distance < 50 bp\n")
            out.write("Rule   : gap = next_gene.start - prev_gene.stop - 1 <= 49\n")
            out.write(f"Input  : {input_file}\n")
            out.write(SEPARATOR + "\n")
            write_results(contig_operons, dataset_name, total_genes, out)
    except OSError as e:
        print(f"[ERROR] Cannot write output file: {e}")
        sys.exit(1)

    print(f"Results written  : {output_file}")


if __name__ == "__main__":
    main()