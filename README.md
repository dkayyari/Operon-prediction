# Operon Prediction from Genome and Metagenome Annotations

![Python](https://img.shields.io/badge/python-3.x-blue)
![License](https://img.shields.io/badge/license-MIT-green)

Predicts **operons** — contiguous groups of co-directional genes likely transcribed
together — from prokaryotic genome and metagenome annotation files, using a simple
intergenic-distance rule.

## Overview

In bacteria and archaea, genes that are close together on the same strand are often
transcribed as a single mRNA (an operon), rather than individually. This project
predicts operon structure directly from genome coordinates:

> Two adjacent genes on the same strand belong to the same operon if the intergenic
> gap between them is **less than 50 bp**
> (`gap = next_gene.start - prev_gene.stop - 1`).

Two scripts apply this rule to two different kinds of input:

- **`ptt_operons.py`** — processes one or more NCBI PTT (Protein Table) files from
  fully assembled genomes, and can batch multiple genomes into a single combined
  report with a cross-genome summary.
- **`metagenome_operons.py`** — processes a GFF3 metagenome annotation, where genes
  are spread across thousands of contigs rather than one chromosome, and reports
  operons per contig.

Both scripts compute the same summary statistics (operon count, size distribution,
strand breakdown, largest/longest operons, genome coverage) and write a structured
text report.

## Results

**Four reference genomes (`ptt_operons.py`):**

| Genome | Genes | Operons | Genes in operons | Coverage |
|---|---:|---:|---:|---:|
| *E. coli* K-12 MG1655 | 4,146 | 784 | 2,270 | 54.8% |
| *B. subtilis* 168 | 4,176 | 771 | 2,294 | 54.9% |
| *Halobacterium* NRC-1 | 2,075 | 396 | 1,017 | 49.0% |
| *Synechocystis* PCC 6803 | 3,170 | 488 | 1,151 | 36.3% |

Roughly half of each genome's genes fall into a predicted operon, consistent with
the known prevalence of operons in prokaryotic genomes — with the cyanobacterium
*Synechocystis* showing noticeably lower operon coverage than the other three.

**Metagenome (`metagenome_operons.py`):**

Across 23,567 genes on 21,401 contigs, 759 operons were predicted spanning 722
contigs, covering 1,657 genes (7.0%). The much lower coverage than the complete
genomes above is expected — metagenome assemblies are fragmented into many short
contigs, so most genes lack a neighboring gene on the same contig to form an operon
with.

Full per-genome and per-contig detail, including size distributions and the
largest/longest operon found in each dataset, is in `results/`.

## Project structure

```
.
├── data/
│   ├── E_coli_K12_MG1655.ptt.gz
│   ├── B_subtilis_168.ptt.gz
│   ├── Halobacterium_NRC1.ptt.gz
│   ├── Synechocystis_PCC6803_uid159873.ptt.gz
│   └── 2088090036.gff                     # Metagenome annotation
├── src/
│   ├── ptt_operons.py                     # Genome PTT files -> operons
│   └── metagenome_operons.py              # Metagenome GFF3 -> operons
├── results/
│   ├── ptt_operon_results.txt
│   └── Metagenome_operon_results.txt
└── LICENSE
```

## Getting started

No external dependencies — just the Python standard library (`gzip`, `os`, `sys`,
`argparse`, `collections`).

```bash
git clone https://github.com/<your-username>/operon-prediction.git
cd operon-prediction
```

**Genome PTT files** — with no arguments, processes every `*.ptt.gz` in `data/`:

```bash
python src/ptt_operons.py
# or specify files explicitly:
python src/ptt_operons.py data/E_coli_K12_MG1655.ptt.gz -o results/ecoli_operons.txt
```

**Metagenome GFF3** — defaults to `data/2088090036.gff`:

```bash
python src/metagenome_operons.py
# or specify a file explicitly:
python src/metagenome_operons.py data/2088090036.gff -o results/my_results.txt
```

## Method notes

- Genes are sorted by start position before prediction (per contig, for the
  metagenome script).
- Wrap-around genes on circular chromosomes are included as-is.
- Malformed or unparseable annotation lines are logged as warnings in the output
  file rather than halting the run.
- Only groups of **2 or more** genes are reported as operons; singleton genes are
  excluded.
- A PTT gene's ID is preferred from the gene-name column, falling back to its
  synonym; a GFF gene's ID prefers `locus_tag` over `ID` in the attributes column.

## Tech stack

- Python (standard library only)

## License

This project is licensed under the [MIT License](LICENSE).
