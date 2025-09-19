<div align="center">
<img src="/assets/icon.svg" width="300"></img>
</div><br><br>

# Myrio

Myrio is a command-line application designed to identify the taxonomy of plants (and potentially other organisms) using amplified sequences (potentially mixed) of barcode genes.

The name Myrio is inspired by the scientific name of the plant _Myriophyllum Spicatum_, commonly known as Eurasian Watermilfoil, an aquatic plant found in the [Léman](https://en.wikipedia.org/wiki/Lake_Geneva).

## Features
* Cross-platform (windows, macOS, and linux are all supported)
* Built with Rust (free from the hassle of installing/using Python code)
* Zero external dependencies, `myrio` won't crash if you haven't installed another binary or library
* Optimized codebase, including but not limited to:
    * Custom sparse vector implementation with efficient operations
    * Parallelism via [`Rayon`](https://github.com/rayon-rs/rayon)
    * Specialized database format able to store pre-computed k-mer counts efficiently
* Flexible output, results can be exported as `.csv` or visualized as a tree in `.txt` format
* Computation over heuristics, relies more on raw parallel computation and memory-efficient design rather than heuristics. For example, Myrio uses full k-mer counts (not just sets, and no minimizers).

## Table of Contents

- [Installation](#Installation)
- [Usage](#Usage)
- [Acknowledgments](#Acknowledgments)

## Installation

To build and install `myrio` on your system, you'll need [rust](https://www.rust-lang.org/learn/get-started) installed. Then run:
``` sh
git clone https://github.com/anesthetice/Myrio
cd Myrio
RUSTFLAGS="-C target-cpu=native" cargo install --path myrio-cli
```

## Usage

### Creating a database (a "tree")

Each reference database corresponds to a single barcode gene (e.g. one database for _matK_, another for _rbcL_, etc.). Databases are generated from a single `.fasta` file.

FASTA entries must contain a `tax={...}` annotation. For example:
``` FASTA
>BOLD_PROCESS_ID=ZPLPP049-13|tax={p:Tracheophyta, c:Magnoliopsida, o:Rosales, f:Rosaceae, g:Prunus, s:Prunus persica}
ATACCCTACCCCATTCATCTGGAAATCTTGGTTCAAACCCTTCGCTATTGGGTGAAAGACGCCTCTTCTTTGCATTTATTACGACTCTTTCTTCACGAGTATTATAATTGGAATAG...
```

The parser is flexible, so the following would also pass:
``` FASTA
>BOLD_PROCESS_ID=ZPLPP049-13|tax={domain: Eukarya, kingdom: Plantae, phylum: Tracheophyta, class: Magnoliopsida, order: Rosales, family: Rosaceae, genus: Prunus, species: Prunus persica}
ATACCCTACCCCATTCATCTGGAAATCTTGGTTCAAACCCTTCGCTATTGGGTGAAAGACGCCTCTTCTTTGCATTTATTACGACTCTTTCTTCACGAGTATTATAATTGGAATAG...

>BOLD_PROCESS_ID=ZPLPP049-13|tax={g:prunus; species: Prunus persica;}
ATACCCTACCCCATTCATCTGGAAATCTTGGTTCAAACCCTTCGCTATTGGGTGAAAGACGCCTCTTCTTTGCATTTATTACGACTCTTTCTTCACGAGTATTATAATTGGAATAG...
```

**Important constraints**:

1. All entries must share the same highest-ranked clade. For example, if the highest rank of the first record is `family: Araliaceae`, then every other record must also have `family: Araliaceae` as their highest-ranked clade (note that the highest rank defined is `Domain`, while the lowest is `Species`).
2. No rank gaps are allowed. For instance, if you specify `family`, you cannot skip `genus` and go directly to `species`.


Once your FASTA database precursor file is ready, you can create a database with:
``` sh
myrio tree new --input BOLD_Plantae_20250831_ITS.fasta --gene "ITS"
# And if we want to pre-compute k-mer counts with `k=18` (highly recommended, will significantly increase database size however)
myrio tree new --input BOLD_Plantae_20250831_ITS.fasta --gene "ITS" -k 18
```

This will create a file called `BOLD_Plantae_20250831_ITS.myrtree`.

If errors are encountered, they will be reported and the problematic entry skipped. For example:

```
Failed to parse taxonomic identity of the record starting on line 356021, Failed to parse string into a list of clade: cannot have rank gaps, expected 6 elements, got 5; string: '>BOLD_PROCESS_ID=MHPAF950-11|tax={p:Tracheophyta, c:Liliopsida, o:Poales, f:Poaceae, s:Poaceae A.guadamuz275}'

Bio-seq parsing error for the record starting on line 393997, Unrecognised character: 'I' (0x49)
```

See [myrio-py/db_gen.py](/myrio-py/db_gen.py) for an example of how to generate a FASTA precursor database (it's a [marimo](https://marimo.io/) notebook)

Pre-built databases are also available on the releases page.

### Running
Example runs:
``` sh
# Input must be a single `.fastq` file.
# `--trees` can be a directory containing multiple `.myrtree` databases.
myrio run --input Berberis_Julianae_matK_rbcL_psbA-trnH_ITS.fastq --trees myrio-db/

# Also, k-mer counts computed (if not already pre-computed) can be cached directly into their respective database file.
myrio run --input Berberis_Julianae_matK_rbcL_psbA-trnH_ITS.fastq --trees myrio-db/ --k-search 19 --cache-counts
# Note that if `--k-search` is not provided, the value is read from the configuration file (`~/.config/myrio/myrio.conf.toml`).

# If you expect more clusters than gene databases, you can set `--nb-clusters`.
myrio run --input Berberis_Julianae_matK_rbcL_psbA-trnH_ITS.fastq --trees myrio-db/BOLD_Plantae_20250831_ITS.myrtree --nb-clusters 4
```

## Acknowledgments

* Special thanks to the [Paoli Lab](https://www.epfl.ch/labs/gr-paoli/) for hosting this project.
* Special thanks to [GenoRobotics](https://www.genorobotics.org/), and especially our team for the 2025 Lemanic Life Sciences Hackathon, which built the [proof-of-concept](https://github.com/GenoRobotics-EPFL/Myrio-Hackathon) for this application.
