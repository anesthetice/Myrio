<div align="center">
<img src="/assets/icon.svg" width="300"></img>
</div><br><br>

# Myrio

Myrio is a command-line application designed to identify the taxonomy of a plant from specific amplified sequences of a sample.

The name Myrio is inspired by the scientific name Myriophyllum spicatum, commonly known as Eurasian watermilfoil, an aquatic plant found in the Léman.

## Installation

TODO

## Usage

### Creating a database (also called 'tree')

Reference databases are distinct by the gene they represent (e.g., the 'matK' gene is represented by a single database, the 'rbcL' gene by a single other, etc.);
they are generated using a single FASTA file as input. Within the FASTA file, entries must contain a `tax={...}` section such as:
``` FASTA
>BOLD_PROCESS_ID=ZPLPP049-13|tax={p:Tracheophyta, c:Magnoliopsida, o:Rosales, f:Rosaceae, g:Prunus, s:Prunus persica}
ATACCCTACCCCATTCATCTGGAAATCTTGGTTCAAACCCTTCGCTATTGGGTGAAAGACGCCTCTTCTTTGCATTTATTACGACTCTTTCTTCACGAGTATTATAATTGGAATAG...
```

The following are also valid, syntax is not the strictest:
``` FASTA
>BOLD_PROCESS_ID=ZPLPP049-13|tax={domain: Eukarya, kingdom: Plantae, phylum: Tracheophyta, class: Magnoliopsida, order: Rosales, family: Rosaceae, genus: Prunus, species: Prunus persica}
ATACCCTACCCCATTCATCTGGAAATCTTGGTTCAAACCCTTCGCTATTGGGTGAAAGACGCCTCTTCTTTGCATTTATTACGACTCTTTCTTCACGAGTATTATAATTGGAATAG...

>BOLD_PROCESS_ID=ZPLPP049-13|tax={g:prunus; species: Prunus persica;}
ATACCCTACCCCATTCATCTGGAAATCTTGGTTCAAACCCTTCGCTATTGGGTGAAAGACGCCTCTTCTTTGCATTTATTACGACTCTTTCTTCACGAGTATTATAATTGGAATAG...
```

Note two important things for this to work:
1. Every entry must have the same highest rank, e.g., you can't have an entry that only specifies a genus and a species while a previous entry specified up to the 'family' rank.
2. Entries cannot contain rank gaps, e.g. if you specify the family, then you cannot skip over the genus rank and only specify the species.


Finally, we can create the actual tree used by 'Myrio' by running the `myrio tree new` subcommand, which should look something like this:
``` sh
myrio tree new --input BOLD_Plantae_20250831_ITS.fasta --gene "ITS"
# And if we want to pre-compute k-mer counts with `k=18` (highly recommended, will significantly increase database size however)
myrio tree new --input BOLD_Plantae_20250831_ITS.fasta --gene "ITS" -k 18
```

This will create a file called `BOLD_Plantae_20250831_ITS.myrtree`, note that any errors encountered will be printed and their associated entry skipped, for example:

```
Failed to parse taxonomic identity of the record starting on line 356021, Failed to parse string into a list of clade: cannot have rank gaps, expected 6 elements, got 5; string: '>BOLD_PROCESS_ID=MHPAF950-11|tax={p:Tracheophyta, c:Liliopsida, o:Poales, f:Poaceae, s:Poaceae A.guadamuz275}'
```

```
Bio-seq parsing error for the record starting on line 393997, Unrecognised character: 'I' (0x49)
```

You can find an example of how to create your own FASTA precursor database at [./myrio-py/db_gen.py](/myrio-py/db_gen.py) (it's a [marimo](https://marimo.io/) notebook)

### Running

``` sh
# Input must be a `.fastq` file, databases can be specified as a directory, meaning any file ending with `.myrtree` within said directory will be used
myrio run --input Berberis_Julianae_matK_rbcL_psbA-trnH_ITS.fastq --trees myrio-db/

# The k-mer counts computed (if not already pre-computed) can also be cached into their respective database file
myrio run --input Berberis_Julianae_matK_rbcL_psbA-trnH_ITS.fastq --trees myrio-db/ --k-search 19 --cache-counts # Note that if `--k-search` is not specified, then `k-search` will be defined by the application's config, which can be found at `~/.config/myrio/myrio.conf.toml`

# Finally, if you expect to have more clusters than gene databases, you can do the following
myrio run --input Berberis_Julianae_matK_rbcL_psbA-trnH_ITS.fastq --trees myrio-db/BOLD_Plantae_20250831_ITS.myrtree --nb-clusters 4
```









## References
