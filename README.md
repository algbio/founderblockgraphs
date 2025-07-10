# rewrite branch
Rewrite of tool `founderblockgraph` focused on memory usage in favor of disk usage. Tested on GCC 15 and based on [SDSL v3](https://github.com/xxsds/sdsl-lite) for compact data structures and [kseq.h](https://github.com/lh3/seqtk)(kseq.h) for FASTA parsing.

# getting started
```console
git clone https://github.com/algbio/founderblockgraphs.git
cd founderblockgraphs
git checkout rewrite
git submodule update --init --recursive
make
```

# usage
See `founderblockgraph -h`.

# citation
The default construction algorithms are from
> Nicola Rizzo, Massimo Equi, Tuukka Norri, Veli Mäkinen.
> [*Elastic founder graphs improved and enhanced*](https://doi.org/10.1016/j.tcs.2023.114269).
> Theoretical Computer Science, 2024.

with many optimizations (better handling of long sequences, parallelization, `--heuristic-subset`, `--ignore-chars`) introduced in
> Nicola Rizzo, Manuel Cáceres, Veli Mäkinen.
> [*Exploiting uniqueness: seed-chain-extend alignment on elastic founder graphs*](https://doi.org/10.1093/bioinformatics/btaf225).
> ISMB 2025.

The non-elastic construction algorithms (`--non-elastic`) are from
> Massimo Equi, Tuukka Norri, Jarno Alanko, Bastien Cazaux, Alexandru I. Tomescu, Veli Mäkinen.
> [*Algorithms and Complexity on Indexing Founder Graphs*](https://doi.org/10.1007/s00453-022-01007-w).
> Algorithmica, 2023.

# todo
 - examples/tests
 - refactor non-elastic algorithms into this branch
 - document EFG tricks related to option `--ignore-chars`, to the start and end of sequences, and to initial and ending runs of gaps
 - allow empty segments when long runs of gaps are present