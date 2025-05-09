# exact matching in repeat-free graphs
Branch of `founderblockgraphs` for rewriting program `locate_patterns`, implementing the *expanded search* algorithm for repeat-free founder graphs (and other graphs where nodes are unique). It expects the input graph to be in [GFA format](https://github.com/GFA-spec/GFA-spec) and contain only forward `L` links, but accepts arbitrary suffix-prefix overlap between the edges  (`0M`, `1M`, etc. up to the minimum length of the connected segments minus 1).

## getting started
Download and compile the project
```console
git clone https://github.com/algbio/founderblockgraphs
cd founderblockgraphs
git checkout repeat-free-locate
git submodule update --init --recursive
make
```

```console
$ ./locate_patterns test/small.gfa test/small_patterns.fasta test/small_patterns.gaf --overwrite
Reading the graph... done.
Adding a supersource to the graph... done.
Indexing the graph... done.
Locate
p1: occurs 1 times
p2: occurs 1 times
p3: occurs 7 times (at most)
```

## todo
- optimize memory by using disk
- flag to only perform decision query, no locate
- flag to check whether input graph is repeat-free
- handle reverse complement nodes and links
- multithreading

## cite
The expanded backward search was initially described in

> Massimo Equi, Tuukka Norri, Jarno Alanko, Bastien Cazaux, Alexandru I Tomescu, Veli Mäkinen.
> [*Algorithms and complexity on indexing founder graphs.*](https://doi.org/10.1007/s00453-022-01007-w)
> Algorithmica, 2023.

## ack
This project uses [GFAKluge](https://github.com/edawson/gfakluge) and [kseq.h](https://github.com/lh3/seqtk) for parsing GFA and FASTA, respectively.
