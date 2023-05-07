# founderblockgraphs
Constructs segment repeat-free founder block graphs from multiple sequence alignments

# getting started

Clone this repository with dependencies:
```
$ git clone --recurse-submodules https://github.com/algbio/founderblockgraphs.git
$ cd founderblockgraphs
```

Build sdsl-lite:
```
$ cd sdsl-lite
$ export GNUMAKEFLAGS=-j9
$ ./install.sh
$ export GNUMAKEFLAGS=
$ cd ..
```

Build this project (`founderblockgraph`, `locate_multiple`, `locate_patterns`):
```
$ make
```
