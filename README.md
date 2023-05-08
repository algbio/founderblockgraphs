# founderblockgraphs
Constructs repeat-free/semi-repeat-free non-elastic/elastic founder graphs from multiple sequence alignments.

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

# usage
```
$ ./founderblockgraphs -h
Usage: founderblockgraph --input=MSA.fasta --output={MSA.index|efg.xgfa} [--elastic]
[--gap-limit=GAPLIMIT] [--threads=THREADNUM] [--graphviz-output=…]
[--output-paths]
Constructs a semi-repeat-free (Elastic) Founder Graph

Input is MSA given in fasta format. Rows with runs of gaps ‘-’ or N’s ≥
GAPLIMIT will be filtered out

  -h, --help                    Print help and exit
  -V, --version                 Print version and exit
      --input=filename          MSA input path
      --output=filename         Index/EFG output path
      --gap-limit=GAPLIMIT      Gap limit (incompatible with --elastic)
                                  (default=`1')
  -e, --elastic                 Min-max-length semi-repeat-free segmentation
                                  (default=off)
      --gfa                     Saves output in xGFA format  (default=off)
  -p, --output-paths            Print the original sequences as paths of the
                                  xGFA graph (requires --gfa).  (default=off)
      --graphviz-output=filename
                                Graphviz output path
      --memory-chart-output=filename
                                Memory chart output path
```

# EFG tricks
The semi-repeat-free property is completely ineffective when one of the sequences is a proper suffix of some other sequence, or when this is the case for a large prefix of the alignment.
This limitation is easily overcome with the following trick: we assume that each sequence starts with a unique `$` character at column 0 (but this character is not actually kept in the graph): pattern matching is **not** affected by this change in the first block.

# todo
 - implement pattern matching on semi-repeat-free EFGs
