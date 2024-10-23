# founderblockgraphs
Constructs repeat-free/semi-repeat-free non-elastic/elastic founder graphs from multiple sequence alignments.

# getting started
Clone this repository with dependencies:
```
$ git clone --recurse-submodules https://github.com/algbio/founderblockgraphs.git
$ cd founderblockgraphs
```

Build this project (`founderblockgraph`, `locate_multiple`, `locate_patterns`):
```
$ make
```

# usage
```
Usage: founderblockgraph --input=MSA.fasta --output={MSA.index|efg.xgfa} [--gfa]
[--elastic] [--gap-limit=GAPLIMIT] [--threads=THREADNUM]
[--graphviz-output=efg.dot] [--output-paths] [--ignore-chars="ALPHABET"]
Constructs a semi-repeat-free (Elastic) Founder Graph

Input is MSA given in fasta format. In standard mode (without --elastic), rows
with runs of gaps ‘-’ or N’s ≥ GAPLIMIT will be filtered out.

  -h, --help                    Print help and exit
  -V, --version                 Print version and exit
      --input=filename          MSA input path
      --output=filename         Index/EFG output path
      --gap-limit=GAPLIMIT      Gap limit (suppressed by --elastic)
                                  (default=`1')
      --graphviz-output=filename
                                Graphviz output path
      --memory-chart-output=filename
                                Memory chart output path
  -e, --elastic                 Min-max-length semi-repeat-free segmentation
                                  (default=off)
      --gfa                     Saves output in xGFA format  (default=off)
  -p, --output-paths            Print the original sequences as paths of the
                                  xGFA graph (requires --gfa)  (default=off)
      --ignore-chars=STRING     Ignore these characters for the indexability
                                  property/pattern matching
  -t, --threads=THREADNUM       Max # threads  (default=`-1')
```

# todo
 - document EFG tricks related to option `--ignore-chars`, to the start and end of sequences, and to initial and ending runs of gaps
 - implement validation of .gfa files
 - implement pattern matching (`locate_multiple`, `locate_patterns`) on EFGs
 - implement min max height segmentation
