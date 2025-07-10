GENGETOPT	?= gengetopt

OPT_FLAGS	?= -O3 -march=native
CFLAGS		+= $(OPT_FLAGS) -std=c99 -Wall
CXXFLAGS	+= $(OPT_FLAGS) -pthread -std=c++20 -Wall -lc -lz
CPPFLAGS	+= -I ./lib/sdsl-lite-v3/include -isystem ./lib/seqtk

founderblockgraph_objects = src/command-line-parsing/founderblockgraph.o src/founderblockgraph.o
founderblockgraph_deps = src/utils.hpp src/algo.hpp src/index.hpp

all: founderblockgraph

clean:
	$(RM) founderblockgraph $(founderblockgraph_objects)

founderblockgraph: $(founderblockgraph_objects) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(founderblockgraph_objects)

founderblockgraph.cc: cmdline.c

%.o: %.cpp $(founderblockgraph_deps)
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

#%.c: %.ggo
#	$(GENGETOPT) --unnamed-opts --input="$<" -F $(basename $(notdir $@)) --output-dir $(dir $@)
