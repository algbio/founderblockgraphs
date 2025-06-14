GENGETOPT	?= gengetopt

OPT_FLAGS	?= -Ofast -march=native
#OPT_FLAGS	?= -O0 -g
CFLAGS		+= $(OPT_FLAGS) -std=c99 -Wall
CXXFLAGS	+= $(OPT_FLAGS) -pthread -std=c++17 -Wall -lc
CPPFLAGS	+= -I ./sdsl-lite-v3/include

founderblockgraph_objects = founderblockgraph_cmdline.o founderblockgraph.o founder_block_index.o


all: founderblockgraph

clean:
	$(RM) founderblockgraph $(founderblockgraph_objects)

founderblockgraph: $(founderblockgraph_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(founderblockgraph_objects)

founderblockgraph.cc: cmdline.c

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<" -F $*
