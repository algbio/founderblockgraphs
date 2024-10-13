GENGETOPT	?= gengetopt

OPT_FLAGS	?= -Ofast -march=native
#OPT_FLAGS	?= -O0 -g
CFLAGS		+= $(OPT_FLAGS) -std=c99 -Wall
CXXFLAGS	+= $(OPT_FLAGS) -pthread -std=c++17 -Wall
CPPFLAGS	+= -I ./sdsl-lite-v3/include

founderblockgraph_objects = founderblockgraph_cmdline.o founderblockgraph.o founder_block_index.o
locate_patterns_objects = locate_patterns_cmdline.o locate_patterns.o founder_block_index.o
locate_multiple_objects = locate_multiple.o founder_block_index.o


all: founderblockgraph locate_patterns locate_multiple

clean:
	$(RM) founderblockgraph locate_patterns $(founderblockgraph_objects) $(locate_patterns_objects) $(locate_multiple_objects)

founderblockgraph: $(founderblockgraph_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(founderblockgraph_objects)

locate_patterns: $(locate_patterns_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(locate_patterns_objects)

locate_multiple: $(locate_multiple_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(locate_multiple_objects)

founderblockgraph.cc: cmdline.c

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<" -F $*
