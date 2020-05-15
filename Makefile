include local.mk

GENGETOPT	?= gengetopt

OPT_FLAGS	?= -O2
CFLAGS		+= $(OPT_FLAGS) -std=c99 -Wall
CXXFLAGS	+= $(OPT_FLAGS) -std=c++17 -Wall
CPPFLAGS	+= -I$(SDSL_ROOT)/include -I$(SDSL_ROOT)/external/libdivsufsort/include
LDFLAGS		+= -L$(SDSL_ROOT)/lib -L$(SDSL_ROOT)/external/libdivsufsort/lib -lsdsl -ldivsufsort -ldivsufsort64


founderblockgraph_objects = cmdline.o founderblockgraph.o founder_block_index.o
locate_patterns_objects = locate_patterns.o founder_block_index.o


all: founderblockgraph locate_patterns

clean:
	$(RM) founderblockgraph locate_patterns $(founderblockgraph_objects) $(locate_patterns_objects)

founderblockgraph: $(founderblockgraph_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(founderblockgraph_objects)

locate_patterns: $(locate_patterns_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(locate_patterns_objects)

founderblockgraph.cc: cmdline.c

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<"
