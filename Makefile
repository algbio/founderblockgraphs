include local.mk

GENGETOPT	?= gengetopt

OPT_FLAGS	?= -O2
CFLAGS		+= $(OPT_FLAGS) -std=c99 -Wall
CXXFLAGS	+= $(OPT_FLAGS) -std=c++17 -Wall
CPPFLAGS	+= -I$(SDSL_ROOT)/include -I$(SDSL_ROOT)/external/libdivsufsort/include
STATIC_LIBRARIES = $(SDSL_ROOT)/external/libdivsufsort/lib/libdivsufsort.a $(SDSL_ROOT)/external/libdivsufsort/lib/libdivsufsort64.a $(SDSL_ROOT)/lib/libsdsl.a


founderblockgraph_objects = founderblockgraph_cmdline.o founderblockgraph.o founder_block_index.o
locate_patterns_objects = locate_patterns_cmdline.o locate_patterns.o founder_block_index.o


all: founderblockgraph locate_patterns

clean:
	$(RM) founderblockgraph locate_patterns $(founderblockgraph_objects) $(locate_patterns_objects)

founderblockgraph: $(founderblockgraph_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(founderblockgraph_objects) $(STATIC_LIBRARIES)

locate_patterns: $(locate_patterns_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(locate_patterns_objects) $(STATIC_LIBRARIES)

founderblockgraph.cc: cmdline.c

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

#%.c: %.ggo
#	$(GENGETOPT) --input="$<" -F $*
