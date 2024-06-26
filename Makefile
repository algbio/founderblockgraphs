GENGETOPT	?= gengetopt

OPT_FLAGS	?= -Ofast -march=native
#OPT_FLAGS	?= -O0 -g
CFLAGS		+= $(OPT_FLAGS) -std=c99 -Wall
CXXFLAGS	+= $(OPT_FLAGS) -pthread -std=c++17 -Wall
SDSL_ROOT	= ./sdsl-lite/build
CPPFLAGS	+= -isystem $(SDSL_ROOT)/include -isystem $(SDSL_ROOT)/external/libdivsufsort/include
STATIC_LIBRARIES = $(SDSL_ROOT)/external/libdivsufsort/lib/libdivsufsort.a $(SDSL_ROOT)/external/libdivsufsort/lib/libdivsufsort64.a $(SDSL_ROOT)/lib/libsdsl.a


founderblockgraph_objects = founderblockgraph_cmdline.o founderblockgraph.o founder_block_index.o
locate_patterns_objects = locate_patterns_cmdline.o locate_patterns.o founder_block_index.o
locate_multiple_objects = locate_multiple.o founder_block_index.o


all: founderblockgraph locate_patterns locate_multiple

clean:
	$(RM) founderblockgraph locate_patterns $(founderblockgraph_objects) $(locate_patterns_objects) $(locate_multiple_objects)

founderblockgraph: $(founderblockgraph_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(founderblockgraph_objects) $(STATIC_LIBRARIES)

locate_patterns: $(locate_patterns_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(locate_patterns_objects) $(STATIC_LIBRARIES)

locate_multiple: $(locate_multiple_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(locate_multiple_objects) $(STATIC_LIBRARIES)

founderblockgraph.cc: cmdline.c

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<" -F $*
