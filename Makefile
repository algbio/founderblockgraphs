include local.mk

CXXFLAGS	+= -std=c++17 -O2 -Wall
CPPFLAGS	+= -I$(SDSL_ROOT)/include -I$(SDSL_ROOT)/external/libdivsufsort/include
LDFLAGS		+= -L$(SDSL_ROOT)/lib -L$(SDSL_ROOT)/external/libdivsufsort/lib -lsdsl -ldivsufsort -ldivsufsort64


founderblockgraph_objects = founderblockgraph.o founder_block_index.o
locate_patterns_objects = locate_patterns.o founder_block_index.o


all: founderblockgraph locate_patterns

clean:
	$(RM) founderblockgraph locate_patterns $(founderblockgraph_objects) $(locate_patterns_objects)

founderblockgraph: $(founderblockgraph_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(founderblockgraph_objects)

locate_patterns: $(locate_patterns_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(locate_patterns_objects)

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<
