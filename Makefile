include local.mk

CXXFLAGS	+= -std=c++17 -O2 -Wall
CPPFLAGS	+= -I$(SDSL_ROOT)/include -I$(SDSL_ROOT)/external/libdivsufsort/include
LDFLAGS		+= -L$(SDSL_ROOT)/lib -L$(SDSL_ROOT)/external/libdivsufsort/lib -lsdsl -ldivsufsort -ldivsufsort64


founderblockgraph_objects = founderblockgraph.o founder_block_index.o


all: founderblockgraph

clean:
	$(RM) founderblockgraph $(founder_block_index)

founderblockgraph: $(founderblockgraph_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(founderblockgraph_objects)

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<
