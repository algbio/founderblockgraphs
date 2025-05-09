GENGETOPT	?= gengetopt
OPT_FLAGS	?= -Ofast -march=native
#OPT_FLAGS	?= -O0 -g
CFLAGS		+= $(OPT_FLAGS) -std=c99 -Wall
CXXFLAGS	+= $(OPT_FLAGS) -std=c++20 -Wall -lc -lz

locate_patterns_objects = src/command-line-parsing/locate_patterns.o src/locate_patterns.o
locate_patterns_deps = lib/gfakluge/src/gfakluge.hpp lib/seqtk/kseq.h src/gafanchor.hpp
locate_patterns_CXXFLAGS = -I ./lib/sdsl-lite-v3/include -I ./lib/gfakluge/src -I ./lib/gfakluge/src/tinyFA -isystem ./lib/gfakluge/src/tinyFA/pliib -I ./lib/seqtk

all: locate_patterns

clean:
	$(RM) locate_patterns $(locate_patterns_objects)

locate_patterns: $(locate_patterns_objects) $(locate_patterns_deps)
	$(CXX) $(CXXFLAGS) -o $@ $(locate_patterns_objects)

src/locate_patterns.o: src/locate_patterns.cpp $(locate_patterns_deps)
	$(CXX) -c $(locate_patterns_CXXFLAGS) $(CXXFLAGS) -o $@ $<

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --unnamed-opts --input="$<" -F $(basename $(notdir $@)) --output-dir $(dir $@)
