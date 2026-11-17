CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
all: bin/arc_tests bin/arc_probe bin/arc_bench
bin/%: %.cpp arcdisp.hpp gen.hpp ../convex/kdisp.hpp
	@mkdir -p bin && $(CXX) $(CXXFLAGS) -o $@ $<
clean:
	rm -rf bin
.PHONY: all clean
