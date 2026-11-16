CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
BIN = bin/kdisp bin/bench bin/hunt bin/selftest bin/fmm

all: $(BIN)
bin/kdisp: kdisp_cli.cpp kdisp.hpp
	@mkdir -p bin && $(CXX) $(CXXFLAGS) -o $@ $<
bin/bench: bench.cpp kdisp.hpp
	@mkdir -p bin && $(CXX) $(CXXFLAGS) -o $@ $<
bin/hunt: hunt.cpp kdisp.hpp randgen.hpp
	@mkdir -p bin && $(CXX) $(CXXFLAGS) -o $@ $<
bin/selftest: selftest.cpp kdisp.hpp randgen.hpp
	@mkdir -p bin && $(CXX) $(CXXFLAGS) -o $@ $<
# the general-dimension clique / matrix-product solver from the repository root
bin/fmm: ../main.cpp ../dispersion.cpp ../clique.cpp ../matrix.cpp
	@mkdir -p bin && $(CXX) -O2 -std=c++17 -o $@ $^

test: all
	bin/selftest 3000 1
	python3 experiments/check_counterexamples.py
clean:
	rm -rf bin
.PHONY: all test clean
