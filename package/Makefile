# Zgqq_evaluator: build the evaluator (C++17, header-only Boost.Odeint included) and run the tests
#   make                 build core/cpp/tridents_bs and core/cpp/iter_quad (quad precision)
#   make test            installation test (about 30 s)
#   make test-quad       test of the quad precision option (about 1 min)
#   make CXX=g++ PYTHON=/path/to/python3 test
CXX ?= c++
PYTHON ?= python3
FLAGS = -std=c++17 -O3 -fno-math-errno -fno-trapping-math -fassociative-math -fno-signed-zeros -DNDEBUG -Wno-deprecated-declarations

# double-double arithmetic relies on exact rounding: no fast-math, no fused contractions
QFLAGS = -std=c++17 -O3 -ffp-contract=off

all: core/cpp/tridents_bs core/cpp/iter_quad

core/cpp/iter_quad: core/cpp/iter_quad.cpp core/cpp/dd.hpp
	$(CXX) $(QFLAGS) -o $@ core/cpp/iter_quad.cpp

core/cpp/tridents_bs: core/cpp/tridents_bs.cpp core/cpp/third_party/boost
	$(CXX) $(FLAGS) -I core/cpp/third_party -o $@ core/cpp/tridents_bs.cpp

core/cpp/third_party/boost: core/cpp/third_party.tar.gz
	cd core/cpp && tar xzf third_party.tar.gz && touch third_party/boost

test: all
	$(PYTHON) tests/run_tests.py

test-quad: all
	$(PYTHON) tests/run_tests.py --quad

clean:
	rm -rf core/cpp/tridents_bs core/cpp/iter_quad core/cpp/third_party core/out/generators_bs_y.bin

.PHONY: all test test-quad clean
