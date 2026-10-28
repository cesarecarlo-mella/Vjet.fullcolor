# Zgqq_evaluator: build the evaluator (C++17, header-only Boost.Odeint included) and run the tests
#   make                 build core/cpp/tridents_bs
#   make test            installation test (about 30 s)
#   make CXX=g++ PYTHON=/path/to/python3 test
CXX ?= c++
PYTHON ?= python3
FLAGS = -std=c++17 -O3 -fno-math-errno -fno-trapping-math -fassociative-math -fno-signed-zeros -DNDEBUG -Wno-deprecated-declarations

all: core/cpp/tridents_bs

core/cpp/tridents_bs: core/cpp/tridents_bs.cpp core/cpp/third_party/boost
	$(CXX) $(FLAGS) -I core/cpp/third_party -o $@ core/cpp/tridents_bs.cpp

core/cpp/third_party/boost: core/cpp/third_party.tar.gz
	cd core/cpp && tar xzf third_party.tar.gz && touch third_party/boost

test: all
	$(PYTHON) tests/run_tests.py

clean:
	rm -rf core/cpp/tridents_bs core/cpp/third_party core/out/generators_bs_y.bin

.PHONY: all test clean
