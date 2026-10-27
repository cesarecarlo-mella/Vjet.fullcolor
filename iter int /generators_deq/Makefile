# tridents_bs (Bulirsch-Stoer evaluator of trident_deq, unchanged) for the generator system
CXX ?= clang++
FLAGS = -std=c++17 -O3 -fno-math-errno -fno-trapping-math -fassociative-math -fno-signed-zeros -DNDEBUG -Wno-deprecated-declarations

cpp/tridents_bs: cpp/tridents_bs.cpp cpp/third_party/boost
	$(CXX) $(FLAGS) -I cpp/third_party -o $@ cpp/tridents_bs.cpp

cpp/third_party/boost: cpp/third_party.tar.gz
	cd cpp && tar xzf third_party.tar.gz && touch third_party/boost
