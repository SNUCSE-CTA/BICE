PROGRAMS = BICE BICE_100k BICE_huge_100k BICE_debug BICE_100k_debug

COMPILER = g++

BICE: main.cpp
	$(COMPILER) -w -std=c++11 -O3 -DNDEBUG main.cpp -o BICE

BICE_100k: main.cpp
	$(COMPILER) -w -std=c++11 -O3 -DNDEBUG main.cpp -o BICE_100k -DSUBGRAPH_MATCHING

BICE_huge_100k: main.cpp
	$(COMPILER) -w -std=c++11 -O3 -DNDEBUG main.cpp -o BICE_huge_100k -DSUBGRAPH_MATCHING -DHUGE_GRAPH

BICE_debug: main.cpp
	$(COMPILER) -w -std=c++11 -O3 -DDEBUG main.cpp -o BICE_debug

BICE_100k_debug: main.cpp
	$(COMPILER) -w -std=c++11 -O3 -DDEBUG main.cpp -o BICE_100k_debug -DSUBGRAPH_MATCHING

clean:
	rm -f $(PROGRAMS)
	rm -f *.o
