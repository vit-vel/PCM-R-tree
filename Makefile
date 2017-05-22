all:
	g++ src/main.cpp -std=c++11 -O2 -o run-rtree
debug:
	g++ src/main.cpp -std=c++11 -Og -o run-rtree-debug
clean:
	rm -f run-rtree
