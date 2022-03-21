all: diet

diet: main.cpp
	g++ -g -o diet main.cpp
valgrind:
	valgrind --tool=callgrind ./diet
	kcachegrind callgrind.out.*
.PHONY: clean
clean:
	rm diet
	rm callgrind.out.*
