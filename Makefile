all: complex simple

complex: complex.cpp
	g++ complex.cpp -o complex

simple : simple.cpp
	g++ simple.cpp -o simple

clean:
	rm -f complex simple
