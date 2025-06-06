all: main simple

main: main.cpp
	g++ main.cpp -o main -fsanitize=address -g -O1

simple : simple.cpp
	g++ simple.cpp -o simple

clean:
	rm -f main simple
