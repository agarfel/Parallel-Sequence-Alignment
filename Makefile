all: main

main: main.cpp
	g++ main.cpp -o main -fsanitize=address -g -O1


clean:
	rm -f main
