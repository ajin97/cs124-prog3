all: kk

kk: kk.o
	g++ kk.o -o kk -O3

kk.o: kk.cpp
	g++ -Wall -c kk.cpp -O3
