CFLAGS= -O3 -g -fpermissive

test_as215: test_as215.cpp as215.o
	g++ -o $@ $^ -lm -lf2c $(CFLAGS)
