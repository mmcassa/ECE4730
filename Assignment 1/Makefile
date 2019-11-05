CC = gcc
CFLAGS = -g -Wall

all: make-graph print-graph floyd-serial floyd-parallel

make-graph: make-graph.c graph.c graph.h
	$(CC) $(CFLAGS) -o make-graph make-graph.c graph.c

print-graph: print-graph.c graph.c graph.h
	$(CC) $(CFLAGS) -o print-graph print-graph.c graph.c

floyd-serial: floyd-serial.c graph.c graph.h
	$(CC) $(CFLAGS) -o floyd-serial floyd-serial.c graph.c

floyd-parallel: floyd-parallel.c graph.c graph.h
	$(CC) $(CFLAGS) -o floyd-parallel floyd-parallel.c graph.c

clean: 
	rm make-graph print-graph floyd-serial floyd-parallel
