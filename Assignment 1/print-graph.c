/* Mitchell Cassaday & William (Luke) Benton
 * ECE4730
 * Fall 2019
 * 4 NOV 2019
 * Goal: Displays a graph from a file
 */

#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

int main(int argc, char* argv[]) {
    int n;
    int **A;
    char *file_name = calloc(100,sizeof(char));


    // Parse out arguemnets
    if (argc != 2) {
        printf("Usage:  print-graph {file_name}\n");
        exit(1);
    }

    file_name = argv[1];

    read_graph(file_name,&n,&A);
    print_graph(n,A);

    return 0;
}