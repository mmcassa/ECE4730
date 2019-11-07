/* Mitchell Cassaday & William (Luke) Benton
 * ECE4730
 * Fall 2019
 * 4 NOV 2019
 * Goal: Apply floyd's algorithm without mpi
 */

#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

int main(int argc, char* argv[]) {
	clock_t t1, t2;
	double result1;
	double result2;
	t1 = clock();
    int i,j,k,n;
    int **A;
    int **floyd;
    char *file_in = calloc(100,sizeof(char));
    char *file_out = calloc(100,sizeof(char));

	
    // Parse out arguemnets
    if (argc != 3) {
        printf("Usage:  print-graph {file_name}\n");
        exit(1);
    }

    file_in = argv[1];
    file_out = argv[2];

    read_graph(file_in,&n,&A);
    //print_graph(n,A);

	t2 = clock();
    floyd = (int **) calloc(n,sizeof(int *));
    for (i=0;i<n;i++) {
        floyd[i] = (int *) calloc(n,sizeof(int));
    }
	
    for(k=0;k<n;k++) {
        for(i=0;i<n;i++) {
            for(j=0;j<n;j++) {
                if (A[i][k] > 0 && A[k][j] > 0) {
                    if ((A[i][k]+A[k][j]) < A[i][j] || A[i][j] == -1) {
                        A[i][j] = A[i][k]+A[k][j]; 
                    }
                }
            }
        }
    }
	t2 = clock() - t2;
	result1 = ((double) t2) / CLOCKS_PER_SEC;
    //printf("\n\n\n");
    //print_graph(n,A);
    write_graph(file_out,n,A);
	t1 = clock() - t1;
	result2 = ((double) t1)/ CLOCKS_PER_SEC;
	
	printf("Total Time: %f\n", result2);
	printf("Floyd's Time: %f\n", result1);
    return 0;
}

