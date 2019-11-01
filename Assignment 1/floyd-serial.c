#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int main(int argc, char* argv[]) {
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
    print_graph(n,A);

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
    printf("\n\n\n");
    print_graph(n,A);
    write_graph(file_out,n,A);
    return 0;
}

