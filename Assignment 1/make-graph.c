#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

int main(int argc, char* argv[]) {
    // Initialize variables
    int i,j,k;
    int n,n_chk = 0;          // Matrix will be size n x n 
    int r,r_chk = 0;          // values within matrix will be 0 to r
    int p,p_chk = 0;          // modulus value
    int u;          // random number generated
    int o_chk = 0;
    char *file_name = calloc(100,sizeof(char));
    int **matrix;   // storage before export values
    char *args[] = {"-n","-r","-p","-o"};
    time_t t;
   
   
   /* Intializes random number generator */
   srand((unsigned) time(&t));


    // Parse out arguemnets
    if (argc != 9) {
        printf("Usage:  make-graph -n {int} -r {int} -p {int} -o {file_name}\n");
        exit(1);
    }
    // Loop through every other arguement
    for(i=1;i<9;i+=2) {

        // Loop to check which command
        for(j=0;j<4;j++) {
            if (strcmp(argv[i],args[j]) == 0)
                break;
        }
        // Determine which arguement is being processed
        switch (j)
            {
            case 0:
                if (n_chk++ == 1) {
                    printf("Usage:  Do not repeat a tag\n");
                    exit(1);
                }
                n = atoi(argv[i+1]);
                break;

            case 1:
                if (r_chk++ == 1) {
                    printf("Usage:  Do not repeat a tag\n");
                    exit(1);
                }
                r = atoi(argv[i+1]);
                break;

            case 2:
                if (p_chk++ == 1) {
                    printf("Usage:  Do not repeat a tag\n");
                    exit(1);
                }
                p = atoi(argv[i+1]);
                break;
            
            case 3:
                if (o_chk++ == 1) {
                    printf("Usage:  Do not repeat a tag\n");
                    exit(1);
                }
                file_name = argv[i+1];
                break;
            
            default:
                printf("Usage:  make-graph -n {int} -r {int} -p {int} -o {file_name}\n");
                exit(1);
                break;
        }
    }

    // Matrix calloc
    matrix = (int **) calloc(n,sizeof(int *));
    for (i=0;i<n;i++) {
        matrix[i] = (int *) calloc(n,sizeof(int));
    }
    // Start the random number generation
    for (i=0;i<n;i++) {
        //matrix[i] = calloc(n,sizeof(int));
        for (j=0;j<n;j++) {
            if (i != j) {
                u = (rand() % (p-1)) + 1;
                if (u <= r) {
                    matrix[i][j] = u;
                } else {
                    matrix[i][j] = -1;
                }

            } else {
                matrix[i][j] = 0;
            }
        }
    }

    // Export graph to file
    write_graph(file_name,n,matrix);
    return 0;
}
