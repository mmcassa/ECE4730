#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char* argv[]) {
    // Initialize variables
    int i,j,k;
    int n;          // Matrix will be size n x n 
    int r;          // values within matrix will be 0 to r
    int p;          // modulus value
    int u;          // random number generated
    char *file_name;
    int **matrix;   // storage before export values
    char *args[] = {"-n","-r","-p","-o"};

    // Parse out arguemnets
    if (argc != 9) {
        printf("Usage:  make-graph -n {int} -r {int} -p {int} -o {file_name}\n");
        exit(1);
    }
    // Loop through every other arguement
    for(i=1;i<9;i+=2) {

        // Loop to check which command
        for(j=0;j<4;j++) {
            if (argv[i] == args[j])
                break;
        }
        // Determine which arguement is being processed
        switch (j)
                {
                case 0:
                    n = atoi(argv[i+1]);
                    break;

                case 1:
                    r = atoi(argv[i+1]);
                    break;

                case 2:
                    p = atoi(argv[i+1]);
                    break;
                
                case 3:
                    
                    break;
                
                default:

                    break;
                }
    }

    return 0;
}
