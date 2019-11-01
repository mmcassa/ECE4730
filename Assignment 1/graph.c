#include "graph.h"
#include <stdlib.h>
#include <stdio.h>


void read_graph(char *file_name, int *n, int ***A) {
    int i, j;
    char temp;
    FILE *file;

    file = fopen(file_name, "r");
    if (file == NULL)
    {
        printf("Error, invalid file name\n");  
        exit(0);
    }
    fread(&temp,1,1,file);
    *n = atoi(&temp);
    fread(&temp,1,1,file); // Skip newline
    
    *A = (int **) calloc(*n,sizeof(int *));
    for(i=0;i<*n;i++) {
        A[0][i] = (int *) calloc(*n,sizeof(int));
        for(j=0;j<*n;j++) {
            fread(&temp,1,1,file);
            if (temp == '-') {
                fread(&temp,1,1,file);
                A[0][i][j] = -1;
            } else {
                A[0][i][j] = atoi(&temp);
            }
            fread(&temp,1,1,file);

        }
    }

    fclose(file);
    
}

void write_graph(char *file_name, int n, int **A) {

    int i, j;
    FILE *file;

    file = fopen(file_name, "w");
    if (file == NULL)
    {
        printf("Error, invalid file name\n");  
        exit(0);
    }

    fprintf(file,"%d\n",n);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            fprintf(file, "%d ", A[i][j]);
        }
    }
    fprintf(file,"\n");
    fclose(file);
}

void print_graph(int n, int **A) {
  
  int i, j;

    // Loop through the number of rows + 2 because the first two rows are special
  for (i = 0; i < (n + 2); i++)
  {
    if (i == 0)
    {
      printf("    |");
    }
    else if (i == 1)
    {
      printf("    |");
    }
    else
    {
      printf("%4d|",i-2);
    }
    // Populate the columns
    for (j = 0; j < n; j++)
    {
      // If the first row
      if (i == 0)
      {
          printf("%4d ", j);
      }
      // If the second row
      else if (i == 1)
      {
        printf("-----");
      }
      else
      {
        printf("%4d ", A[i-2][j]);
      }
    }
    printf("\n");
  }

}
