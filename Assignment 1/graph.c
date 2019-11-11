/* Mitchell Cassaday & William (Luke) Benton
 * ECE4730
 * Fall 2019
 * 4 NOV 2019
 * Goal: Helper functions
 */

#include "graph.h"
#include <stdlib.h>
#include <stdio.h>

#define MAXCHAR 15

void read_graph(char *file_name, int *n, int ***A) {
    int i, j;
    //char temp[MAXCHAR];
    FILE *file;

    file = fopen(file_name, "r");
    if (file == NULL)
    {
        printf("Error, invalid file name\n");  
        exit(0);
    }
    // Read in n
    fscanf (file, "%d\n", n);
    
    *A = (int **) calloc(*n,sizeof(int *));
    for(i=0;i<*n;i++) {
        A[0][i] = (int *) calloc(*n,sizeof(int));
        for(j=0;j<*n;j++) {
            fscanf(file, "%d ",&(A[0][i][j]));
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

void print_graph2(int n,int m, int **A) {
  
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
    for (j = 0; j < m; j++)
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