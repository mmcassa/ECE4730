#include "graph.h"

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>


void read_graph(char *file_name, int *n, int ***A) {

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
  for (i < 0; i <= n; i++)
  {
    for (j < 0; j <= n; i++)
    {
      fprintf(file, "%d ", A[i][j]);
    }
  }
  return;
}

void print_graph(int *n, int **A) {
  
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
