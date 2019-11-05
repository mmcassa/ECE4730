#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv)
{
  int  rank;
  int  size;
  char sbuf[4][8];
  void *rbuf;
  MPI_Status status;

  rbuf = (void *)calloc(1, sizeof(void));
  strcpy(sbuf[0], "msg1");
  strcpy(sbuf[1], "msg2");
  strcpy(sbuf[2], "msg3");
  strcpy(sbuf[3], "msg4");

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
    MPI_Send(sbuf[rank], 8, MPI_CHAR, size-1, 99, MPI_COMM_WORLD);
    printf("Message sent from Rank: %d to Rank: %d: %s\n", rank, size-1, sbuf[rank]);
  }
  else
  {
    MPI_Send(sbuf[rank], 8, MPI_CHAR, rank - 1, 99, MPI_COMM_WORLD);
    printf("Message sent from Rank: %d to Rank: %d: %s\n", rank, rank-1, sbuf[rank]);
  }

  if (rank != 3)
  {
    MPI_Recv(rbuf, 8, MPI_CHAR, rank + 1, 99, MPI_COMM_WORLD, &status);
    printf("Message received on Rank: %d from Rank: %d: %s\n", rank, rank+1, (char *)rbuf);
  }
  else
  {
    MPI_Recv(rbuf, 8, MPI_CHAR, 0, 99, MPI_COMM_WORLD, &status);
    printf("Message received on Rank: %d, from Rank: 0: %s\n", rank, (char *)rbuf);
  }

  MPI_Finalize();

  return 0;
}
