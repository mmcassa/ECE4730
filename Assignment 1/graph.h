#ifndef _graph_h_
#define _graph_h_

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) \
                     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
                     
void read_graph(char *file_name, int *n, int ***A);
void write_graph(char *file_name, int n, int **A);
void print_graph(int n, int **A);
#endif 