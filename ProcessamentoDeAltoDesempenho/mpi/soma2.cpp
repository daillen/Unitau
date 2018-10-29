#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

void print(int n, int m, int **arr)
{
  for(int i=0; i<n; i++) {
    for(int j=0; j<m; j++) {
      printf("%d ", arr[i][j]);
    }
    printf("\n");
  }
}

int main(int argc, char * argv[])
{
  int **a, **b, **c, *buffer, *buffera, *bufferb;
  int qtde, i, j, k, local, N, M;
  
  
  int mpi_rank, mpi_size, mpi_tag=100;
  MPI_Status  mpi_status;
  
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  if(mpi_rank == 0) {
    scanf("%d %d", &N, &M);
    // Somente 0 possui os vetores a, b e c globais
    a = (int**)malloc(N*sizeof(int*));
    b = (int**)malloc(N*sizeof(int*));
    c = (int**)malloc(N*sizeof(int*));

    // Zero carrega a, b e c
    for(i=0; i<N; i++) {
      a[i] = (int*)malloc(M*sizeof(int));
      b[i] = (int*)malloc(M*sizeof(int));
      c[i] = (int*)malloc(M*sizeof(int));
    }

    for(int i=0; i<N; i++) {
      for(int j=0; j<M; j++) {
	scanf("%d", &a[i][j]);
      }
    }
    
    for(int i=0; i<N; i++) {
      for(int j=0; j<M; j++) {
	scanf("%d", &b[i][j]);
      }
    }

    //print(N, M, a);
    //printf("\n");
    //print(N, M, b);

    int n_local = N / mpi_size;
    int remainder = N % mpi_size;
    
    buffera = (int*)malloc((N+1)*sizeof(int));
    bufferb = (int*)malloc((N+1)*sizeof(int));
    int j = 0, idx = 0;
    for(int i=1; i<mpi_size; i++) {
      int n_to_send = n_local + (remainder-- > 0 ? 1 : 0);
      int local_size = n_to_send * M;
      MPI_Send(&local_size, 1, MPI_INT, i, mpi_tag, MPI_COMM_WORLD);

      idx = 0, j=0;
      for(int cnt=0; cnt < n_to_send; j++, cnt++) {
	for(int k=0; k<M; k++, idx++) {
	  buffera[idx] = a[j][k];
	}
      }

      MPI_Send(&buffera, local_size, MPI_INT, i, mpi_tag, MPI_COMM_WORLD);
      
      j = 1; idx = 0;
      for(int cnt=0; cnt < n_to_send; j++, cnt++) {
	for(int k=0; k<M; k++, idx++) {
	  bufferb[idx] = b[j][k];
	  printf("%d ", buffer[idx]);
	}
      }
      MPI_Send(&bufferb, local_size, MPI_INT, i, mpi_tag, MPI_COMM_WORLD);
      
    }
  } else {
   
    MPI_Recv(&N, 1, MPI_INT, 0, mpi_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    int *a, *b, *c;
    a = (int*)malloc((N+1)*sizeof(int));
    b = (int*)malloc((N+1)*sizeof(int));
    c = (int*)malloc((N+1)*sizeof(int));
    for(int i = 0; i<N; i++){
      printf("%d\n", a[i]);
    }

    MPI_Recv(&a, N, MPI_INT, 0, mpi_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&b, N, MPI_INT, 0, mpi_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    for(int i=0; i<N; i++) {
      printf("%d\n", a[i]);
    }

    for(int i=0; i<N; i++) {
      c[i] = a[i] + b[i];
    }
 
    for(int i=0; i<N; i++) {
      printf("%d\n", c[i]);
    }
  
    MPI_Send(&c, N, MPI_INT, 0, mpi_tag, MPI_COMM_WORLD);
  }

  if (mpi_rank == 0) {
    int n_local = N / mpi_size;
    int remainder = N % mpi_size;

    int j = 0, idx, cnt;
    for(int i=1; i<mpi_size; i++) {
      int n_to_recv = n_local + (remainder-- > 0 ? 1 : 0);
      int local_size = n_to_recv * M;

      MPI_Recv(&buffer, local_size, MPI_INT, i, mpi_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      idx = 0;
      for(int cnt=0; cnt < n_to_recv; j++, cnt++) {
	for(int k=0; k<M; k++, idx++) {
	  c[j][k] = buffer[idx];
	}
      }

    }
    
    for(; j < N; j++) {
      for(int k=0; k<M; k++) {
	c[j][k] = a[j][k] + b[j][k];
      }
    }
    
    print(N, M, c);
  }

  MPI_Finalize();
}
