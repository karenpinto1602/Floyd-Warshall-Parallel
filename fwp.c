#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
//#include <time.h>
const int INFINITY = 1000000;

void Read_matrix(int localmat[], int n, int myrank, int p, MPI_Comm comm); //this function reads matrix data from local matrix.

void Print_matrix(int localmat[], int n, int myrank, int p, MPI_Comm comm); //prints the matrix

void Floyd(int localmat[], int n, int myrank, int p, MPI_Comm comm); //implemets Floyd algorithm in MPI

int Owner(int k, int p, int n);

void Copy_row(int localmat[], int n, int p, int ro_wk[], int k);

void Printrow(int localmat[], int n, int myrank, int i); //prints row of the matrix

int main(int argc, char* argv[]) {

//clock_t t;
//t = clock();
        int n;
        int* localmat;
        MPI_Comm comm;
        int p, myrank;
        MPI_Init(&argc, &argv);
        comm = MPI_COMM_WORLD;
        MPI_Comm_size(comm, &p);
        MPI_Comm_rank(comm, &myrank);
        if (myrank == 0) {
                printf("How many vertices:-\n");
                scanf("%d", &n);
        }
        MPI_Bcast(&n, 1, MPI_INT, 0, comm);
        localmat = malloc(n*n/p*sizeof(int));
        if (myrank == 0)
                printf("Enter the localmatrix\n");
        Read_matrix(localmat, n, myrank, p, comm);

        if (myrank == 0)
                printf("We got\n");
        Print_matrix(localmat, n, myrank, p, comm);
        if (myrank == 0)
                printf("\n");
        Floyd(localmat, n, myrank, p, comm);

        if (myrank == 0)
                printf("The solution is:\n");
        Print_matrix(localmat, n, myrank, p, comm);

        free(localmat);
        MPI_Finalize();
//t = clock()-t;
//double time_taken = ((double)t)/CLOCKS_PER_SEC;
//if(myrank == 0){
//printf("\nTime taken is %f seconds\n",time_taken);
//}
        return 0;
}

void Read_matrix(int localmat[], int n, int myrank, int p, MPI_Comm comm) {
        int i, j;
        int* tempmat = NULL;
        if (myrank == 0) {
                tempmat = malloc(n*n*sizeof(int));
        for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                        scanf("%d", &tempmat[i*n+j]);
        MPI_Scatter(tempmat, n*n/p, MPI_INT,localmat, n*n/p, MPI_INT, 0, comm);
        free(tempmat);
        } else {
                MPI_Scatter(tempmat, n*n/p, MPI_INT,localmat, n*n/p, MPI_INT, 0, comm);
        }
}

void Print_matrix(int localmat[], int n, int myrank, int p, MPI_Comm comm) {
        int i, j;
        int* tempmat = NULL;
        if (myrank == 0) {
                tempmat = malloc(n*n*sizeof(int));
                MPI_Gather(localmat, n*n/p, MPI_INT,
                tempmat, n*n/p, MPI_INT, 0, comm);
                for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++)
                                if (tempmat[i*n+j] == INFINITY)
                                        printf("i ");
                                else
                                        printf("%d ", tempmat[i*n+j]);
                        printf("\n");
                }
                free(tempmat);
                } else {
                MPI_Gather(localmat, n*n/p, MPI_INT,
                tempmat, n*n/p, MPI_INT, 0, comm);
        }
}

void Floyd(int localmat[], int n, int myrank, int p, MPI_Comm comm) {
        int globalk, locali, globalj, temp;
        int root;
        int* rowk = malloc(n*sizeof(int));
        for (globalk = 0; globalk < n; globalk++) {
                root = Owner(globalk, p, n);
                if (myrank == root)
                        Copy_row(localmat, n, p, rowk, globalk);
                MPI_Bcast(rowk, n, MPI_INT, root, comm);
                for (locali = 0; locali < n/p; locali++){
                for (globalj = 0; globalj < n; globalj++) {
                                temp = localmat[locali*n + globalk] + rowk[globalj];
                                if (temp < localmat[locali*n+globalj])
                                        localmat[locali*n + globalj] = temp;
                        }
        }}
        free(rowk);
}

int Owner(int k, int p, int n) {
        return k/(n/p);
}

void Copy_row(int localmat[], int n, int p, int rowk[], int k) {
        int j;
        int localk = k % (n/p);
        for (j = 0; j < n; j++)
                rowk[j] = localmat[localk*n + j];
}

