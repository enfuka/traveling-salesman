#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <string.h>
#include <sys/times.h>
#define mat_dim 16
double start, stop, used, mf;
double ftime(void);
double ftime(void)
{
    struct tms t;
    times(&t);
    return (t.tms_utime + t.tms_stime) / 100.0;
}
int main(int argc, char *argv[])
{
    int i, j, k, l, m, prID, Root = 0, localn;
    int **A, **B, **C, *a, *b, *c, *Temp_A, *Temp_B, *Temp_C;
    int dims[2], periods[2], coords[2], exdims[2], rank, Size, sqrtProcs;
    MPI_Comm comm_2d, comm_row, comm_col;
    MPI_Status status;
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &(Size));
    MPI_Comm_rank(MPI_COMM_WORLD, &(rank));

    sqrtProcs = (int)sqrt((double)Size);
    dims[0] = dims[1] = sqrtProcs;
    periods[0] = periods[1] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &(comm_2d));
    MPI_Cart_coords(comm_2d, rank, 2, coords); // coords: 2D coordinates of the process eg. rank 0 -> (0,0) rank 1 -> (0,1) rank 2 -> (1,0) rank 3 -> (1,1)

    exdims[0] = 0;
    exdims[1] = 1;
    MPI_Cart_sub(comm_2d, exdims, &(comm_row));
    exdims[0] = 1;
    exdims[1] = 0;
    MPI_Cart_sub(comm_2d, exdims, &(comm_col));
    if (rank == Root)
    {
        A = (int **)malloc(mat_dim * sizeof(int *));
        B = (int **)malloc(mat_dim * sizeof(int *));
        C = (int **)malloc(mat_dim * sizeof(int *));
        for (i = 0; i < mat_dim; i++)
        {
            A[i] = (int *)malloc(mat_dim * sizeof(int));
            B[i] = (int *)malloc(mat_dim * sizeof(int *));
            C[i] = (int *)malloc(mat_dim * sizeof(int *));
            for (j = 0; j < mat_dim; j++)
            {
                A[i][j] = 1;
                B[i][j] = 2;
            }
        }
        start = ftime();
    }
    localn = mat_dim / sqrtProcs;
    a = (int *)malloc(localn * localn * sizeof(int));
    b = (int *)malloc(localn * localn * sizeof(int));
    /* converting into 1-d array */
    Temp_A = (int *)malloc(sizeof(int) * mat_dim * mat_dim);
    Temp_B = (int *)malloc(sizeof(int) * mat_dim * mat_dim);
    if (rank == Root)
    { /* Matrix A and B arrangements*/
        for (k = 0; k < sqrtProcs; k++)
        {
            for (l = 0; l < sqrtProcs; l++)
            {
                prID = k * sqrtProcs + l;
                for (i = 0; i < localn; i++)
                    for (j = 0; j < localn; j++)
                    {
                        Temp_A[(prID * localn * localn) + (i * localn) + j] = A[k * localn + i][l * localn + j];
                        Temp_B[(prID * localn * localn) + (i * localn) + j] = B[k *
                                                                                    localn +
                                                                                i][l * localn + j];
                    }
            }
        }
    }

    MPI_Barrier(comm_2d);
    MPI_Scatter(Temp_A, localn * localn, MPI_INT, a, localn * localn, MPI_INT, 0,
                comm_2d);
    MPI_Scatter(Temp_B, localn * localn, MPI_INT, b, localn * localn, MPI_INT, 0,
                comm_2d);

    /* Initial arrangement of Matrices */
    if (coords[0] != 0)
        MPI_Sendrecv_replace(a, localn * localn, MPI_INT, ((coords[1] + sqrtProcs - coords[0]) % sqrtProcs), 0, ((coords[1] + coords[0]) % sqrtProcs), 0, comm_row,
                             &status);
    if (coords[1] != 0)
        MPI_Sendrecv_replace(b, localn * localn, MPI_INT, ((coords[0] + sqrtProcs - coords[1]) % sqrtProcs), 0, ((coords[0] + coords[1]) % sqrtProcs), 0, comm_col,
                             &status);
    c = (int *)malloc(localn * localn * sizeof(int));
    for (i = 0; i < localn * localn; i++)
        c[i] = 0;

    /* Canon's Algo */
    for (l = 0; l < sqrtProcs; l++)
    {
        m = 0;
        for (i = 0; i < localn; i++)
        {
            for (j = 0; j < localn; j++)
            {
                for (k = 0; k < localn; k++)
                    c[m] += a[i * localn + k] * b[k * localn + j];
                m++;
            }
        }
        /* Move by one position Blocks of A leftwards and that of B upwards with
        wraparound */
        MPI_Sendrecv_replace(a, localn * localn, MPI_INT, ((coords[1] + sqrtProcs - 1) % sqrtProcs), 0, ((coords[1] + 1) % sqrtProcs), 0, comm_row, &status);
        MPI_Sendrecv_replace(b, localn * localn, MPI_INT, ((coords[0] + sqrtProcs - 1) % sqrtProcs), 0, ((coords[0] + 1) % sqrtProcs), 0, comm_col, &status);
    }
    /* 1-d array buffer for C matrix */
    if (rank == Root)
        Temp_C = (int *)malloc(sizeof(int) * mat_dim * mat_dim);
    MPI_Barrier(comm_2d);
    MPI_Gather(c, localn * localn, MPI_INT, Temp_C, localn * localn, MPI_INT, Root,
               comm_2d);
    if (rank == Root)
    { /* Arrange output matrix C */
        for (k = 0; k < sqrtProcs; k++)
        {
            for (l = 0; l < sqrtProcs; l++)
            {
                prID = k * sqrtProcs + l;
                for (i = 0; i < localn; i++)
                    for (j = 0; j < localn; j++)
                        C[k * localn + i][l * localn + j] = Temp_C[(prID * localn *
                                                                    localn) +
                                                                   (i * localn) + j];
            }
        }
        /*
        //Print results...........
        printf("Matrix A \n");
        for(i = 0; i < mat_dim; i++) {
        for(j = 0; j < mat_dim; j++)
        printf ("%7d ", A[i][j]);
        printf ("\n");
        }
        printf("\n");
        printf("Matrix B \n");
        for(i = 0; i < mat_dim; i++){
        for(j = 0; j < mat_dim; j++)
        printf("%7d ", B[i][j]);
        printf("\n");
        }
        printf("\n");
        printf("Matrix C \n");
        for(i = 0; i < mat_dim; i++){
        for(j = 0; j < mat_dim; j++)
        printf("%7d ",C[i][j]);
        printf("\n");
        }*/
        stop = ftime();
        used = stop - start;
        mf = ((mat_dim * 2.0) / used / 1000000.0) * mat_dim * mat_dim;
        printf("\n Elapsed time : %10.2f\n", used);
        printf("\n DP MFLOPS : %10.2f\n", mf);
    }
    MPI_Finalize();
    return 0;
}