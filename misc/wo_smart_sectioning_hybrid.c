// run using $ cd ..../ProgrammingEnv/project/
//           $ gcc -Wall clustering.c functions.c -o clustering -lm
// the -lm is to link to the math library
//           $ ./clustering > output.txt

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <memory.h>
#include "mpi.h"
#include "omp.h"
#include "functions.h"

//#define DIETAG (N+5000);
static void manager(double *A);
static void worker(int chunk_size);
static void sequential_naive(double *A, int num_rows);

int main(int argc, char *argv[]) {
    double A[NUM_VECS * N]; 
    double x_i[N];
    int i, j;
    
    // Sequential 
    srand(200);
    for (i = 0; i < NUM_VECS; i++) {
        for (j = 0; j < N; j++) {
            init_vec(x_i, N);
            A[i * N + j] = *x_i;
        }
    }
    
    //int t = omp_get_num_procs();
    //omp_set_num_threads(t);     


    int my_rank;
    int num_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    double start_time = 0;
    double end_time = 0;   
    
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    //printf("Number of processes: %d\n", num_procs);
    
    if (num_procs > 1) {
        int chunk_size = NUM_VECS / (num_procs - 1); // not including process 0
        if (!my_rank) {
            manager(A);
        }else {
            worker(chunk_size);
        }
    }else {
        if (!my_rank) {
            sequential_naive(A, NUM_VECS);
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    MPI_Finalize();
    
    if (!my_rank) {
        double time_interval = end_time - start_time;
        printf("Total time (sec): %f\n", time_interval);
    }
    
    //printf("Total time (sec): %f for %d threads\n", time_interval, t);
    return 0;    
}



static void manager(double *A) {
    //int m = NUM_VECS;//+1;
    //double merged_pairs[m * 4];    
    int i;   
    int num_procs;
    MPI_Status status;
    int f_row, f_col, c_row, c_col, first_row;


    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    int chunk_size = NUM_VECS / (num_procs - 1); // not including process 0
    double A_chunk[chunk_size * N]; 

    print_oned_mat("\n--------------------------------------------\n\nMatrix A: \n", 
                   A, NUM_VECS, N);
    printf("\n\n");

    printf("Sending address of first row of each chunk\n");
    for (i = 1; i < num_procs; i++) {
        first_row = (i - 1) * chunk_size;

        c_row = 0;
        for (f_row = first_row; f_row < (first_row + chunk_size); f_row++) {
            c_col = 0;
            for (f_col = 0; f_col < N; f_col++) {
                A_chunk[c_row*N + c_col] = *(A + f_row*N + f_col);
                c_col++;
            }
            c_row++;
        }
        // buf, count, datatype, dest, tag, comm
        MPI_Send(A_chunk, chunk_size * N, MPI_DOUBLE, i, chunk_size, MPI_COMM_WORLD);
    }

    // each worker process clusters a chunk of data rows and returns a 
    // matrix of the vectors selected to merge into a group with each iteration within the process
    // so each worker process returns a matrix of size chunk_size * N
    // we want to print each of these matrices    

    //double mp_chunk[chunk_size * 4];
    double A_chunk_root[1*N]; 
    double last_call[num_procs * N];   
    //double last_call[(num_procs-1) * N];
    printf("Receiving merged pairs from workers\n");
    for (i = 1; i < num_procs; i++) {
        
        MPI_Recv(A_chunk_root, 1 * N, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        first_row = (i - 1);
        for (int j = 0; j < N; j++ ) {
            last_call[i * N + j] = A_chunk_root[j];
        }
        /*MPI_Recv(mp_chunk, chunk_size * 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        //print_oned_mat("Worker results: \n", mp_chunk, chunk_size, 4);

        first_row = (i - 1) * chunk_size;
        c_row = 0;
        for (f_row = first_row; f_row < (first_row + chunk_size); f_row++) {
            c_col = 0;
            for (f_col = 0; f_col < 4; f_col++) {
                merged_pairs[f_row*4 + f_col] = mp_chunk[c_row*4 + c_col];
                c_col++;
            }
            c_row++;
        }*/
    }

    sequential_naive(last_call, num_procs);

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    printf("Sending dietags for %d workers\n", num_procs-1); // to get all workers to exit
    for (i = 0; i < num_procs; i++) {
        MPI_Send(0, 0, MPI_DOUBLE, i, 5000, MPI_COMM_WORLD);
    }

    // want to have some kind of iterator to go with the prints (so we can be sure of the order)
    //if (NUM_VECS <= 15) {print_oned_mat("\nResults: \n", merged_pairs, NUM_VECS, 4);}

    
    printf("\n||-----Done parallel part-----||.\n\n");    
}


static void worker(int chunk_size) {
    double A[chunk_size * N];
    int my_rank;
    MPI_Status status;  
    double merged_pairs[chunk_size * 4]; // check not chunk_size+1
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); //grab this process's rank

    for (;;) {
        MPI_Recv(A, chunk_size * N, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == 5000) {
            return;
        }

        //int vec_chunk = chunk_size; // == status.MPI_TAG;
        double Dists[chunk_size * chunk_size]; 
        resetMatrix(chunk_size, chunk_size, Dists);
        double A_new[chunk_size * N];

        printf("\nMatrices for process %d: \n", my_rank);
        if (NUM_VECS <= 15) {print_oned_mat("A_chunk: \n", A, chunk_size, N);}
        
        collectDistances(chunk_size, A, Dists); // collect pairwise distances of all the vectors
        if (NUM_VECS <= 15) {print_oned_mat("\nDists: \n", Dists, chunk_size, chunk_size);}

        double min_val = 1000;
        double coords[3] = {0, 0, min_val};

        minDist(min_val, chunk_size, Dists, coords, chunk_size);

        int val0 = (int) coords[0];
        int val1 = (int) coords[1];

        double new_vec[2];
        int j;
        for (j = 0; j < N; j++) {
            new_vec[j] = (A[val0 * N + j] + A[val1 * N + j]) / 2;
        }
        printf("\n\nMinimum distance %.3f found at (%d, %d), level %d, process %d\n\tFrom [%.2f, %.2f] and [%.2f, %.2f] to [%.2f, %.2f]\n", 
                coords[2], val0, val1, chunk_size, my_rank,
                *(A + val0*N + 0), *(A + val0*N + 1), *(A + val1*N + 0), *(A + val1*N + 1), new_vec[0], new_vec[1]);
        printf("\n--------------------------------------------\n\n");

        //memset(str_full, 0, sizeof(str_full));
        updateActiveMatrix(chunk_size, coords, A, A_new);
        int x_old = chunk_size;
        int x = chunk_size - 1;
        addToMergedList(merged_pairs, A, val0, val1, x);

        copyNewToOld(A_new, A, x_old, x);
        resetMatrix(x, x, Dists);
            
        while (x > 1) {
            if (NUM_VECS <= 15) {print_oned_mat("New A: \n", A, x, N);}
            collectDistances(x, A, Dists);
            //print_oned_mat("\nDists: \n", Dists, x, x);

            min_val = 1000;
            minDist(min_val, x, Dists, coords, chunk_size);

            val0 = (int) coords[0];
            val1 = (int) coords[1];

            for (j = 0; j < N; j++) {
                new_vec[j] = (A[val0 * N + j] + A[val1 * N + j]) / 2;
            }
            printf("\n\nMinimum distance %.3f found at (%d, %d), level %d, process %d\n\tFrom [%.2f, %.2f] and [%.2f, %.2f] to [%.2f, %.2f]\n", 
                coords[2], val0, val1, x, my_rank,
                *(A + val0*N + 0), *(A + val0*N + 1), *(A + val1*N + 0), *(A + val1*N + 1), new_vec[0], new_vec[1]);

            printf("\n--------------------------------------------\n\n");

            resetMatrix(x, x, A_new);
            updateActiveMatrix(x, coords, A, A_new);    

            x_old = x;
            x = x - 1;
            addToMergedList(merged_pairs, A, val0, val1, x);
            copyNewToOld(A_new, A, x_old, x);
            resetMatrix(x, x, Dists);
        }

        print_oned_mat("\n---------New A: \n", A, x, N);

        printf("Returning root of this subtree to manager\n");
        MPI_Send(A, 1 * N, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
        //printf("Returning merged pairs to manager\n");
        //MPI_Send(merged_pairs, chunk_size * 4, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
    }
    
    // the prints are all out of order, want to fix that
    //print_oned_mat("\Worker results: \n", merged_pairs, chunk_size, 4);
    //printf("Returning merged pairs to manager\n");
    //MPI_Send(merged_pairs, chunk_size * 4, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
}



static void sequential_naive(double *A, int num_rows) {    
    int m = num_rows+1;
    double merged_pairs[m * 4];    
    double Dists[num_rows * num_rows]; 
    double A_new[num_rows * N];
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); //grab this process's rank

    print_oned_mat("\n--------------------------------------------\n\nSequential Matrix A: \n", 
                   A, num_rows, N);

    //parallelCollectDistances(num_rows, A, Dists);
    collectDistances(num_rows, A, Dists); // collect pairwise distances of all the vectors
    if (num_rows <= 15) {print_oned_mat("\nDists: \n", Dists, num_rows, num_rows);}

    double min_val = 1000;
    double coords[3] = {0, 0, min_val};

    minDist(min_val, num_rows, Dists, coords, num_rows);

    int val0 = (int) coords[0];
    int val1 = (int) coords[1];

    double new_vec[2];
    int j;
    for (j = 0; j < N; j++) {
        new_vec[j] = (A[val0 * N + j] + A[val1 * N + j]) / 2;
    }

    printf("\n\nMinimum distance %.3f found at (%d, %d), level %d, process %d\n\tFrom [%.2f, %.2f] and [%.2f, %.2f] to [%.2f, %.2f]\n", 
               coords[2], val0, val1, num_rows, my_rank,
               *(A + val0*N + 0), *(A + val0*N + 1), *(A + val1*N + 0), *(A + val1*N + 1), new_vec[0], new_vec[1]);
    printf("\n--------------------------------------------\n\n");

    //memset(str_full, 0, sizeof(str_full));

    updateActiveMatrix(num_rows, coords, A, A_new);
    int x_old = num_rows;
    int x = num_rows - 1;
    addToMergedList(merged_pairs, A, val0, val1, x);

    copyNewToOld(A_new, A, x_old, x);
    resetMatrix(x, x, Dists);

    while (x > 1) {
        if (num_rows <= 15) {print_oned_mat("New A: \n", A, x, N);}
        //parallelCollectDistances(x, A, Dists);
        collectDistances(x, A, Dists);
        //print_oned_mat("\nDists: \n", Dists, x, x);

        min_val = 1000;
        minDist(min_val, x, Dists, coords, num_rows);

        val0 = (int) coords[0];
        val1 = (int) coords[1];

        for (j = 0; j < N; j++) {
            new_vec[j] = (A[val0 * N + j] + A[val1 * N + j]) / 2;
        }
        printf("\n\nMinimum distance %.3f found at (%d, %d), level %d, process %d\n\tFrom [%.2f, %.2f] and [%.2f, %.2f] to [%.2f, %.2f]\n", 
               coords[2], val0, val1, x, my_rank,
               *(A + val0*N + 0), *(A + val0*N + 1), *(A + val1*N + 0), *(A + val1*N + 1), new_vec[0], new_vec[1]);

        printf("\n--------------------------------------------\n\n");

        resetMatrix(x, x, A_new);
        updateActiveMatrix(x, coords, A, A_new);    

        x_old = x;
        x = x - 1;
        addToMergedList(merged_pairs, A, val0, val1, x);
        copyNewToOld(A_new, A, x_old, x);
        resetMatrix(x, x, Dists);
    }

    // the prints are all out of order, want to fix that
    if (num_rows <= 15) {print_oned_mat("\nSequential results: \n", merged_pairs, m, 4);}

    printf("\n||-----Done sequential-----||.\n\n");
}






