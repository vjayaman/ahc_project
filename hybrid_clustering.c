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
    
    // Initialize matrix A of random values, dimension NUM_VECS X N
    srand(200);
    for (i = 0; i < NUM_VECS; i++) {
        for (j = 0; j < N; j++) {
            init_vec(x_i, N);
            A[i * N + j] = *x_i;
        }
    }
    
    int my_rank;
    int num_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    double start_time = 0;
    double end_time = 0;   
    
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    // if myjob has more than one process specified, run the parallel 
    // implementation, otherwise, the sequential
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
    return 0;    
}



static void manager(double *A) {
    int num_procs, f_row, f_col, c_row, c_col;
    
    MPI_Status status;  
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int chunk_size = NUM_VECS / (num_procs - 1); // not including process 0
    double A_chunk[chunk_size * N]; 

    // Print the initial matrix of vectors for clustering
    print_oned_mat("\n--------------------------------------------\n\nMatrix A: \n", 
                   A, NUM_VECS, N);
    printf("\n\n");

    // Identify the bounds of the matrix (both in x and y directions)
    double minX = minCoord(NUM_VECS, N, 0, A);
    double maxX = maxCoord(NUM_VECS, N, 0, A);
    double minY = minCoord(NUM_VECS, N, 1, A);
    double maxY = maxCoord(NUM_VECS, N, 1, A);
    printf("Min x: %f, Max x: %f, Min y: %f, Max y: %f\n", minX, maxX, minY, maxY);

    // Simple sort of the vectors of A by their y value
    int i,j, k;
    double temp[N]; 
    for (i = 0; i < NUM_VECS; i++){
        for (k = i+1; k < NUM_VECS; k++) {
            if (A[i * N + 1] > A[k * N + 1]) {  // sorting by y values
                for (j = 0; j < N; j++) {
                    temp[j] = A[i * N + j];
                    A[i * N + j] = A[k * N + j];
                    A[k * N + j] = temp[j];
                }
            }
        }
    }
    
    // This was when the grid mapping was considered
    minX = minX - 1; // so we can always use (< x <=), even for boundaries
    
    // Copy a chunk of the matrix and send to any worker process
    int first_row;
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

        //print_oned_mat("\nA_chunk\n", A_chunk, chunk_size, N);
        //printf("\n");   

        // buf, count, datatype, dest, tag, comm
        MPI_Send(A_chunk, chunk_size * N, MPI_DOUBLE, i, chunk_size, MPI_COMM_WORLD);
    }
    
    // each worker process clusters a chunk of data rows and returns the root 
    // of the subtree for that chunk
    // manager then clusters these subtrees
    
    double A_chunk_root[1*N]; 
    double last_call[num_procs * N];   
    
    // manager receives a vector (subtree root) from each worker and 
    // saves into a matrix for sequential clustering
    for (i = 1; i < num_procs; i++) {
        MPI_Recv(A_chunk_root, 1 * N, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        for (int j = 0; j < N; j++ ) {
            last_call[i * N + j] = A_chunk_root[j];
        }
    }
    print_oned_mat("\n--------------------------------------------\n\nLastcall Matrix A: \n", 
                    last_call, num_procs-1, N);     
    sequential_naive(last_call, num_procs-1);

    // Confirming the number of workers
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // to get all workers to exit
    printf("\nSending dietags for %d workers\n", num_procs-1); 
    for (i = 0; i < num_procs; i++) {
        MPI_Send(0, 0, MPI_DOUBLE, i, 5000, MPI_COMM_WORLD);
    }

    printf("\n||-----Done parallel part-----||.\n\n");    
}


static void worker(int chunk_size) {
    double A[chunk_size * N];
    int my_rank, j;
    MPI_Status status;  
    double merged_pairs[chunk_size * 4]; // check not chunk_size+1
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); //grab this process's rank

    for (;;) {
        // worker receives a chunk of vectors
        MPI_Recv(A, chunk_size * N, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == 5000) {return;}

        // if the number of workers does not equal the number of vectors:
        if (chunk_size > 1) {
            // garbage removal
            double Dists[chunk_size * chunk_size]; 
            resetMatrix(chunk_size, chunk_size, Dists);
            double A_new[chunk_size * N];
            
            // collect pairwise distances of all the vectors
            collectDistances(chunk_size, A, Dists); 
            //print_oned_mat("\nDists: \n", Dists, chunk_size, chunk_size);

            // arbitrary minimum value; coords will contain indices of minimum 
            // value in the Dists matrix as well as the minimum value itself
            double min_val = 1000;
            double coords[3] = {0, 0, min_val};
            minDist(min_val, chunk_size, Dists, coords, chunk_size);

            int val0 = (int) coords[0];
            int val1 = (int) coords[1];

            // identify the new center of this merged pair and print for result processing
            double new_vec[2];
            for (j = 0; j < N; j++) {
                new_vec[j] = (A[val0 * N + j] + A[val1 * N + j]) / 2;
            }

            printf("\n\nMinimum distance %.3f found at (%d, %d), level %d, process %d\tFrom [%.2f, %.2f] and [%.2f, %.2f] to [%.2f, %.2f]\n", 
                coords[2], val0, val1, chunk_size, my_rank,
                *(A + val0*N + 0), *(A + val0*N + 1), *(A + val1*N + 0), *(A + val1*N + 1), new_vec[0], new_vec[1]);
            printf("\n--------------------------------------------\n\n");            

            // update the matrix of active vectors for clustering 
            // i.e. add the new center and remove its elements from A
            // x is the level of clustering (starts large and decrements to 2 before we return a result)
            // we have a list of merged pairs as well, and reset the Dists matrix (so there is no leftover data)
            updateActiveMatrix(chunk_size, coords, A, A_new);
            int x_old = chunk_size;
            int x = chunk_size - 1;
            addToMergedList(merged_pairs, A, val0, val1, x);
            copyNewToOld(A_new, A, x_old, x);
            resetMatrix(x, x, Dists);
            
            while (x > 1) {
                //if (NUM_VECS <= 15) {print_oned_mat("New A: \n", A, x, N);}
                collectDistances(x, A, Dists);
                print_oned_mat("\nDists: \n", Dists, x, x);

                min_val = 1000;
                minDist(min_val, x, Dists, coords, chunk_size);

                val0 = (int) coords[0];
                val1 = (int) coords[1];

                for (j = 0; j < N; j++) {
                    new_vec[j] = (A[val0 * N + j] + A[val1 * N + j]) / 2;
                }
                printf("\n\nMinimum distance %.3f found at (%d, %d), level %d, process %d\tFrom [%.2f, %.2f] and [%.2f, %.2f] to [%.2f, %.2f]\n", 
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
        }           

        //printf("\nReturning root of this subtree to manager\n");
        MPI_Send(A, 1 * N, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
    }
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

    collectDistances(num_rows, A, Dists); // collect pairwise distances of all the vectors
    print_oned_mat("\nDists: \n", Dists, num_rows, num_rows);

    double min_val = 1000;
    double coords[3] = {0, 0, min_val};

    minDist(min_val, num_rows, Dists, coords, num_rows);

    int val0 = (int) coords[0];
    int val1 = (int) coords[1];
    int j;
    double new_vec[2];
    for (j = 0; j < N; j++) {
        new_vec[j] = (A[val0 * N + j] + A[val1 * N + j]) / 2;
    }

    printf("\n\nMinimum distance %.3f found at (%d, %d), level %d, process %d\tFrom [%.2f, %.2f] and [%.2f, %.2f] to [%.2f, %.2f]\n", 
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
        //if (num_rows <= 15) {print_oned_mat("New A: \n", A, x, N);}
        collectDistances(x, A, Dists);
        //print_oned_mat("\nDists: \n", Dists, x, x);

        min_val = 1000;
        minDist(min_val, x, Dists, coords, num_rows);

        val0 = (int) coords[0];
        val1 = (int) coords[1];

        for (j = 0; j < N; j++) {
            new_vec[j] = (A[val0 * N + j] + A[val1 * N + j]) / 2;
        }
        printf("\n\nMinimum distance %.3f found at (%d, %d), level %d, process %d\tFrom [%.2f, %.2f] and [%.2f, %.2f] to [%.2f, %.2f]\n", 
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







