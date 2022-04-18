#include "functions.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

static void sequential(double *A) {    
    int m = NUM_VECS+1;
    double merged_pairs[m * 4];    
    double Dists[NUM_VECS * NUM_VECS]; 
    double A_new[NUM_VECS * N];
    
    print_oned_mat("\n--------------------------------------------\n\nMatrix A: \n", 
                       A, NUM_VECS, N);

    parallelCollectDistances(NUM_VECS, A, Dists);
    //collectDistances(NUM_VECS, A, Dists); // collect pairwise distances of all the vectors
    print_oned_mat("\nDists: \n", Dists, NUM_VECS, NUM_VECS);

    double min_val = 1000;
    double coords[3] = {0, 0, min_val};

    minDist(min_val, NUM_VECS, Dists, coords, NUM_VECS);

    int val0 = (int) coords[0];
    int val1 = (int) coords[1];

    double new_vec[2];
    int j;
    for (j = 0; j < N; j++) {
        new_vec[j] = (A[val0 * N + j] + A[val1 * N + j]) / 2;
    }
    printf("\n\nMinimum distance %.3f found at (%d, %d)\n\tFrom [%.2f, %.2f] and [%.2f, %.2f] to [%.2f, %.2f]\n", 
               coords[2], val0, val1, 
               *(A + val0*N + 0), *(A + val0*N + 1), *(A + val1*N + 0), *(A + val1*N + 1), new_vec[0], new_vec[1]);
    printf("\n--------------------------------------------\n\n");

    //memset(str_full, 0, sizeof(str_full));

    updateActiveMatrix(NUM_VECS, coords, A, A_new);
    int x_old = NUM_VECS;
    int x = NUM_VECS - 1;
    addToMergedList(merged_pairs, A, val0, val1, x);

    copyNewToOld(A_new, A, x_old, x);
    resetMatrix(x, x, Dists);

    while (x > 1) {
        print_oned_mat("New A: \n", A, x, N);
        parallelCollectDistances(x, A, Dists);
        //collectDistances(x, A, Dists);
        //print_oned_mat("\nDists: \n", Dists, x, x);

        min_val = 1000;
        minDist(min_val, x, Dists, coords, NUM_VECS);

        val0 = (int) coords[0];
        val1 = (int) coords[1];

        for (j = 0; j < N; j++) {
            new_vec[j] = (A[val0 * N + j] + A[val1 * N + j]) / 2;
        }
        printf("\n\nMinimum distance %.3f found at (%d, %d)\n\tFrom [%.2f, %.2f] and [%.2f, %.2f] to [%.2f, %.2f]\n", 
               coords[2], val0, val1, 
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
    print_oned_mat("\nSequential results: \n", merged_pairs, m, 4);
    printf("\n||-----Done sequential-----||.\n\n");
}
