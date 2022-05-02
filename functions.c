#include "functions.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Functions */

// if (x > y) return y else x
double min(double x, double y) {
    return (x > y) ? y : x;
}

// if (x > y) return x else y
double max(double x, double y) {
    return (x > y) ? x : y;
}

double fabs(double x);
double pow(double x, double exp);
double sqrt(double x);

// Fills a vector with random integers in the range [0, MAX_RAND)
// From vec_sum example provided by prof
void init_vec(double *vec, int len)
{
    int i;
    for (i = 0; i < len; i++)
    {
        vec[i] = rand() % MAX_RAND;
    }    
}

// Prints the given vector to stdout
// From vec_sum example provided by prof
void print_vec(const char *label, double *vec, int len)
{
#if PRINT_VECS
    printf("%s", label);
    
    int i;
    for (i = 0; i < len; i++)
    {
        printf("%.2f ", vec[i]);
    }
    //printf("\n\n");
#endif
}

// Prints the given 2D matrix to stdout
// From vec_sum example provided by prof
void print_matrix(const char *label, int dim1, int dim2, double (*matrix)[dim2])
{
    printf("%s", label);
    
    int i, j;
    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim2; j++) {
            printf("%.3f \t", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    //printf("\n\n");
}

// reset matrix values to 0 (make sure no garbage values are retained)
// useful for initialization as well
void resetMatrix(int dim1, int dim2, double *matrix) {
    int i, j;
    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim2; j++) {
            matrix[i * dim2 + j] = 0;
        }
    }
}

// identify the maximum value along one column of the A matrix (for bounds)
double maxCoord(int dim1, int dim2, int j, double *matrix) {
    int i;
    double maxval = 0;
    for (i = 0; i < dim1; i++) {
        maxval = max(maxval, matrix[i * dim2 + j]);
    }
    //printf("min %d: %f, max %d: %f\n", j, minval, j, maxval);
    return maxval;
}

// identify the minimum value along one column of the A matrix (for bounds)
double minCoord(int dim1, int dim2, int j, double *matrix) {
    int i;
    double minval = 5000;
    for (i = 0; i < dim1; i++) {
        minval = min(minval, matrix[i * dim2 + j]);
    }
    return minval;
}

// Prints the given 2D contiguous array to stdout
void print_oned_mat(const char *label, double *A, int dim1, int dim2) {
    int i,j;
    printf("%s", label);
    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim2; j++) {
            printf("%.2f\t", *(A + i*dim2 + j));
        }
        printf("\n");
    }
}

// collect the pairwise distances for vectors in A
// save into the upper triangular matrix of Dists (not including the diagonal)
void collectDistances(int y, double *A, double *Dists) {
    double sum;
    int row0, row1, i2;
    for (row0 = 0; row0 < y; row0++) {
        for (row1 = row0+1; row1 < y; row1++) {
            sum = 0;
            //#pragma omp parallel for private(i2) reduction(+:sum)
            for (i2 = 0; i2 < N; i2++) {
                sum += pow(A[row0 * N + i2] - A[row1 * N + i2], 2);
            }
            Dists[row0 * y + row1] = sqrt(sum);
        }
    }
}


// identify the minimum pairwise distance in matrix Dists and save the 
// indices for where it was found
// the distance must be greater that 0.005
void minDist(double min_val, int x, double *Dists, double *coords, int dim2) {
    int i1, i2;
    double mv = min_val;
    for (i1 = 0; i1 < x; i1++) {
        for (i2 = i1+1; i2 < x; i2++) {
            if (Dists[i1 * x + i2] > 0.005) {
                if ((Dists[i1 * x + i2] < mv)) {
                    mv = Dists[i1 * x + i2];
                    coords[0] = i1;
                    coords[1] = i2;
                    //printf("coords0: %f, coords1: %f, min val: %f\n", coords[0], coords[1], mv);
                    coords[2] = mv;
                }
            }
        }
    }
}

// update the A matrix: add the new center vector and remove the old elements 
// that are now grouped under the center vector
void updateActiveMatrix(int num_vecs, double *coords, double *A, double *A_new) { // double *indices
    int k1 = 0;
    int i, j;
    for (i = 0; i < num_vecs; i++) {
        if ((i != (int)coords[0]) && (i != (int)coords[1])) {
            for (j = 0; j < N; j++) {
                A_new[k1 * N + j] = A[i * N + j];
            }
            k1++;
        }
    }

    int val0 = (int)coords[0];
    int val1 = (int)coords[1];
    for (j = 0; j < N; j++) {
        A_new[k1 * N + j] = (A[val0 * N + j] + A[val1 * N + j]) / 2;
        //printf("(A[%d][j] + A[%d][j]) / 2 = (%.2f + %.2f)/2 = %.2f\n", 
        //       (int)coords[0], (int)coords[1], A[(int)coords[0]][j], A[(int)coords[1]][j], A_new[k1][j]);
    }
}

// copy 2D contiguous array values from one to another
void copyNewToOld(double *A_new, double *A, int x_old, int x_new) {
    int i, j;
    for (i = 0; i < x_new; i++) {
        for (j = 0; j < N; j++) {
            A[i * N + j] = A_new[i * N + j];
        }
    }
}


// https://stackoverflow.com/questions/12700497/how-to-concatenate-two-integers-in-c
double concatenate(double x, double y) {
    double pow = 10;
    while(y >= pow) {pow *= 10;}
    return x * pow + y;
}

// list of merged pairs, adding to a list of these for result processing
void addToMergedList(double *merged_pairs, double *A, int val0, int val1, int x) {
    int i;
    for (i = 0; i < 2; i++) {
        merged_pairs[(x * N*2) + i] = A[val0 * N + i];
        merged_pairs[(x * N*2) + i+2] = A[val1* N + i];
    }
}


void concatStrings(char *str0, char *str1, char *str2) {
    strcat(str0, str1);
    strcat(str0, "-");
    strcat(str0, str2);
    strcat(str0, ")");
}


void concatPair(char *str0, double val0, double val1) {

    char str1[12];
    char str2[12];
    sprintf(str1, "%d", (int) val0);
    sprintf(str2, "%d", (int) val1);

    strcat(str0, str1);
    strcat(str0, "-");
    strcat(str0, str2);
    strcat(str0, ")");
}
