#ifndef functions
#define functions
#include <stdio.h>

/* Constants */

#define N 2
#define PRINT_VECS 1  // flag so we can turn off printing when N is large
#define MAX_RAND 10   // max value of elements generated for array
#define NUM_VECS 2000
//#define NTHREADS 8

/* Prototypes */
double pow(double x, double exp);
double floor(double x);
double sqrt(double x);
void init_vec(double *vec, int len);
void print_vec(const char *label, double *vec, int len);
void print_matrix(const char *label, int dim1, int dim2, double (*matrix)[dim2]);
void print_oned_mat(const char *label, double *A, int dim1, int dim2);
void minDist(double min_val, int x, double *Dists, double *coords, int dim2); 
void updateActiveMatrix(int num_vecs, double *coords, double *A, double *A_new);
//void seqCollectDistances(int y, double *A, double *Dists);
void collectDistances(int y, double *A, double *Dists);
void copyNewToOld(double *A_new, double *A, int x_old, int x_new);
void resetMatrix(int dim1, int dim2, double *matrix); 
double concatenate(double x, double y);
void addToMergedList(double *merged_pairs, double *A, int val0, int val1, int x); 
void concatPair(char *str0, double val0, double val1);
void concatStrings(char *str0, char *str1, char *str2);
double min(double x, double y);
double max(double x, double y);
double maxCoord(int dim1, int dim2, int j, double *matrix);
double minCoord(int dim1, int dim2, int j, double *matrix);
//memset(str_full, 0, sizeof(str_full));

#endif
