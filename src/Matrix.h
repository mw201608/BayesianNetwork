#ifndef MATRIX_H
#define MATRIX_H
#include "nrutil.h"
#define TINY 1.0e-20

class Matrix {
public:
  Matrix(){}
public:
  static void   matrixmult(float** m1, int dx1, int dy1,
			  float** m2, int dx2, int dy2,
			  float** m) ;
  static void   covariance(float** m1, int dx1, int dy1,float** m, float* avg);
  static void   transpose(float** m1, int dx, int dy, float** m);
  static float  mean(float* m,int dy);
  static void   mean(float**m, int dx, int dy, float* avg);
  static float  std(float* m,int dy);
  static void   detAndInverse(float** m, int dx, double* det, float** mi);
  static void   lubksb(float **a, int n, int *indx, float* b);
  static void   ludcmp(float **a, int n, int *indx, float *d);
  static void   svdcmp(float **a, int m, int n, float *w, float **v);
  static float  pythag(float a, float b);
  static void   svbksb(float **u, float* w, float **v, int m, int n, float b[], float x[]);
  static void   print(float **m, int dx, int dy);
  static void   subtract(float* x, int dim, float mu, float* dx);
  static void   subtract(float* x, int dim, float* mu12, float* dx);
  static void   subtract(float**, int dx, int dy, float*, float**);
  static void   add(float**, int dx, int dy, float, float**);
  static float  dot(float* x, int dim, float* y);

  static void   svdfit(float** x, float y[], int ndata, float a[], int ma);
};


#endif
