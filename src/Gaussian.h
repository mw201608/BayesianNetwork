#ifndef GAUSSIAN_H
#define GAUSSIAN_H

//local definitions
#include "Matrix.h"
#include "common.h"

class Gaussian {
 public: 
  Gaussian(){}

 public: 
  static float bic(int , float*, float alpha);  //marginal 
  static float bic(int dim, float* x, int np, float** x2, float alpha ); //conditional
};
#endif
