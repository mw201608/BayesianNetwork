/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A package for Bayesian network reconstruction
%
% Author:  Dr. Jun Zhu, Amgen Inc., Thousand Oaks, CA, 2002
%          Dr. Jun Zhu, Rosetta Informatics, a wholly owned subsidiary of Merck & CO., Seattle, WA, 2004, 2008
%
% Acknowledge:  Thanks to Amgen Inc. and Merck & CO. for their generous supports
%
% If you use this package, please cite the following references
%
% (1) Zhu J, Lum PY, Lamb J, GuhaThakurta D, Edwards SW, et al. An integrative genomics approach to the
%     reconstruction of gene networks in segregating populations. Cytogenet Genome Res 105: 363-374 (2004)
% (2) Zhu J, Wiener MC, Zhang C, Fridman A, Minch E, Lum PY, Sachs JR, & Schadt EE Increasing the power to
%     detect causal associations by combining genotypic and expression data in segregating populations
%     PLoS Comput Biol 3, e69. (2007)
% (3) Zhu, J., Zhang, B., Smith, E.N., Drees, B., Brem, R.B., Kruglyak, L., Bumgarner, R.E. and Schadt, E.E.
%     Integrating large-scale functional genomic data to dissect the complexity of yeast regulatory networks,
%     Nat Genet, 40, 854-861 (2008)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <math.h>
#include "Gaussian.h"
#include "common.h"

/** A Gaussian (normal) distribution.
  */

/* conditional Gaussian 
 * The marginal means of the child and parent variables,
  * respectively, as mu(1) and mu(2), and the respective
  * marginal variances as Sigma(11) and Sigma(22), and the
  * covariance as <tt>Sigma(12)</tt>, 
  * then the conditional mean mu(1|2)
  * and conditional variance Sigma(1|2) are as follows.
  * 
  *     mu(1|2) = mu(1) + Sigma(12) Sigma(22)^{-1} (X(2)-mu(2))
  *     Sigma(1|2) = Sigma(11) - Sigma(12) Sigma(22)^{-1} Sigma(21)
  *
  *     Sigma(12)Sigma(22)^{-1} is the regression coefficient which will
  *     be calculated svdfit 
  *     
  */
float Gaussian::bic(int dim, float* x, int np, float** x2, float alpha ) 
{

    // parameter need to estimate maximum likelihood 
    // (mean, variance and covariance)
    int npara = (np+1)*np/2+2;
    float* mu;
    NEW(mu, (1+np), float);
    float** sigma= matrix(0,np,0,np);
    float** data = matrix(0, dim-1, 0, np);
    for (int i=0; i<dim; i++) {
      data[i][0] = x[i];
    }
    for (int j=0; j<dim; j++) {
      for (int i=1; i<=np; i++) {
	data[j][i] =x2[i-1][j];
      }
    }

    //calculate covariance matrix
    Matrix::covariance(data, dim, np+1, sigma, mu);
    float mu2[np];
    for (int i=0; i<np; i++) {
      mu2[i] = mu[i+1];
    }
    float** sigma12 = matrix(0, 0, 0, np-1);
    for (int i=0; i<np; i++) {
      sigma12[0][i] = sigma[0][i+1];
    }
    /*
    float** sigma22 = matrix(0, np-1, 0, np-1);
    for (int i=0; i<np; i++) {
      for (int j=0; j<np; j++) {
	sigma22[i][j] = sigma[i+1][j+1];
      }
    }
 
    double sigma22_det;
    float** sigma22_inverse = matrix(0, np-1, 0, np-1);
    Matrix::detAndInverse(sigma22, np, &sigma22_det, sigma22_inverse);
    */

    /*
     *calculate conditional mean
     *mu(1|2) = mu(1) + Sigma(12) Sigma(22)^{-1} (X(2)-mu(2))
     */
    float** mu12 = matrix(0,0, 0,  dim-1);
    float** t = matrix(0, 0, 0, np-1);    
    //Matrix::matrixmult(sigma12, 1, np, sigma22_inverse, np, np, t);
    float*  ta; NEW(ta, np, float);
    float** X2 = matrix(0, np-1, 0, dim-1);
    Matrix::subtract(x2, np, dim, mu2, X2);
    float* X; NEW(X, dim, float); 
    Matrix::subtract(x, dim,  mu[0], X);
    Matrix::svdfit(X2, X, dim, ta, np);
    
    for (int i=0; i<np; i++) {
      t[0][i] =ta[i];
    }
    free(ta);
   
    Matrix::matrixmult(t, 1, np, X2, np, dim, mu12);
    Matrix::add(mu12, 1, dim, mu[0], mu12);
    
    /*
     * conditional variance Sigma(1|2) are as follows.
     * Sigma(1|2) = Sigma(11) - Sigma(12) Sigma(22)^{-1} Sigma(21)
     */
    float** sigma21 = matrix(0, np-1, 0, 0);
    Matrix::transpose(sigma12, 1, np, sigma21);
    float** st = matrix(0,0,0,0);
    Matrix::matrixmult(t, 1, np, sigma21, np, 1, st);
    float Sigma12 = sigma[0][0]-st[0][0];

    //
    float* dx;
    NEW(dx, dim, float);
    Matrix::subtract(x, dim, mu12[0], dx);
    float pp = -Matrix::dot(dx, dim, dx)/(2*Sigma12)
               -dim*log(sqrt(Sigma12))- dim*0.5*log(2*3.1415);
    
    //bic score
    //float b = pp-0.5*npara*log(dim);
    float b = pp-alpha*npara;

    //free space here
    free(mu);
    free_matrix(sigma, 0,np,0,np);
    free_matrix(data, 0, dim-1, 0, np);
    free_matrix(sigma12, 0, 0, 0, np-1);
    free_matrix(mu12, 0,0, 0,  dim-1);
    free_matrix(t, 0, 0, 0, np-1);  
    free_matrix(X2, 0, np-1, 0, dim-1);
    free(X);
    free_matrix(sigma21, 0, np-1, 0, 0);
    free_matrix(st, 0, 0, 0, 0);
    free(dx);

    return b;
}

/** Computes the density of a 1-dimensional Gaussian with the given
  * mean and standard deviation (not the variance).
  */
float Gaussian::bic(int dim, float *x, float alpha)
{
    // parameter need to estimate maximum likelihood (mean and std)
    int npara = 2;
    float mu = Matrix::mean(x, dim);
    float std = Matrix::std(x, dim);
    float* v;
    NEW(v, dim, float);
    //
    float dot = 0;
    for (int i=0; i<dim; i++) {
      v[i] = x[i]-mu;
      dot += v[i]*v[i];
    }

    //ML 
    //float pp = -0.5*Matrix::dot(v, dim, v)/(std*std)-dim*log(std)- dim*0.5*log(2*3.1415);
    float pp = -0.5*dot/(std*std)-dim*log(std)- dim*0.5*log(2*3.1415);
    //bic score
    //float b = pp-0.5*npara*log(dim);
    float b = pp-alpha*npara;

    //free space
    free( v);

    return b;
}
