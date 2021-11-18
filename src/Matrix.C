#include <stdio.h>
#include <fstream>
#include <math.h>

#include "Matrix.h"
/*
 * a collection of matrix operation
 */


void Matrix::matrixmult(float** m1, int dx1, int dy1,
		       float** m2, int dx2, int dy2,
		       float** m) 
{
  /* initializing */
  for(int x=0; x<dx1; x++) {
    for (int y=0; y<dy2; y++) {
      m[x][y]=0;
    }
  }

  /* calculating */
  for(int x=0; x<dx1; x++) {
    for(int y=0; y<dy2; y++) {
      for(int i=0; i<dy1; i++) {
	m[x][y]+=m1[x][i]*m2[i][y];
      }
    }
  }
}

void Matrix::covariance(float** m1, int dx, int dy,float** m, float* avg)
{
  //calculate average
  for(int y=0; y<dy; y++) {
    avg[y]=0;
  }
  for (int x=0; x<dx; x++) {    
    for(int y=0; y<dy; y++) {
      avg[y]+=m1[x][y];
    }
  }
  for(int y=0; y<dy; y++) {
    avg[y]/=dx;
  }
      
  float** m2= matrix(0,dx-1,0,dy-1);
  for(int x=0; x<dx; x++) {
    for(int y=0; y<dy; y++) {
      m2[x][y]=m1[x][y]-avg[y];
    }
  }
  float** mt= matrix(0,dy-1,0,dx-1);
  transpose(m2,dx,dy,mt);

  matrixmult(mt,dy,dx,m2,dx,dy,m);
  for(int i=0;i<dy; i++) {
    for(int j=0;j<dy; j++) {
      m[i][j]/=dx;
    }
  }

  free_matrix(m2,0,dx-1,0,dy-1);
  free_matrix(mt,0,dy-1,0,dx-1);
}

void Matrix::transpose(float** m1, int dx, int dy, float** m)
{
  for(int x=0; x<dx; x++) {
    for(int y=0; y<dy; y++) {
      m[y][x]=m1[x][y];
    }
  }
}

float Matrix::mean(float* m,int dy) {
  float mf=0;

  for(int i=0; i<dy; i++) mf+=m[i];
  mf=mf/dy;
  return mf;
}

float Matrix::std(float* m,int dy) {
  float mstd=0;
  float avg=mean(m,dy);

  for(int i=0; i<dy; i++) {
    mstd+=(m[i]-avg)*(m[i]-avg);
  }
  mstd=sqrt(mstd/(dy-1));
  return mstd;
}

void Matrix::detAndInverse(float** m, int dx, double* det, float** mi) {
  int* indx= new int[dx];
  float* col= new float[dx];

  float  sign;
  ludcmp(m,dx,indx,&sign);
  //calculate determinant
  *det=sign;
  for(int i=0; i<dx; i++) *det*=m[i][i];

  //calculate inverse
  for(int j=0; j<dx; j++) {
    for (int i=0; i<dx; i++) col[i]=0;
    col[j]=1.0;
    lubksb(m,dx,indx,col);
    for(int i=0; i<dx; i++) mi[i][j]=col[i];
  }

  delete col;
  delete indx;
}

void Matrix::lubksb(float **a, int n, int *indx, float b[])
{
  int i,ip,j;
  float sum;
  
  for (i=0;i<n;i++) {//lower part
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    for (j=0;j<=i-1;j++) sum -= a[i][j]*b[j];
    b[i]=sum;
  }
  
  for (i=n-1;i>=0;i--) { //upper part
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


void Matrix::ludcmp(float **a, int n, int *indx, float *d)
{
  int i,imax,j,k;
  float big,dum,sum,temp;
  float *vv;    //vv stores the implicit scaling of each row

  vv=Vector(0,n-1);
  *d=1.0;
  //no row interchange yet. Loop over rows to get the implicit scaling info
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;       // save the scaling
  }
  //loop over columns of Crout's method
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;  //initialize for the search for largest pivot element
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {//interchange rows
      for (k=0;k<n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n-1) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  free_vector(vv,0,n-1);
}

void Matrix::print(float** m, int dx, int dy) {
  for(int x=0; x <dx; x++) {
    for(int y=0; y<dy; y++) {
      printf("%-9.4g  ",m[x][y]);
    }
    printf("\n");
  }
}

void  Matrix::subtract(float* x, int dim, float y, float* res) {
  for (int i=0; i<dim; i++) {
    res[i] = x[i]- y;
  }
}

void  Matrix::subtract(float* x, int dim, float* y, float* res) {
  for (int i=0; i<dim; i++) {
    res[i] = x[i]- y[i];
  }
}

void Matrix::subtract(float** x, int dx, int dy, float* y, float** res) {
  for (int i=0; i<dx; i++) {
    for (int j=0; j<dy; j++) {
      res[i][j] = x[i][j] - y[i];
    }
  }
}

void Matrix::add(float** x, int dx, int dy, float y, float** res) {
  for (int i=0; i<dx; i++) {
    for (int j=0; j<dy; j++) {
      res[i][j] = x[i][j] - y;
    }
  } 
}
  
float Matrix::dot(float* x, int dim, float* y) {
  float res= 0;
  for (int i=0; i<dim; i++) {
    res += x[i]*y[i];
  }

  return res;
}

/*
 *general linear fit by svdcmp
 * y=ax
 * ma: dimension of a
 * ndata: number of data points for fitting
 * a:  coeffients for return
 *
 */
void Matrix::svdfit(float** x, float y[], int ndata, float a[], int ma)
{
  float TOL = 1.0e-5;
  int j,i;
  float wmax,thresh;

  float** u = matrix(0, ndata-1, 0, ma-1);
  Matrix::transpose(x, ma, ndata, u);
  float* w = Vector(0, ma-1);
  float** v = matrix(0,ma-1, 0, ma-1);
  
  Matrix::svdcmp(u,ndata,ma,w,v);
  //set weight
  wmax = 0;
  for (j=0;j<ma;j++)
    if (w[j] > wmax) wmax=w[j];
  thresh=TOL*wmax;
  for (j=0;j<ma;j++)
    if (w[j] < thresh) w[j]=0.0;
  
  //set coefficents
  Matrix::svbksb(u,w,v,ndata,ma,y,a);

  /*
  cout<<"fitted "<<ndata<<" points\n";
  float sum;
  for (i=0;i<ndata;i++) {
    for (sum=0.0,j=0;j<ma;j++) sum += a[j]*x[j][i];
    cout<<y[i]<<" "<<x[0][i]<<" "<<sum<<endl;
  }
  */

  //free space
  free_matrix(u, 0, ndata-1, 0, ma-1);
  free_vector(w, 0, ma-1);
  free_matrix(v, 0, ma-1, 0, ma-1);
}

float Matrix::pythag(float a, float b)
{
        float absa,absb;
        absa=fabs(a);
        absb=fabs(b);
        if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
        else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void Matrix::svdcmp(float **a, int m, int n, float w[], float **v)
{
  int flag,i,its,j,jj,k,l,nm;
  float anorm,c,f,g,h,s,scale,x,y,z,*rv1;

  rv1=Vector(0,n-1);
  g=scale=anorm=0.0;
  for (i=0;i<n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a[k][i]);
      if (scale) {
	for (k=i;k<m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i < m && i != n-1) {
      for (k=l;k<n;k++) scale += fabs(a[i][k]);
      if (scale) {
	for (k=l;k<n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<n;k++) rv1[k]=a[i][k]/h;
	for (j=l;j<m;j++) {
	  for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l;k<n;k++) a[i][k] *= scale;
      }
    }
    anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    if (i < n) {
      if (g) {
	for (j=l;j<n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=IMIN(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<n;j++) {
	for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<m;j++) a[j][i] *= g;
    } else for (j=i;j<m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=0;l--) {
	nm=l-1;
	if ((float)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((float)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((float)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=Matrix::pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=Matrix::pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=Matrix::pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=Matrix::pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=0;jj<m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free_vector(rv1,0,n-1);
}


void Matrix::svbksb(float **u, float w[], float **v, int m, int n, float b[], float x[])
{
  int jj,j,i;
  float s,*tmp;

  tmp=Vector(0,n-1);
  for (j=0;j<n;j++) {
    s=0.0;
    if (w[j]) {
      for (i=0;i<m;i++) s += u[i][j]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (j=0;j<n;j++) {
    s=0.0;
    for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
    x[j]=s;
  }
  free_vector(tmp,0,n-1);
}
