#ifndef COMMON_H
#define COMMON_H
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>

#define NEW(x,n,t) if ((x=(t*) malloc((n)*sizeof(t)))==NULL) {\
                      fprintf(stderr,"NEW:Out of Memory.Require=%d*%d\n",\
                      (n),sizeof(t)); exit(1);}

#define NEWP(x,n,t) if ((x=(t**) malloc((n)*sizeof(t*)))==0) { \
                       fprintf(stderr,"NEWP:Out of Memory.Require=%d*%d",\
                       (n),sizeof(t*));exit(1);}

#define FREEP(x, n) 	{for(int i=0;i < (n);i++) \
                           if((x)[i] != NULL) free((x)[i]);\
                        if((x) != NULL) free((x));(x) = NULL;}

#define FREEPP(x,ii,jj)	{for(int i=0;i <(ii);i++){\
			   for(int j=0;j<(jj);j++) \
                             free((x)[i][j]);\
                           free((x)[i]);}\
			(x)=NULL;}



#define min(A,B)  ((A) < (B) ? (A) : (B))
#define max(A,B)  ((A) > (B) ? (A) : (B))


bool MC(float s1, float s2, float T);
int min_index(float* s, int size);

void initrandom();
float randomf(float maxval);
int randomi(int maxval);
int round(float);

float           parabolic(float xa, float xb, float xc, 
			  float fa, float fb, float fc);


string itos (int i);

#endif 
