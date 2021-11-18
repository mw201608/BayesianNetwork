#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include "common.h"

void initrandom() {
   int t = time(NULL);
   printf("init seed : %d\n",t);
   srand(t);
}

float randomf(float maxval) {
   
   int i = rand();
   float r = ((float) i)/((float) RAND_MAX);
   return r*maxval;
   
}

int randomi(int maxval) {
   
   int i = (int) randomf(maxval);
   return min(i,maxval-1);
   
}

int round(float f)
{
  int i= (int)(f+0.5);
  return i;
}
/* Monte Carlo (Metropolis) */
bool MC(float s1, float s2, float T ) {

  if(s1>s2) return true;

  if(s1<s2 && T<0 ) return false;
  
  float p=randomf(T);
  //cout<<T<<" randomf(T) "<<p<<" "<<(s2-s1)<<endl;
  if(p< exp(s1-s2)) return true;
  else return false;
}

int min_index(float* s, int size) 
{
  int   index=0;
  float sm=s[index];

  for(int i=1; i<size; i++) {
    if(sm>s[i]) {
      index=i;
      sm=s[index];
    }
  }
  return index;
}

float parabolic(float xa, float xb, float xc, 
			  float fa, float fb, float fc)
{
  float r,q,u;
  float TINY=1.0e-10;

  r=(xb-xa)*(fb-fc);
  q=(xb-xc)*(fb-fa);

  if(fabs(q-r) < TINY) {
    if((q-r) >0) {
      u=xb-((xb-xc)*q-(xb-xa)*r)/(2.0*TINY);
    }
    else {
      u=xb-((xb-xc)*q-(xb-xa)*r)/(-2.0*TINY);
    }
  }
  else {
    u=xb-((xb-xc)*q-(xb-xa)*r)/(2.0*(q-r));
  }

  return u;

}

string itos (int i) {
  ostringstream oss;
  oss<<i;
  return oss.str();
}
