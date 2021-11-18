#ifndef EQTL_H
#define EQTL_H
using namespace std;
#include <vector>
class eQTL {

public: 
  eQTL() {complexity = -1; peak=-1;}

  void add(float loc, float lod, float R2);
  float getLocus(int c);
  float getLod(int c);
  float getR2(int c);
  int   size() { return locus.size();}
  void  printMode(float cutoff);
  float getMode(int chrom);
  float getPeak();
  int   getComplexity(float cutoff);

private : 
  int  complexity ;
  float peak; 
  vector<float> locus;
  vector<float> lods;
  vector<float> R2;    //percentage of variance explained by the genetics
  vector<float> modes;  
};
#endif
