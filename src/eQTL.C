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

#include "eQTL.h"
#include <math.h>
#include "common.h"


void eQTL::add(float loc, float lod, float r2) {
  locus.push_back(loc);
  lods.push_back(lod);
  R2.push_back(r2);
}
 
float eQTL::getLocus(int c) {
  return locus[c];
}

float eQTL::getLod(int c) {
  return lods[c];
}

float eQTL::getR2(int c) {
  return R2[c];
}
void eQTL::printMode(float cutoff) {
  if(complexity<0) getComplexity(cutoff);
 
  for(int i=1; i<modes.size(); i++) {
    if(modes[i]>cutoff) {
      cout<<i<<": "<<log(modes[i])<<" ";
    }
  }
  cout<<endl;
}

float eQTL::getMode(int chrom) {
  if(complexity<0) getComplexity(2);
  return modes[chrom];
}

float eQTL::getPeak() {
  if(peak<0) getComplexity(2);
  return peak;
}



int eQTL::getComplexity(float cutoff) {

  //check whether it has been initialized
  if(complexity >=0) return complexity;
  modes.push_back(0);
  for(int i=0; i<locus.size(); i++) {
    if(int(locus[i]/1000)>modes.size()-1) {
      modes.push_back(lods[i]);
    }
    else if(lods[i] > modes[int(locus[i]/1000)]) {
      modes[int(locus[i]/1000)] = lods[i];
    }
    if(lods[i]>peak) peak=lods [i];
  }

  complexity = 0;
  for(int i=1; i<modes.size(); i++) {
    if(modes[i]>cutoff) complexity++;
  }

  /* if there is no significant QTL (complexity=0), 
   * then set complexity to 12 so that it can almost
   * everyone's child.
   */
  if(complexity==0) complexity = 12;

  return complexity;
}

