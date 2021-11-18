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

#ifndef BN_H
#define BN_H

#include <math.h>
#include <vector>
#include <map>
#include "Node.h"
#include "Relationship.h"

/**
 *$Revision: 1.12 $
 *$Date: 2004/05/21 00:18:35 $
 *Author:  Jun Zhu
 *Creation Date: 3/23/01
 */

class BN {

public:
  BN();

public:
  int init(string biffile);
  void initSpace();
  void initPrior();
  void initPrior(float);
  void initPrior_qtlmatch();
  void initPrior_causal();
  void initPrior_complexity();
  void initPrior_peak();
  void updateNodesInfo();
  void setName(string name) {this->name =  name;}
  void setDataDim(int dim) {dataDim = dim; updateAlpha();}
  void setPseudoFlag(bool flag) {pseudoFlag = flag;}
  void setAlpha(float alpha) { this->alpha = alpha;}
  void setACutoff(float acutoff) {this->acutoff = acutoff;}
  void setQCutoff(float qcutoff) {this->qcutoff = qcutoff;}
  void updateAlpha() { alpha = alpha*log((double)dataDim);}
  void readdata();
  void readCorr(string corrfile);
  void permute();  //permute index

  //nodes
  void addNode(Node n) { nodes.push_back(n);}
  Node* getNodeByName(string name);
  int   getNodeIdByName(string name);
  Node* getNode(int i) { return &nodes[i];}
  void setSize(int size) {this->size = size;}
  int  getSize() {return size;}

  //trylist
  void setMaxTryList(int maxN) {maxTrylist = maxN; 
                               cout<<"maxTrylist "<<maxTrylist<<endl;}

  void loadTryList(string trylistFile);
  void updateTryList(int nodeId);
  void outputTryList(string trylistFile);

  //structure
  void  initStructure_reachable();
  void  learnStructure();
  void  enumerateStructures();
  void  enumerateStructures(string enumerateFile);
  double likelihood() ;
  double likelihood(int i);
  void inference_old(string inferencefile, int iDim);
  void inference(string inferencefile, int iDim, string outfile);
  void predictChange(vector<string> changeNode, vector<int> changeState);
  void predictOutcome(vector<string> changeNode, vector<int> changeState,
		      string, int);
  void predictOutcome(int );
  int pickState(int i);
  int simulateState(int i, bool* simulated, int* state) ;
  int inferenceState(int i, int j, int );
  float BIC();
  float BIC(int);
  float BIC(int, bool);
  void  printStructure();
  void  toDotInput(string outfile);
  void  toXMLBIF(string xmlfile);
  void  toBIF(string BIFile);
  void  toBNT(string BNTfile);
  void  refine(float t);
  int** blockpath(int**, int, int);
  int** shortestPath();
  int** readPath(string pathfile);
  void  printShortestPath(int**);
  Relationship weakestRelation(int, int, int**);
  Relationship weakestRelation(int);
  void  findMarkovBlanket(int );

  //relationship
  bool hasRelation(int pId, int cId);
  void addRelation(int pId, int cId);
  void addPotentialRelation(int pId, int cId);
  void addFixedRelation(int pId, int cId);
  void removeRelation(int pId, int cId);
  void reverseRelation(int pId, int cId);
  void updatePList(int cId, int pId, int operation);
  bool hasCycle(int n);
  int  getNumOfParent(int nodeId);

  //prior
  void  setPrior(Relationship* r);
  float getPrior(Node* p, Node* c);
  float getPrior(int pId, int cId);
  double getScalingPrior(int pId, int cId);
  double getPeakPrior(int pId, int cId);
  short getCC(int pId, int cId);
  void  addFromNode(int id);
  void  setReachable(Relationship* r);
  float adjustReachable();
  void  setCloseness(Relationship* r);

  //deal with node
  void updateNodeTable(int i, bool);
  void printNodeTable(int i);
  void calculateNodeLikelihood(int i);
  float matcheQTL(eQTL* p, eQTL* c);
  float matcheQTL(eQTL* p, eQTL* c, int k, float* sumw);
  float cosine(eQTL* p, eQTL* c);

private:
  bool              updateNodesInfoF;
  bool              pseudoFlag;
  int               size;  //number of nodes in the network
  int               dataDim; //number of data points for each node
  float**           data;  // data
  int               numParams;  // number of paramerters to describe the network
  float*            eQTLPeak;
  vector <Node>     nodes;
  map <string, int> nodeIds;
  map <string, float> trylist;  //try list, the result is bic score. The key is
                                //no parent: "n" n is node id
                                //with parents: "n:p1:p2:..:pk", n is node id
                                // and p1,p2, ... pk are its parent id.

  bool*             discreteF;  // whether a node is discrete.
  int*              nDiscrete;  //number of notions that each node has
  int*              numNodeParam; // number of parameters to describe a node
  int*              permIdx;     //permuted index
  bool**            relations;
  bool**            potentialRelations;
  bool**            fixedRelations;
  float**           closeness;
  float**           reachable;
  vector <int>      fromNodes; //starting points for testing reachable nodes
  int*              numChildren;
  short**           pList; //parent list
  short*            numParent;
  int*              numPotentialParent;
  short**           cc;
  float**           priors;
  float**           corr;    // correlation coefficients
  float             acutoff; //association cutoff
  float             qcutoff; //QTL score cutoff
  string            name;
  double            llike;  //log likelihood        
  double            alpha;  // use in BIC to penalize complexity 
  int               maxTrylist;
  const static int MAXNUMPARENTS = 3;
};
#endif
