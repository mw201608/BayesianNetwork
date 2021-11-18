#ifndef NODE_H
#define NODE_H

#include <string>
//#include <vector.h>
#include <vector>
#include "eQTL.h"
#include "common.h"
class Node {

public: 
  Node(string name) {
      this->name = name;
      eqtl = eQTL();
      sumv=0;
      sumv2=0;
      dirtyF = true;
      table =  NULL;
      tablecount = NULL;
  }

  //Node's name
  void   setName(string name) {this->name = name; }
  string getName() ;

  //set node type discrete/continuous
  void setType(string nodetype);
  bool isDiscreteNode() {return isDiscrete;}
  //the node can have discrete values. Each value has a meaning
  void   addNotion(string notion);
  string getNotion(int i) {return notions[i];}
  // number of notions is also used for counting how many discrete values
  // the node can have
  int    nNotions() {return notions.size();}

  //eQTL
  eQTL*  getEQTL() {return &eqtl;}

  //data
  void   addData(float dvalue);
  void   printData();
  vector<float>* getData() {return &data;}
  float  getData(int n) {return data[n];}
  void initOrder();
  vector<float>* getOrder() {return &order;}
  float getOrder(int n) {return order[n];}
  float sum();
  float sumSquare();

  //statistics of data
  void   calculateLikelihood();
  void   setLikelihood(float liken) {likelihood=liken; dirtyF=false;}
  double getLikelihood() {return likelihood;}
  float  mutualInfo(Node* n);
  float  cmi(Node* n, Node* m);
  float  correlation(Node* n);
  float  spearCorrelation(Node* n);
  void   calculatePrior();
  float  getEntropy();
  float  getPrior(int n) { return prior[n];}
  void   setTableCount(float**, int**, int);
  int    getTableSize();
  float  getTable(int i, int j);
  int    getTableCount(int i, int j);

  //flag to check whether the likelihood is updated
  bool   isDirty() {return dirtyF; }
  void   setDirtyF() {dirtyF=true;}

private : 
  string  name;
  bool    isDiscrete; 
  vector  <float> prior;
  vector  <int> count;
  vector  <float> data;
  vector  <float> order;
  vector  <int> transitionPoints;
  vector  <string> notions; // the meaning of each value;
  eQTL    eqtl;
  float   sumv;
  float   sumv2;

  bool    dirtyF;
  double  likelihood;
  float** table; //conditional probability
  int**   tablecount; 
  int     tablesize;
};
#endif
