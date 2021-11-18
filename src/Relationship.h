#ifndef RELATIONSHIP_H
#define RELATIONSHIP_H

#include "Node.h"
class Relationship {

public:
  Relationship() {}
  Relationship(int p, int c) {parent=p; child=c;}

  void  setChild(int c) { child=c;}
  int   getChild() ;
  void  setParent(int p) {parent=p;}
  int   getParent() ;
  void  setPrior(float p);
  float getPrior();

private:
  int parent;
  int child;
  float prior;
};
#endif
