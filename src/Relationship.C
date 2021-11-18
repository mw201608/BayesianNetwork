#include "Relationship.h"

int Relationship::getChild() {
  return child;
}

int Relationship::getParent() {
  return parent;
}

void Relationship::setPrior(float p) {
  prior = p;
}

float Relationship::getPrior() {
  return prior;
}
