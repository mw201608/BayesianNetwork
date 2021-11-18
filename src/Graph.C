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

#include <vector>

#include "Graph.h"
#include "common.h"

/*
 * a collection of graph operation
 */

void  Graph::init_ancestor_matrix(bool** A, bool** dag, int n)
{
  int* order;
  NEW(order, n, int);
  topological_sort(dag, order, n);
  for (int j=0; j<n; j++) {
    update_node(A, order[j], dag, n);
  }
  free(order);
}

// Return the nodes in topological order (parents before children).
void Graph::topological_sort(bool** dag, int* order, int n)
{
  int indeg[n];

  //in-degree
  for (int i=0; i<n; i++) {
    indeg[i] = 0;
  }

  vector<int> zero_indeg;
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (dag[j][i]) {
	indeg[i]++;
      }
    }
    if (indeg[i] ==0) {
      zero_indeg.push_back(i);
    }
  }

  //initialize order
  for (int i=0; i<n; i++) {
    order[i] = 0;
  }

  int t=0;
  while(zero_indeg.size()>0) {
    int v = zero_indeg[0];
    zero_indeg.erase(zero_indeg.begin());
    order[t] = v;
    t = t + 1;
    for(int c=0; c<n; c++) {
      if (dag[v][c]) {
	indeg[c]--;
	if (indeg[c] == 0) {
	  zero_indeg.push_back(c);
	}
      }
    }
  } 
}

//update ancester matrix
void Graph::update_node(bool** A, int j, bool** dag, int n)
{
  for (int i=0; i<n; i++) {
    A[i][j] = false;
  }

  
  for (int i=0; i<n; i++) {
    if (dag[i][j]) {
      A[i][j] = true;
      //parent's parents
      for(int k=0; k<n; k++) {
	if(A[k][i]) { 
	  A[k][j] = true;
	}
      }
    }
  }
}


//update ancester matrix after adding a link
void Graph::do_addition(bool** A, int p, int c, int n)
{
  A[p][c] = 1; 
  for(int i=0; i<n; i++) {
    if(A[i][p]) {
      A[i][c] = 1;
    }
  }

  for (int j=0; j<n; j++) {
    if(A[c][j]) {
      for (int i=0; i<n; i++) {
	if(A[i][c]) {
	  A[i][j] = true;
	}
      }
    }
  }
}

//update ancestor matrix after a removal
void Graph::do_removal(bool** A, int p, int c, bool** dag, int n)
{
  int* order;
  NEW (order, n, int);
  topological_sort(dag, order, n);

  bool desc[n];
  //find all the descendants of c, and put them in topological order
  for (int i=0; i<n; i++) {
    if (A[c][i]) {
      desc[i]  = true;
    }
    else {
      desc[i] = false;
    }
  }

  update_node(A, c, dag, n);
  for (int i=0; i<n; i++) {
    if(desc[order[i]]) {
      update_node(A, order[i], dag, n);
    }
  }
  
  free(order);
}



