#ifndef GRAPH_H
#define GRAPH_H

class Graph {
public:
  Graph(){}
public:
  static void  topological_sort(bool** dag, int* order, int dim);
  static void  init_ancestor_matrix(bool** A, bool** dag, int dim);
  static void  update_node(bool** A, int j, bool** dag, int dim);
  static void  do_addition(bool** A, int p, int c, int dim);
  static void  do_removal(bool** A, int p, int c, bool** dag, int dim);
};

#endif
