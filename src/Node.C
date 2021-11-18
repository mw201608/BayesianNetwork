#include <math.h>
#include "Node.h"
#include "common.h"

string Node::getName() {
  return name;
}

void Node::setType(string nodeType) {
  if (nodeType == "discrete") {
    isDiscrete = true;
  }
  else if(nodeType == "continuous") {
    isDiscrete = false;
  }
  else {
    cout<<"unknown data type "<<nodeType <<" for node "<<name<<endl;
    exit(0);
  }
}

void Node::addData(float d) {
  data.push_back(d);
}

void Node::addNotion(string notion) {
  notions.push_back(notion);
}

void Node::calculateLikelihood() {
  likelihood=0;

  for(int i=0; i<nNotions(); i++) {
      likelihood+=count[i]*log(prior[i]);
  }
  
  dirtyF=false;
}

void Node::setTableCount(float** _table, int** _tablecount, int _tablesize) {
    //free old space first
    if(table != NULL) {
        FREEP(table, tablesize);
    }
    if(tablecount !=NULL) {
        FREEP(tablecount, tablesize);
    }
    table = _table;
    tablecount = _tablecount;
    tablesize = _tablesize;
}

float Node::getTable(int i, int j) {
    return table[i][j];
}

int Node::getTableSize() {
  return tablesize;
}

int Node::getTableCount(int i, int j) {
  return tablecount[i][j];
}

float Node::mutualInfo(Node* n) {
  if(!isDiscrete) {
    return fabs(correlation(n));
  }
  float p[nNotions()][n->nNotions()];
  //set pseudo counts
  for(int i=0; i<nNotions(); i++) {
    for(int j=0; j<n->nNotions(); j++) {
      p[i][j]=0.01;
    }
  }
  //count
  for(int d=0; d<data.size(); d++) {
    p[(int)data[d]][(int)(n->getData(d))]+=1.0;
  }
  //calculate probability
  float sum=0;
  for(int i=0; i<nNotions(); i++) {
    for(int j=0; j<n->nNotions(); j++) {
      sum+=p[i][j];
    }
  }  
  for(int i=0; i<nNotions(); i++) {
    for(int j=0; j<n->nNotions(); j++) {
      p[i][j]=p[i][j]/sum;
    }
  }    
  //calculate entropy
  float hxy=0;
  for(int i=0; i<nNotions(); i++) {
    for(int j=0; j<n->nNotions(); j++) {
      hxy-=p[i][j]*log(p[i][j]);
    }
  }   
  float hx=0;
  for(int i=0; i<nNotions(); i++) {
    hx-=prior[i]*log(prior[i]);
  }
  float hy=0;
  for(int i=0; i<n->nNotions(); i++) {
    hy-=n->getPrior(i)*log(n->getPrior(i));
  }  
  
  //using abs to take calculation error
  return fabsf(hx+hy-hxy);
}

/*
 * calculate conditional mutual information (cmi)
 * I(X;Y|Z) = sum(P(xi,yj,zk)log(p(xi,yj|zk)/(P(xi|zk)*p(yj|zk)))
 */
float Node::cmi(Node* n, Node* m) {
  float p[nNotions()][n->nNotions()][m->nNotions()];
  //set pseudo counts
  for(int k=0; k<nNotions(); k++) {
    for(int i=0; i<n->nNotions(); i++) {
      for(int j=0; j<m->nNotions(); j++) {
	p[k][i][j]=0.01;
      }
    }
  }

  //count
  for(int d=0; d<data.size(); d++) {
    p[(int)data[d]][(int)n->getData(d)][(int)m->getData(d)]+=1.0;
  }

  //calculate probability
  float sum=0;
  for(int k=0; k<nNotions(); k++) {
    for(int i=0; i<n->nNotions(); i++) {
      for(int j=0; j<m->nNotions(); j++) {
	sum+=p[k][i][j];
      }
    }
  }  

  //joint probability p(i,j,k)
  float J[nNotions()][n->nNotions()][m->nNotions()];
  for(int k=0; k<nNotions(); k++) {
    for(int i=0; i<n->nNotions(); i++) {
      for(int j=0; j<m->nNotions(); j++) {
	J[k][i][j]=p[k][i][j]/sum;
      }
    }    
  }

  //conditional joint probability p(i,j|k)
  float pk[nNotions()][n->nNotions()][m->nNotions()];
  for(int k=0; k<nNotions(); k++) {
    float s=0;
    for(int i=0; i<n->nNotions(); i++) {
      for(int j=0; j<m->nNotions(); j++) {
	s +=p[k][i][j];
      }
    } 
    //
    for(int i=0; i<n->nNotions(); i++) {
      for(int j=0; j<m->nNotions(); j++) {
	pk[k][i][j]=p[k][i][j]/s;
      }
    } 
  }

  //marginal conditional probability p(i|k) and p(j|k)
  float pi[nNotions()][n->nNotions()], pj[nNotions()][m->nNotions()];
  for(int k=0; k<nNotions(); k++) {
    sum = 0;
    for(int i=0; i<n->nNotions(); i++) {
      for(int j =0; j<m->nNotions(); j++) {
	sum += p[k][i][j];
      }
    }
    //pi[k][i]
    for (int i=0; i<n->nNotions(); i++) {
      float s = 0;
      for(int j=0; j<m->nNotions(); j++) {
	s +=p[k][i][j];
      }
      pi[k][i] = s/sum;
    }
    //pj[k][j]
    for (int j=0; j<m->nNotions(); j++) {
      float s = 0;
      for(int i=0; i<n->nNotions(); i++) {
	s +=p[k][i][j];
      }
      pj[k][j] = s/sum;
    }
  }

  sum=0;
  for(int k=0; k<nNotions(); k++) {
    for(int i=0; i<nNotions(); i++) {
      for(int j=0; j<nNotions(); j++) {
	sum +=J[k][i][j]*log(pk[k][i][j]/(pi[k][i]*pj[k][j]));
      }
    }
  }

  return sum;
}


float Node::getEntropy() {
  //calculate entropy
  float hx=0;
  for(int i=0; i<nNotions(); i++) {
    hx-=prior[i]*log(prior[i]);
  }  

  return hx;
}

void Node:: calculatePrior() {
    int numChoice = nNotions();
    float pr[numChoice];
    //set pseudo counts
    for(int i=0; i<numChoice; i++) {
      //pr[i]=0.01;
       pr[i]=0.00001;
    }
    //calculate the MLE of prior
    for(int i=0; i<data.size(); i++) { 
        pr[(int)data[i]]+=1.0;
    }
    float sum=0; 
    for(int i=0; i<numChoice; i++) {
        sum+=pr[i];
    }

    //save the result
    for(int i=0; i<numChoice; i++) {
        count.push_back((int)pr[i]);
        pr[i]=pr[i]/sum;
	prior.push_back(pr[i]);
    } 
}

void Node:: printData() {
    cout<<name<<"\t";
    for(int i=0; i<data.size(); i++) { 
      cout<<data[i]<<"\t";
    }
    cout<<"\n";
}

//Pearson correlation
float Node::correlation(Node* n) {
  float d =0;
  float dd;
  for(int i=0; i<data.size(); i++) {
      d+=data[i]*n->getData(i);
  }

  d=(d-sum()*n->sum()/data.size());

  dd=sqrt((sumSquare()-sum()*sum()/data.size())*
	  (n->sumSquare()-n->sum()*n->sum()/data.size()));
  if(dd==0) {
    d=0;
  }
  else {
    d=d/dd;
  }
  return d;
}

float Node::sum() {
  if(sumv!=0) return sumv;

  sumv=0;
  for(int i=0 ; i<data.size(); i++) {
    sumv+=data[i];
  }

  return sumv;
}

float Node::sumSquare() {
  if(sumv2!=0) return sumv2;

  sumv2=0;
  for(int i=0 ; i<data.size(); i++) {
    sumv2+=data[i]*data[i];
  }
  return sumv2;
}
