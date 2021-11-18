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

//general include information
#include <map>
#include <fstream>
#include <iostream>
#include <strstream>
//Include files for xml parser
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
//local include files
#include "BNHandler.h"
#include "BN.h"
#include "Gaussian.h"
#include "Graph.h"
#include "common.h"

BN::BN() {
  relations = NULL;
  potentialRelations = NULL;
  fixedRelations = NULL;
  closeness = NULL;
  reachable = NULL;
  pList = NULL;
  numParent = NULL;
  numChildren = NULL;
  numPotentialParent = NULL;
  priors = NULL;
  corr = NULL;
  updateNodesInfoF = true; //need update nodes info after scanning all nodes
  maxTrylist = 150000000;
}

/*
 * read in a file in BIF format
 * and use it to initialize the Bayesian network.
 */
int BN::init(string biffile) {
    // Initialize the XML utility
    try
    {
         XMLPlatformUtils::Initialize();
    }

    catch (const XMLException& toCatch)
    {
         cerr << "Error during initialization! :\n"
              << toCatch.getMessage() << endl;
         return 1;
    }

    //  Create a SAX parser object. 
    SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();

    #ifdef VALIDATION
    // always report validation errors
    parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
    parser->setFeature(XMLUni::fgXercesDynamic, false);
    #else
    //  Do not report validation errors
    parser->setFeature(XMLUni::fgSAX2CoreValidation, false);
    #endif
    

    bool doNamespaces    = true;
    bool doSchema        = true;
    bool schemaFullChecking = false;
    bool namespacePrefixes = false;
    
    parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, doNamespaces);
    parser->setFeature(XMLUni::fgXercesSchema, doSchema);
    parser->setFeature(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking);
    parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes, namespacePrefixes);

    //
    //  Create the handler object and install it as the document and error
    //  handler for the parser. Then parse the file and catch any exceptions
    //  that propogate out
    //

    int errorCount = 0;
    try
    {
        BNHandler handler(this);
	parser->setContentHandler(&handler);
        parser->setErrorHandler(&handler);
        parser->parse(biffile.c_str());
        errorCount = parser->getErrorCount();
    }

    catch (const XMLException& toCatch)
    {
        cerr << "\nAn error occurred\n  Error: "
             << toCatch.getMessage()
             << "\n" << endl;
        XMLPlatformUtils::Terminate();
	return 1;
    }

    //  Delete the parser itself.  Must be done prior to calling Terminate.
    delete parser;

    // And call the termination method
    XMLPlatformUtils::Terminate();

    return 0;
}

/*
 * after reading in all nodes information, we need to update nodes' info.  
 * Also register all nodes' ID so that we can get node's ID from name
 */
void BN::updateNodesInfo() {
    if(!updateNodesInfoF) return;
    updateNodesInfoF = false;

    //after scanning nodes, do some sanitary check
    if(size != nodes.size()) {
      cout<<"Error: network size is "<<size
	  <<", but nodes size is "<<nodes.size()<<endl;
      exit(0);
    }
    //after scanning nodes, we update how many discrete values 
    //each node has
    numParams=0;
    for(int i=0; i<size; i++) {
      nDiscrete[i] = nodes[i].nNotions();
      numNodeParam[i] = nDiscrete[i]-1;
      numParams += numNodeParam[i];
    }

    //register node's id
    for(int i=0; i<size; i++) {
      nodeIds[nodes[i].getName()] = i;
    }
}

/* locate space needed */
void BN::initSpace() {
    if(priors == NULL || relations == NULL || 
       potentialRelations == NULL || 
       fixedRelations == NULL || closeness==NULL || reachable == NULL ||
       numParent == NULL || numChildren == NULL ||
       numPotentialParent == NULL ||
       nDiscrete == NULL || numNodeParam == NULL) {

        NEWP(data, size, float); 
	NEWP(priors, size, float);
	NEWP(corr, size, float);
	NEWP(cc, size, short);
	NEWP(relations, size, bool);
	NEWP(potentialRelations, size, bool);
	NEWP(fixedRelations, size, bool);
	//NEWP(closeness, size, float);
	NEWP(reachable, size, float);
	NEWP(pList, size, short);
	NEW(numParent, size, short);
	NEW(numChildren, size, int);
	NEW(numPotentialParent, size, int);
	NEW(discreteF, size, bool);
	NEW(nDiscrete, size, int);
	NEW(numNodeParam, size, int);
	NEW(eQTLPeak, size, float);	
	NEW(permIdx, size, int);
	
	//initialize permutation index
	for(int i=0; i<size; i++) {
	  permIdx[i] = i;
	}

	for(int i=0; i<size; i++) {
	  numParent[i] = 0;
	  numChildren[i] = 0;
	  numPotentialParent[i] = 0;
	}

	//initialize relations
	for(int i=0; i<size; i++) {
	    NEW(data[i], dataDim, float);
	    NEW(priors[i], size, float);
	    NEW(corr[i], size, float);
	    NEW(cc[i], size, short);
	    NEW(relations[i], size, bool);
	    NEW(potentialRelations[i], size, bool);
	    NEW(fixedRelations[i], size, bool);
	    //NEW(closeness[i], size, float);
	    NEW(reachable[i], size, float);
	    NEW(pList[i], MAXNUMPARENTS, short);
	    for(int j=0; j<size; j++) {
	        //-1 means it has not been set
		priors[i][j] = -1;
		//0 mean it has not been set
		cc[i][j] = 0;
		relations[i][j] = false;
		potentialRelations[i][j] = false;
		fixedRelations[i][j] = false;
		//closeness[i][j] = 0;
		reachable[i][j] = 1;
	    }
	}
    }
    cout<<"end of initSpace"<<endl;
}
/*initialize priors to non-informative
 */
void BN::initPrior() {
  float mi, me;
  for(int i=0; i<size; i++) {
      for(int j=i+1; j<size; j++) {
	priors[i][j] = log(1.0/size);
	priors[j][i] = log(1.0/size);
      }
  }
}

/*initialize prior according to
 *eQTL's overlapping result.
 */
void BN::initPrior(float te) {
  float me;
  //float scale = 1.0/(float(size)*te);
  float scale = 1.0;
  for(int i=0; i<size; i++) {
      cerr<<i<<" "<<te<<endl;
      for(int j=i+1; j<size; j++) {
	//cout<<i<<" "<<j<<endl;
	//only the priors that have not been set by init structure will
	//be set here.
	if(priors[i][j] <0 || priors[j][i] <0) {
	  //me =  matcheQTL(nodes[i].getEQTL(), nodes[j].getEQTL());
	  me =1;
	  //cout<<"me "<<me<<endl;
	  //make the correlation non-symmetric 
	  double rij  = getScalingPrior(i,j);
	  double rji  = getScalingPrior(j,i);
	  double s = rij+rji;
	  rij = rij/s;
	  rji = rji/s;

	  if(priors[i][j] <0 && rij*me*fabsf(reachable[i][j]) >=te)  {
	    if(reachable[j][i]>1) {
	      priors[i][j] = -100000;
	    }
	    else {
	      if(rij*me*scale*fabsf(reachable[i][j])==0) {
		priors[i][j] = -10000;
	      }
	      else {
		priors[i][j] = log(rij*me*scale*fabsf(reachable[i][j]));
	      }
	    }
	  }
	  else if(priors[i][j] <0 ||priors[i][j] ==0) {
	    priors[i][j] = -100000;
	  }
	  else {
	    priors[i][j] = log(priors[i][j]);
	  }
	  if(priors[j][i] <0 && rji*me*fabsf(reachable[j][i]) >=te)  {
	    if(reachable[i][j] >1) {
	      priors[j][i] = -100000;
	    }
	    else {
	      if(rji*me*scale*fabsf(reachable[j][i])==0) {
		priors[j][i] = -10000;
	      }
	      else {
		priors[j][i] = log(rji*me*scale*fabsf(reachable[j][i]));
	      }
	    }
	  }
	  else if(priors[j][i] <0 || priors[j][i] ==0) {
	    priors[j][i] = -100000;
	  }
	  else {
	    priors[j][i] = log(priors[j][i]);
	  }
	}
	else {
	  //re-scale to log scale
	  if(priors[i][j]==0) {
	    priors[i][j] = -100000;
	  }
	  else {
	    priors[i][j] = log(priors[i][j]);
	  }
	  if(priors[j][i]==0) {
	    priors[j][i] = -100000;
	  }
	  else {
	    priors[j][i] = log(priors[j][i]);
	  }
	}
      }
  }
}


/*initialize prior according to
 *eQTL's overlap
 */
void BN::initPrior_qtlmatch() {
  float me;
  for(int i=0; i<size; i++) {
      for(int j=i+1; j<size; j++) {
	  me =  matcheQTL(nodes[i].getEQTL(), nodes[j].getEQTL());
	  //in case there is no significant QTL in one of them
	  if(me==0 & (nodes[i].getEQTL()->getComplexity(qcutoff)>=12 ||
	              nodes[j].getEQTL()->getComplexity(qcutoff)>=12)) {
	    me =0.135;
	  }
	  else if (me<0.01) {
	    me = 0.01;
	  }
	  priors[i][j] = log(me);
	  priors[j][i] = log(me);
      }
  }
}

/*initialize prior according to
 *eQTL LOD score.
 */
void BN::initPrior_causal() {
  for(int i=0; i<size; i++) {
      for(int j=i+1; j<size; j++) {
	  double rij  = getScalingPrior(i,j);
	  double rji  = getScalingPrior(j,i);
	  double s = 0.5*(rij+rji);
	  rij = rij/s;
	  rji = rji/s;
	  priors[i][j] = log(rij);
	  priors[j][i] = log(rji);
      }
  }
}

/*initialize prior according to peak
 *eQTL LOD score.
 */
void BN::initPrior_peak() {

  //initialization
  for(int i=0; i<size; i++) {
      for(int j=0; j<size; j++) {
	   priors[i][j] = 0;
      }
  }

  //get Peak postion
  for(int i=0; i<size; i++) {
    eQTLPeak[i] = nodes[i].getEQTL()->getPeak();
  }

  for(int i=0; i<size; i++) {
    if(eQTLPeak[i] >qcutoff) { 
      for(int j=i+1; j<size; j++) {
	if(eQTLPeak[j]>qcutoff) {
	  double rij  = getPeakPrior(i,j);
	  double rji ;
	  if (rij !=0) {
	    rji = 1/rij;
	    priors[i][j] = log(rij);
	    priors[j][i] = log(rji);
	  }
	}
      }
    }
  }
}

/*initialize prior according to
 *eQTL complexity
 */
void BN::initPrior_complexity() {
  for(int i=0; i<size; i++) {
      cerr<<i<<endl;
      for(int j=i+1; j<size; j++) {
	//make the correlation non-symmetric 
	  double r ;
	  if(nodes[i].getEQTL()->getComplexity(qcutoff)+
	     nodes[j].getEQTL()->getComplexity(qcutoff)==0) {
	    r = 0.5;
	  }
	  else {
	    if(nodes[i].getEQTL()->getComplexity(qcutoff)==0) {
	      r = double(nodes[j].getEQTL()->getComplexity(qcutoff))/
		(2*nodes[j].getEQTL()->getComplexity(qcutoff)+1);
	    }
	    else if(nodes[j].getEQTL()->getComplexity(qcutoff)==0){
	      r = double(nodes[i].getEQTL()->getComplexity(qcutoff)+1)/
		(2*nodes[i].getEQTL()->getComplexity(qcutoff)+1);
	    }
	    else {
	      r = double(nodes[j].getEQTL()->getComplexity(qcutoff))/
		  (nodes[i].getEQTL()->getComplexity(qcutoff)+
		   nodes[j].getEQTL()->getComplexity(qcutoff));
	    }
	  }
	  double rij = r;
	  double rji = 1-r;
	  priors[i][j] = log(rij);
	  priors[j][i] = log(rji);
      }
  }
}

//pull data from individual node into a central matrix for fast access
void BN::readdata() {
  for (int i=0; i<size; i++) {
    for(int j=0; j<dataDim; j++) {
      data[i][j] = nodes[i].getData(j);
    }
    discreteF[i] = nodes[i].isDiscreteNode();
  }
}

//read in correlation matrix
void BN::readCorr(string corrfile) {
    //read in data
    ifstream cfile(corrfile.c_str()); 
    if(!cfile) {
        cout <<"Can not open correlation file "<<cfile<<"\n";
	exit(-1);
    }
    else {
	//read in data
        int maxlen = 81920;
	char line[maxlen];
	int i=0;
	while(cfile.getline(line,maxlen,'\n') ) {
	    istrstream istr(line);
	    for(int j=0; j<size; j++) {
	        istr>>corr[i][j];
	    }
	    i++;
	}
	cfile.close(); 
    }
}


/* permute index */
void BN::permute() {
  int p, n, k;
  int m = size*5;
  for (int i=0; i<size/2; i++) {
    p=randomi(m)%size;
    n= randomi(m)%size;
    //switch 
    k=permIdx[p];
    permIdx[p]=permIdx[n];
    permIdx[n]=k;
  }

}

/*
 * load trylist
 */
void BN::loadTryList(string trylistFile) {
  ifstream tfile(trylistFile.c_str());
  if(!tfile) {
    cout <<"Can not open trylist file "<<trylistFile<<"\n";
    exit(-1);
  }
  else {
    int maxlen=1024;
    char line[maxlen];
    //read in trylist
    while(tfile.getline(line,maxlen,'\n') ) {
      istrstream istr(line);
      string key;
      istr>>key;
      float pp;
      istr>>pp;
      if(trylist.count(key)==0) {
	trylist[key] = pp;
      }
    }
  }
  tfile.close();
}

float BN::BIC() {
  float b=0;
  for (int i=0; i<size; i++) {
    b += BIC(i, true);
  }
  return b;
}

float BN::BIC(int nodeId) {
  return BIC(nodeId, false);
}

float BN::BIC(int nodeId, bool saveF) {
  float b;
  string key= itos(nodeId);
  
  int dx= numParent[nodeId];
  if (dx>0) {
    for (int i=0; i<numParent[nodeId]; i++) {
	key +=":"+itos(pList[nodeId][i]);
    }
  }
  if(trylist.count(key)==0) {
    if (discreteF[nodeId]) { //discrete
        b=-alpha*numNodeParam[nodeId]+likelihood(nodeId);
    }
    else { //continuous
      if(numParent[nodeId]==0) { //no parent 
	b = Gaussian::bic(dataDim, data[nodeId], alpha);
      }
      else { //with parents
	float** datax = matrix(0, dx-1, 0, dataDim-1);
	int np=0;
	for (int i=0; i<size; i++) {
	  if(relations[i][nodeId]) {
	    for (int k=0; k<dataDim; k++) {
	      datax[np][k] = data[i][k];
	    }
	    np++;
	  }
	}
        b = Gaussian::bic(dataDim, data[nodeId], dx, datax, alpha); 
	free_matrix(datax, 0,dx-1, 0, dataDim-1);
      }
    }
    //if(saveF || trylist.size()<150000000) {
    if(saveF || trylist.size()<maxTrylist) {
      trylist[key] = b;
    }
  }
  else {
    b = trylist[key];
  }
  return b;
}
/*
 * update trylist
 */
void BN::updateTryList(int nodeId) {
  float b;  //bic score
  //no children
  BIC(nodeId);

  //with parents
  for (int i=0; i<size; i++) {
    if(potentialRelations[i][nodeId]) {
      addRelation(i, nodeId);
      BIC(nodeId);
      /*
      //multiple parents
      for (int j=i+1; j<size; j++) {
	if(potentialRelations[j][nodeId]) {
	  addRelation(j, nodeId);
	  BIC(nodeId);
	  //recover
	  removeRelation(j, nodeId);
	}
      }
      */
      //recover
      removeRelation(i, nodeId);
    }
  }
}

void BN::outputTryList(string updateTrylistFile) {
  //output trylist
  ofstream out(updateTrylistFile.c_str(),ios::out);
  for (map<string,float>::iterator it = trylist.begin(); 
       it != trylist.end(); it++) {
      out << it->first << "\t" << it->second <<endl;
  }
  out.close();
}


/*
 * weighted average correlated across all chromosome
 */
float BN::matcheQTL(eQTL* ep, eQTL* ec) {

  int nChrom = int(ep->getLocus(ep->size()-1)/1000);
  
  float sum=0;
  float sumw=0;
  for(int k=1; k<=nChrom; k++) {
    sum += matcheQTL(ep, ec, k, &sumw);
  }
  return sum/(sumw);
}


float BN::matcheQTL(eQTL* ep, eQTL* ec, int K, float* sumw) {
  //eQTL could be in different length, make sure
  //match comparing the same thing
  int i = 0;
  double tmp;
  double weight =0;
  double m=0;
  double sump=0;
  double sumc=0;
  float cutoff = qcutoff;
  while(i<ep->size()) {
    int k = int(ep->getLocus(i)/1000);
    if(k==K ) {
	tmp = ep->getLod(i)*ec->getLod(i);
	if(weight < tmp) {
	  weight = tmp;
	}
	//if(ep->getLod(i)>qcutoff && ec->getLod(j)>qcutoff) {
	if(ep->getLod(i)>cutoff && ec->getLod(i)>cutoff) {
	  m += tmp;
	}
	sump  += ep->getLod(i)*ep->getLod(i);
	sumc += ec->getLod(i)*ec->getLod(i);
    }
    i++;
  }
   
  double sum =0;
  if(ep->getMode(K)>10) {
      weight = weight*10.0/ep->getMode(K);
  }
  if(ec->getMode(K)>10) {
    weight = weight*10.0/ec->getMode(K);
  }
	      
  //cout<<k<<" "<<weight[k]<<" "<<m[k]<<" "<<sqrt(sump[k]*sumc[k])
  //	<<" "<<m[k]/sqrt(sump[k]*sumc[k])<<endl;
  sum = weight*m/sqrt(sump*sumc);
  *sumw += weight;

  return sum;
}


double BN::getScalingPrior(int pnode, int cnode) {
  eQTL* ep = nodes[pnode].getEQTL();
  eQTL* ec = nodes[cnode].getEQTL();
  int nchrom = int(ep->getLocus(ep->size()-1)/1000);
  double ne[nchrom], nd[nchrom];

  for (int k=0; k<=nchrom; k++) {
      ne[k]=0;
      nd[k]=0;
  }

  //eQTL could be in different length, make sure
  //match comparing the same thing
  int i, j, n;
  i = 0;
  j = 0;
  while(i<ep->size() && j<ec->size()) {
    n = int(ep->getLocus(i)/1000);
    if(ep->getLocus(i) == ec->getLocus(j)) {
      if(ep->getLod(i) >qcutoff && ne[n]<ep->getLod(i)) {
	ne[n] = ep->getLod(i);
	/*every thing that below 2 is not reliable */
	if(ec->getLod(j)<2) {
	    nd[n] = 2;
	}
	else {
	    nd[n] = ec->getLod(j);
	}
      }
      i++; j++;
    }
    else if(ep->getLocus(i) > ec->getLocus(j)) {
      j++;
    }
    else {
      i++;
    }
  }
  double sum=1;

  for (int k=1; k<=nchrom; k++) {
    if (nd[k]>0) {
      sum = sum*ne[k]/nd[k];
    }
  }
  
  //cout<<pnode<<" "<<cnode<<" "<<sume<<" "<<sumd<<endl;
  return sum;
}

/* get Peak prior */
double BN::getPeakPrior(int pnode, int cnode) {
  eQTL* ep = nodes[pnode].getEQTL();
  eQTL* ec = nodes[cnode].getEQTL();
  int nchrom = int(ep->getLocus(ep->size()-1)/1000);

  //peak position on each chromosome
  double p[nchrom], c[nchrom];
  for (int k=0; k<=nchrom; k++) {
    p[k]=ep->getMode(k);
    c[k]=ec->getMode(k);
  }

  double sump=0;
  float sumw = 0;
  for (int k=1; k<=nchrom; k++) {
    if (p[k]>4.3 & c[k]>qcutoff) {
      sump +=matcheQTL(ep,ec,k, &sumw)*p[k]/c[k];
    }
  }
  double sumc=0;
  for (int k=1; k<=nchrom; k++) {
    if (c[k]>4.3 & p[k]>qcutoff) {
      sumc +=matcheQTL(ec,ep,k, &sumw)*c[k]/p[k];
    }
  }
    
  double sum= 1;
  if(sump+sumc>0) {
    sum = 2*sump/(sump+sumc);
  }

  return sum;
}


float BN::cosine(eQTL* ep, eQTL* ec) {
  long double sump=0;
  long double sumc=0;
  long double sum=0;
  //eQTL could be in different length, make sure
  //match comparing the same thing
  int i, j;
  i = 0;
  j = 0;
  while(i<ep->size() && j<ec->size()) {
    if(ep->getLocus(i) == ec->getLocus(j)) {
        sum += ep->getLod(i)*ec->getLod(j);
	sump  += ep->getLod(i)*ep->getLod(i);
	sumc += ec->getLod(j)*ec->getLod(j);
	i++; j++;
    }
    else if(ep->getLocus(i) > ec->getLocus(j)) {
      sum += ec->getLod(j);
      sump +=1;
      sumc += ec->getLod(j)*ec->getLod(j);
      j++;
    }
    else {
      sum += ep->getLod(i);
      sump +=ep->getLod(i)*ep->getLod(i);
      sumc +=1;
      i++;
    }
  }

  return sum/sqrt(sump*sumc);

}

/**********************************************************
 *
 *    priors
 *
 *********************************************************/
void BN::setPrior(Relationship* r) {
  #ifdef DEBUG
    cout<< "pId="<<r->getParent()<<endl;
    cout<< "cId="<<r->getChild()<<endl;
  #endif
  priors[r->getParent()][r->getChild()] = r->getPrior();
}

/* if a negative number is returned, then there is no prior information
 * is set for the relationship
 */
float BN::getPrior(Node* p, Node* c) {
  return priors[nodeIds[p->getName()]][nodeIds[c->getName()]];
}

float BN::getPrior(int pId, int cId) {
  return priors[pId][cId];
}

short BN::getCC(int pId, int cId) {
  //calculate CC first
  if(cc[pId][cId] ==0) {
    if(nodes[pId].correlation(&nodes[cId]) >0) {
      cc[pId][cId] = 1;
      cc[cId][pId] =1;
    }
    else {
      cc[pId][cId] = -1;
      cc[cId][pId] = -1;
    }
  }
  return cc[pId][cId];
}

/**********************************************************
 *
 *    reachable information
 * positive value: causal relationship
 * negative value: close relationship
 *
 *********************************************************/
void BN::setReachable(Relationship* r) {
    #ifdef DEBUG
    cout<< "pId="<<r->getParent()<<endl;
    cout<< "cId="<<r->getChild()<<endl;
    #endif
    //reachable information can affect prior information
    //if(priors[r->getParent()][r->getChild()]==0) {
    //  priors[r->getParent()][r->getChild()]=-1; 
    //}

    //always keep the highest number
    if(fabsf(reachable[r->getParent()][r->getChild()]) < fabsf(r->getPrior())) {
      reachable[r->getParent()][r->getChild()] = r->getPrior();
    }
}

/**********************************************************
 *
 *   closeness information
 *
 *********************************************************/
void BN::setCloseness(Relationship* r) {
    #ifdef DEBUG
    cout<< "pId="<<r->getParent()<<endl;
    cout<< "cId="<<r->getChild()<<endl;
    #endif
    //always keep the highest number
    if(closeness[r->getParent()][r->getChild()] < r->getPrior()) {
      closeness[r->getParent()][r->getChild()] = r->getPrior();
    }
}

void BN::addFromNode(int id) {
  bool containsF = false;
  for(int i=0; i<fromNodes.size(); i++) {
    if(fromNodes[i]==id) {
      containsF = true;
    }
  }
  if(!containsF) {
    fromNodes.push_back(id);
  }
}

// given a network and a set of from-nodes, penalize or reward reachabilities
float BN::adjustReachable() {
  float adjustment=0;
  map<int, int> waitList;
  map<int, int> checkedList;
  
  for(int i=0; i<fromNodes.size(); i++) {
    waitList[fromNodes[i]]=1;
  }
  
  while(waitList.size()>0) {
    int nodetmps = waitList.begin()->first;
    checkedList[nodetmps] =1;
    for(int i=0; i<size; i++) {
      if(relations[nodetmps][i]) {
	  if(checkedList.find(i)==checkedList.end() ) {
	    //add to waiting list
	    waitList[i]=1;
	  } 
	  adjustment+=reachable[nodetmps][i];
      }
    }
  }
  return adjustment;
}

/************************************************************
 *
 *        relations
 *
 ************************************************************/
void BN::addRelation(int pId, int cId) {
  if(!relations[pId][cId]) {
    relations[pId][cId] = true;
    numParent[cId]++;
    numChildren[pId]++;
    updatePList(cId, pId, 1);

    //update number of parameters needed
    numParams -= numNodeParam[cId];
    numNodeParam[cId] *= nDiscrete[pId];
    numParams += numNodeParam[cId];

    // set dirty flag
    nodes[cId].setDirtyF();
  }
}

void BN::addFixedRelation(int pId, int cId) {
   fixedRelations[pId][cId] = true;
   addRelation(pId, cId);
}

void BN::addPotentialRelation(int pId, int cId) {
    potentialRelations[pId][cId] = true;
    numPotentialParent[cId]++;
}

void BN::removeRelation(int pId, int cId) {
  relations[pId][cId] = false;
  numParent[cId]--;
  numChildren[pId]--;
  updatePList(cId, pId, -1);

  //update number of parameters needed
  if(discreteF[cId]) {
    numParams -= numNodeParam[cId];
    numNodeParam[cId] = numNodeParam[cId]/nDiscrete[pId];
    numParams += numNodeParam[cId];

    //set dirty flag
    nodes[cId].setDirtyF();
  }
}

bool BN::hasRelation(int pId, int cId) {
  return relations[pId][cId];
}


void BN::reverseRelation(int pId, int cId) {
  removeRelation(pId, cId);
  addRelation(cId, pId);
}

/********************************************
 * update parent list for a node
 ********************************************/
void BN::updatePList(int cId, int pId, int op) {
  int i,j;
  if(op<0) {//remove a parant
    for (i=0; i<=numParent[cId]; i++) {
      if(pList[cId][i] == pId) {
	//shift down
	for (j=i+1; j<numParent[cId]+1; j++) {
	  pList[cId][j-1] = pList[cId][j];
	}
      }
    }
  }
  else { //add a parant
    if(pId<pList[cId][0]) {
      j=0;
    }
    else if(pId>pList[cId][numParent[cId]-2]) {
      j= numParent[cId] -1;
    }
    else {
      i=0;
      while (pId<pList[cId][i]) {
	i++;
      }
      j=i;
    }

    //shift up
    for (i=numParent[cId]-1; i>j; i--) {
      pList[cId][i]=pList[cId][i-1];
    }
    pList[cId][j] = pId;
  }
}

/************************************************************
 *
 *   Node
 *
 ************************************************************/
int BN::getNodeIdByName(string name) {
    return nodeIds[name];
}

Node* BN::getNodeByName(string name) {
  Node* node=NULL;
  bool matched=false;
  for(vector<Node>::iterator i=nodes.begin(); 
      !matched && i!=nodes.end(); i++) {
    if(i->getName() == name ) {
      matched=true;
      node= &*i;
    }
  }
  return node;
}

void BN::printNodeTable(int i) {
    //find out how much space is needed
    int np3=1;
    for(int j=0; j<size; j++) {
      if(relations[j][i]) {
	cout<<nodes[j].getName()<<" ";
        np3=np3*nDiscrete[j];
      }
    }
    cout<<endl;

    for(int j=0; j<np3; j++) {
        cout<<j<<" ";
    	for(int k=0; k<nDiscrete[i]; k++) {
	  cout<< nodes[i].getTable(j, k)<<" ";
	}
	cout<<endl;
    }
}


void BN::updateNodeTable(int i, bool debug) {
    
    //find out how much space is needed and where is the mid point
    int np3=1;
    int constj=0; 
    for(int j=0; j<numParent[i]; j++) {
        np3=np3*nDiscrete[pList[i][j]];
	constj=constj*nDiscrete[pList[i][j]]+1;
    }

    //initialize table
    float** table;
    int** count;
    NEWP(table, np3, float);
    NEWP(count, np3, int);
    float pseudocount = dataDim/float(np3*nDiscrete[i]);
    if(pseudocount>sqrt((float)dataDim)) {
      pseudocount=sqrt((float)dataDim)*0.5;
    }
    if(!pseudoFlag) {
      pseudocount = 0;
    }
    
    for(int j=0; j<np3; j++) {
        NEW(table[j], nDiscrete[i], float);
	NEW(count[j], nDiscrete[i], int);
	for(int k=0; k<nDiscrete[i]; k++) {
	    //set pseudo counts
	    table[j][k]=pseudocount+0.01;
	    //init count
	    count[j][k]=0;
	}
    }
	
    //enforce result to be 1 when all input are ones
    table[constj][1] += pseudocount;

    /* if the relationship is one-to-one, enforce linear function */
    /*
    if(numParent[i] ==1) {
      for(int n=0; n<size; n++) {
	if(relations[n][i]) {
	  if(getCC(n, i) >0) {
	    table[0][0] +=pseudocount*0.5;
	    table[nDiscrete[n]-1][nDiscrete[i]-1] +=pseudocount*0.5;
	  }
	  else {
	    table[0][nDiscrete[i]-1] +=pseudocount*0.5;
	    table[nDiscrete[n]-1][0] +=pseudocount*0.5;
	  }
	}
      }
    }
    */

    //calculate the MLE of conditional probability
    for(int d=0; d<dataDim; d++) {
        int j=0;
	for(int n=0; n<numParent[i]; n++) {
	    j=j*nDiscrete[pList[i][n]]
	      +(int)data[pList[i][n]][d];
	}
	table[j][(int)data[i][d]]+=1.0;
	count[j][(int)data[i][d]]++;
    }

    //normalize
    for(int j=0; j<np3; j++) {
        float sum=0;
	for(int k=0; k<nDiscrete[i]; k++) {
	    sum+=table[j][k];
	}
	if(debug) {
	  cout<<j<<" ";
	  for(int k=0; k<nDiscrete[i]; k++) {
	    cout<<table[j][k]<<" ";
	  }
	  cout<<endl;
	}
	for(int k=0; k<nDiscrete[i]; k++) {
	    table[j][k]=table[j][k]/sum;
	    //cout<<table[j][k]<<" ";
	}
    }

    /*
    for(int j=0; j<np3; j++) {
	for(int k=0; k<nDiscrete[i]; k++) {
	    cout<<table[j][k]<<" ";
	}
	cout<<endl;
    }
    */

    nodes[i].setTableCount(table, count, np3);
}

void BN::calculateNodeLikelihood(int i) {
    float liketmp=0;
    //use conditional probability
    for(int j=0; j<nodes[i].getTableSize(); j++) {
	for(int n=0; n<nDiscrete[i]; n++) {
	  liketmp+=nodes[i].getTableCount(j,n)*log(nodes[i].getTable(j, n));
	}
    }

    //add prior info
    if(numParent[i]>0) {
      for(int j=0; j<numParent[i]; j++) {
	liketmp+=getPrior(pList[i][j], i);
      }
    }
    nodes[i].setLikelihood(liketmp);
}

/* 
 * calculate likelihood of the network
 *
 */
double BN::likelihood() {
  llike=0; 
  for(int i=0; i<size; i++) {
    llike+= likelihood(i) ;
  }
  /*
  if(fromNodes.size()>0) {
    float adjustment = adjustReachable();
    #ifdef DEBUG
    cout<<"adjustment: "<<adjustment<<endl;
    #endif
    llike+=adjustment;
  }
  */
  #ifdef DEBUG
  cout<<"Likelihood: "<<llike<<endl;
  #endif
  return llike;
}

double BN::likelihood(int i) {
  if(nodes[i].isDirty()) {
      if(numParent[i]>0) {
	  updateNodeTable(i, false);
	  calculateNodeLikelihood(i);
      }
      else {
	nodes[i].calculateLikelihood();
      }      
  }

  return nodes[i].getLikelihood();
}

/* given a state, calculate the probability of the state */
void BN::inference_old(string file, int iDim) {
    int maxlen = 81920;
    char line[maxlen];
    int  state[iDim][size];

    //read in data
    ifstream dfile(file.c_str()); 
    if(!dfile) {
        cout <<"Can not open data file "<<file<<"\n";
	exit(-1);
    }
    else {
	//read in data
	while(dfile.getline(line,maxlen,'\n') ) {
	    // comment line starting with #, so skip it
	    if(line[0] != '#') {
	        istrstream istr(line);
		string sName;
		int sValue;
		istr>>sName;
		for(int i=0; i<iDim; i++) {
		  istr>>sValue;
		  state[i][getNodeIdByName(sName)] = sValue;
		}
	    }
	}
	dfile.close();   
    }

    for(int n=0; n<iDim; n++) {
      float ilike = 0;
      for(int i=0; i<size; i++) {
	if(numParent[i]==0) {
	  ilike+=log(nodes[i].getPrior(state[n][i]));
	}   
	else {
	  //get parent state
	  int ps =0;
	  for(int j=0; j<numParent[i]; j++) {
	      ps=ps*nDiscrete[pList[i][j]]+state[n][pList[i][j]];
	  }
	  ilike+=log(nodes[i].getTable(ps, state[n][i]));
	}
      }
      cout<<n<<" "<<ilike<<endl;
    }
}

/* given a state, calculate the probability of the state */
void BN::inference(string file, int iDim, string outfile) {
    int maxlen = 819200;
    char line[maxlen];
    short** state;
    
    NEWP(state, iDim, short);
    for(int i=0; i<iDim; i++) {
      NEW(state[i], size, short);
    }

    //posterior
    double plike[iDim];
    double pae[size][3];

    bool sig[3][size];

    //read in data
    cout<<"loading data...\n";
    ifstream dfile(file.c_str()); 
    if(!dfile) {
        cout <<"Can not open data file "<<file<<"\n";
	exit(-1);
    }
    else {
	//read in data
        int nline=0;
	while(dfile.getline(line,maxlen,'\n') ) {
	    // comment line starting with #, so skip it
	    if(line[0] != '#') {
	        istrstream istr(line);
		int sValue;
		for(int i=0; i<size; i++) {
		  istr>>sValue;
		  //cout<<nline<<" "<<i <<" "<<sValue<<endl;
		  state[nline][i] = sValue;
		}
		nline++;
		if (nline%100==0) {
		  cerr<<"\015"<<nline;
		}
	    }
	}
	cerr<<endl;
	dfile.close();   
    }

    //posterior
    cout<<"calculating posteriors...\n";

    //get base
    float blike =0;
    int n=0;
    for(int i=0; i<size; i++) {
      if(numParent[i]==0) {
	blike+=log(nodes[i].getPrior(state[n][i]));
      }   
      else {
	//get parent state
	int ps =0;
	for(int j=0; j<numParent[i]; j++) {
	    ps=ps*nDiscrete[pList[i][j]]+state[n][pList[i][j]];
	}
	blike+=log(nodes[i].getTable(ps, state[n][i]));
      }
    }
    cout<<"baseLike "<<blike<<endl;
    //go through data and calculate the joint probabilities
    cout<<"calculating likelihood...\n";
    ofstream out(outfile.c_str(),ios::out);
    for(int n=0; n<iDim; n++) {
      if(n%100==0) {
	cerr<<"\015"<<n;
      }
      float ilike = 0;
      for(int i=0; i<size; i++) {
	if(numParent[i]==0) {
	  ilike+=log(nodes[i].getPrior(state[n][i]));
	}   
	else {
	  //get parent state
	  int ps =0;
	  for(int j=0; j<numParent[i]; j++) {
	      ps=ps*nDiscrete[pList[i][j]]+state[n][pList[i][j]];
	  }
	  ilike+=log(nodes[i].getTable(ps, state[n][i]));
	}
      }
      out<<ilike<<endl;
      plike[n] = exp(ilike-blike);
      
    }
    cerr<<endl;
    out.close();
    return;

    //calculate the posterior and output the result
    for (int i=0; i<size; i++) {
      cerr<<"\015"<<i;
      //initialization
      for (int k=0; k<size; k++) {
	  for (int s=0; s<3; s++) {
	    sig[s][k] =0;
	  } 
	}
      for (int j=0; j<3; j++) {
	int ne=0;
	double pe =0;
	for (int k=0; k<size; k++) {
	  for (int s=0; s<3; s++) {
	    pae[k][s] =0;
	  } 
	}

	for (int n=0; n<iDim; n++) {
	  if(state[n][i] ==j) {//perturbation (evidence)
	    pe += plike[n];
	    ne++;
	    for (int k=0; k<size; k++) {
	      pae[k][state[n][k]] +=plike[n];
	    } 
	  }
	}
	for (int k=0; k<size; k++) {
	  for (int s=0; s<3; s++) {
	    pae[k][s] /=pe;
	  } 
	}
	//select perturbation signature
	double std ;
	double cutoff = 4.6;
	for (int k=0; k<size; k++) {
	  for (int s=0; s<3; s+=2) {
	    std = sqrt(ne*nodes[k].getPrior(s)*(1-nodes[k].getPrior(s)));
	    if((pae[k][s]-nodes[k].getPrior(s))/std >cutoff) {
	      sig[j][k] = 1;
	    } 
	  }
	}
      }
      //summarize perturbation signature
      for (int k=0; k<size; k++) {
	if((sig[0][k]||sig[2][k]) && !sig[1][k]) {
	  cout<<1<<"\t";
	} 
	else {
	  cout<<0<<"\t";
	}
      }
      cout<<endl;
    }
    out.close();
}




/* given a node, and it state, simulate 10,000 states with the constraint of the given nodes, then predict the outcome */
void BN::predictOutcome(int maxN) {
    int  state[size];
    bool simulated[size];
    int  count[size][3];  //state count
    int  N=0;

    //initialization
    for(int i=0; i<size; i++) {
      for (int j=0; j<3; j++) {
	count[i][j]=0;
      }
    }

    while(N<maxN) {
      cerr<<N<<endl;
      //initialization
      for(int i=0; i<size; i++) {
	simulated[i] = false;
      }
      
      //simulation
      bool done =false;
      while(!done) {
	done =true;
	for (int i=0; i<size; i++) {
	  if(!simulated[i]) {
	    done=false;
	    if (numParent[i] ==0) {
	      //no parent
	      float m = randomf(1.0); 
	      float msum=0;
	      for (int j=0; j<nDiscrete[i]; j++) {
		if(m>msum && m<msum+ nodes[i].getPrior(j)) {
		  state[i]=j;
		}     
		msum += nodes[i].getPrior(j);
	      }
	      simulated[i] = true;
            }
	    else { //with parents
	      //get parent state
	      int ps =0;
	      for(int j=0; ps>=0 &j<size; j++) {
		if(relations[j][i]) {
		  ps=ps*nDiscrete[j]+state[j];
		  if(!simulated[j]) {
		    ps = -100000;
		  }
		}
	      }
	      
	      if (ps >=0) {
		float m = randomf(1.0); 
		float msum=0;
		for (int j=0; j<nDiscrete[i]; j++) {
		  if(m>msum && m<msum+ nodes[i].getTable(ps, j)) {
		    state[i]=j;
		  }     
		  msum += nodes[i].getTable(ps,j);
		}
		simulated[i] = true;
	      }
	    }
	  }
	}
      }
      //check input nodes
      bool accepted =true;
      /*for(int i=0; i<changeNodes.size(); i++) {
	int id= getNodeIdByName(changeNodes[i]);
	if (state[id] != changeStates[i]) {
	  accepted = false;
	}
      }*/
      if(accepted) {
	N++;
	/*for(int i=0; i<size; i++) {
	  count[i][state[i]]++;
	  }*/
	for(int i=0; i<size; i++) {
	  cout<<state[i]<<" ";
	}
	cout<<endl;
      }
    }

    /*output result
    for(int i=0; i<size; i++) {
      cout<<nodes[i].getName()<<" ";
      for (int j=0; j<nDiscrete[i]; j++) {
	cout<<count[i][j]<<" "<<nodes[i].getPrior(j)<<" ";
      }
      cout<<endl;
      }*/
}

/* pick a state based on prior */
int BN::pickState(int i) {
  float m = randomf(1.0); 
  float msum=0;
  int s;
  for (int j=0; j<nDiscrete[i]; j++) {
    if(m>msum && m<msum+ nodes[i].getPrior(j)) {
      s = j;
    }     
    msum += nodes[i].getPrior(j);
  }
  return s;
}

/* pick a state based on parents state */
int BN::simulateState(int i, bool* simulated, int* state) {
  int s;
  //get parent state
  int ps =0;
  for(int j=0; j<size; j++) {
    if(relations[j][i]) {
      if(simulated[j] ){
	ps=ps*nDiscrete[j]+state[j];
      }
      else {
	ps=ps*nDiscrete[j]+pickState(j);
      }
    }
  }	      
  float m = randomf(1.0); 
  float msum=0;
  for (int j=0; j<nDiscrete[i]; j++) {
    if(m>msum && m<msum+ nodes[i].getTable(ps, j)) {
      s=j;
    }     
    msum += nodes[i].getTable(ps,j);
  }
  return s;
}

/* give a child state , inference the parant state */
int BN::inferenceState(int pi, int ci, int cs) {
    int ps;   //parent state
    
    //initialize table
    float* table;
    int* count;
    NEW(table, nDiscrete[pi], float);
    NEW(count, nDiscrete[pi], int);

    for(int k=0; k<nDiscrete[pi]; k++) {
      count[k]=0;
    }
	
    //calculate the MLE of conditional probability
    for(int d=0; d<dataDim; d++) {
      if((int) data[ci][d] ==cs) {
	count[(int)data[pi][d]]++;
      }
    }

    //normalize
    float sum=0;
    for(int k=0; k<nDiscrete[pi]; k++) {
      sum+=count[k];
    }

    for(int k=0; k<nDiscrete[pi]; k++) {
      table[k]=count[k]/sum;
    }

    //inference
    float m = randomf(1.0); 
    float msum=0;
    for (int j=0; j<nDiscrete[pi]; j++) {
      if(m>msum && m<msum+ table[j]) {
	ps = j;
      }     
      msum += table[j];
    }
    return ps;
}


/* given a node, and it state, simulate states with the constraint of the given nodes, then predict the outcome */
void BN::predictOutcome(vector<string> changeNodes, vector<int> changeStates,
			string pathfile, int maxN) {
    int  count[size][3];
    int  state[size];
    bool simulated[size];
    int  N=0;
    int** path;

    //initializing counts
    for (int i=0; i<size; i++) {
      for (int j=0; j<3; j++) {
	count[i][j] =0;
      }
    }

    //get in path 
    cerr<<"reading pathfile ...";
    path = readPath(pathfile);
    cerr<<"done"<<endl;
    int cid = getNodeIdByName(changeNodes[0]);
    //simulation 
    while(N<maxN) {
      permute();
      cerr<<N<<endl;
      //initialization
      for(int i=0; i<size; i++) {
	if(path[cid][i] ==0) {
	  state[i]=pickState(i);
	  simulated[i] = true;
	}
	else { 
	  simulated[i] = false;
	}
      }
      //input
      state[cid] = changeStates[0];

      //simulation
      bool done =false;
      int D=0;
      while(!done) {
	D++;
	cerr<<N<<": dist "<<D<<endl;
	done =true;
	for (int I=0; I<size; I++) { 
	  int i= permIdx[I];
	  if(!simulated[i]) {
	    done=false;
	    for (int J=0; !simulated[i] & J<size; J++) {
	      int j = permIdx[J];
	      if(relations[j][i] & simulated[j]) {
		//based on parents states
		state[i] =simulateState(i, simulated, state);
		simulated[i] = true;
	      }
	      else if(relations[i][j] & simulated[j]) {
		//based on child states
		state[i] = inferenceState(i,j, state[j]);
		simulated[i] = true;
	      }
	    }
	  }
	}
      }
      //count
      N++;
      for(int i=0; i<size; i++) {
	cout<<state[i]<<"\t";
	count[i][state[i]]++;
      }
      cout<<endl;
    }

    /*
    //standard deviation of count
    float avg = 0.333*maxN;
    float std = sqrt(maxN*0.3333*0.6667);
    float cutoff= 4.3*std;

    for(int i=0; i<size; i++) {
      if(count[i][0]-avg>0) {
	cout<<nodes[i].getName()<<" "<<0<<endl;
      }
      if(count[i][2]-avg>cutoff) {
	cout<<nodes[i].getName()<<" "<<2<<endl;
      }
    } 
    */

}

/* given a node, and it state, calculate other genes */
void BN::predictChange(vector<string> changeNodes, vector<int> changeStates) {

    int  state[size];
    bool fixedNodes[size];
    bool checked[size];

    //initialization
    for(int i=0; i<size; i++) {
      //set initial state
      state[i] =1;
      fixedNodes[i] = false;
      checked[i] = false;
    }

    //set input nodes
    for(int i=0; i<changeNodes.size(); i++) {
      int id= getNodeIdByName(changeNodes[i]);
      fixedNodes[id] = true;
      state[id]= changeStates[i];
    }

    //iteratively predict changes
    vector<int> waitingList;
    
    //set waiting list
    for(int pid=0; pid<size; pid++) {
      if(fixedNodes[pid]) {
	for(int i=0; i<size; i++) {
	  if(!fixedNodes[i] && relations[pid][i]) {
	    waitingList.push_back(i);
	  }
	}
      }
    }

    while(waitingList.size()>0) {
      #ifdef DEBUG
      cout<<"waiting ";
      for(int i=0; i<waitingList.size(); i++) {
	cout<<nodes[waitingList[i]].getName()<<" ";
      }
      cout<<endl;
      #endif
      int i= waitingList[0];
      checked[i] = true;
      waitingList.erase(waitingList.begin());
      #ifdef DEBUG
      cout<<"testing "<<nodes[i].getName()<<endl;
      #endif
      //get parent state
      int ps =0;
      for(int j=0; j<size; j++) {
	if(relations[j][i]) {
	  ps=ps*nDiscrete[j]+state[j];
	}
      }

      /* deterministic */
      float m= max(max(nodes[i].getTable(ps, 0),
		       nodes[i].getTable(ps, 1)),
		   nodes[i].getTable(ps, 2));
      int ms;
      if(m==nodes[i].getTable(ps, 1)) ms =1;
      else if(m==nodes[i].getTable(ps, 0)) ms =0;
      else if(m==nodes[i].getTable(ps, 2)) ms =2;

      /*
      if(ms==1) {//check trend
	if(nodes[i].getTable(ps, 0)>1/nodes[i].nNotions() &&
	   nodes[i].getTable(ps, 0)-nodes[i].getTable(ps, 2)>0.2) {
	  ms=0;
	}
	if(nodes[i].getTable(ps, 2)>1/nodes[i].nNotions() &&
	   nodes[i].getTable(ps, 2)-nodes[i].getTable(ps, 0)>0.2) {
	  ms=2;
	}
      }
      */


      if(state[i] != ms) {
	if(state[i] !=1 ) {
	  cerr<<"warning "<<nodes[i].getName()<<" \n"; 
	}
	state[i] =ms;
	//update waitinglist
	for(int j=0; j<size; j++) {
	  if(!fixedNodes[j] && relations[i][j]) {
	    waitingList.push_back(j);
	  }
	}
      }

      
      /* stochastic 
      float m = randomf(1.0); 
      int ms;
      float msum=0;
      for (int j=0; j<nDiscrete[i]-1; j++) {
	if(m>msum && m<msum+ nodes[i].getTable(ps, j) ) {
	  ms=j;
	}     
	msum += nodes[i].getTable(ps, j);
      }
      if(m>msum) {
	ms= nDiscrete[i]-1;
      }
      
      //update waitinglist
      for(int j=0; j<size; j++) {
	if(!fixedNodes[j] && relations[i][j]) {
	  if(!checked[j] || state[i] !=ms) {
	    waitingList.push_back(j);
	  }
	}
      }
      state[i] = ms;

      */
    }

    //output states
    for(int j=0; j<size; j++) {
      if(!fixedNodes[j] && state[j] !=1) {
	cout<<nodes[j].getName()<<" "<<state[j]<<endl;
      }
      //cout<<nodes[j].getName()<<" "<<state[j]<<endl;     
    }
}

/* 
 * calculate BIC value 
 * log(E(likelihood))-0.5*k*logN
 * where k is number of independent parameters for the network
 * N is number of data point
 *
 */
//double BN::BIC() {
  //return likelihood();

  // calculation parameter size
  //we assume all nodes have the same notion size.
  /*int n=nodes[0].nNotions();
  for(int i=0; i<size; i++) {
    int k = n-1;
    for(int j=0; j<numParent[i]; j++) {
      k = k*n;
    }
    bic += k;
  }
  */

  /* alpha actually is alpha*log(dataDim) to save time */
  //double bic = -alpha*numParams+likelihood();

  //return bic;
//}

/*
 * read in the shortest path file
 */
int** BN::readPath(string pathfile) {
  int** dist;
  NEWP(dist, size, int);
  for(int i=0; i<size; i++) {
    NEW(dist[i], size, int);
  }

  int maxlen = 81920;
  char* line;
  NEW(line, maxlen, char);

  //read in path file
  ifstream pfile(pathfile.c_str());
  if(!pfile) {
    cout <<"Can not open path file "<<pathfile<<"\n";
    exit(-1);
  }
  else {
    int i=0;
    //read in data
    while(pfile.getline(line,maxlen,'\n') ) {
      istrstream istr(line);
      for (int j=0; j<size; j++) {
	istr>>dist[i][j];
      }
      i++;
    }

    pfile.close();
  }
  return dist;
}

/*
 * fine all the shortest pathes using repeated squaring matrix multiplication
 */
int** BN::shortestPath() {
  int** dist;
  bool** path;
  NEWP(dist, size, int);
  NEWP(path, size, bool);
  for(int i=0; i<size; i++) {
    NEW(dist[i], size, int);
    NEW(path[i], size, bool);
  }

  //inialize
  for(int i=0; i<size; i++) {
    for(int j=0; j<size; j++) {
      if(relations[i][j]) {
	dist[i][j] = 1;
	path[i][j] = true;
	//cout<<nodes[i].getName()<<" "<<nodes[j].getName()<< " "<<dist0[i][j] <<endl;
      }
      else {
	dist[i][j] = size+1;
	path[i][j] = false;
      }
    }
  }

  bool change = true;
  for(int p=2; change; p++) {
    cerr<< "p="<<p<<endl;
    change = false;
    for(int i=0; i<size; i++) {
        //cout<<"\015"<<i;
	for(int j=0; j<size; j++) {
	  if(i!=j) {
	    //find a path
	    for(int k=0; k<size; k++) {
		if(path[i][k] && path[k][j] &&
                   dist[i][j]> dist[i][k]+dist[k][j]) {
		  dist[i][j] = dist[i][k]+dist[k][j];
		  path[i][j] = true;
		  //cout<<nodes[i].getName()<<" "<<nodes[j].getName()<< " "
		  //  <<dist1[i][j] <<endl;
		  change = true;
		}
	    }
	  }
	}
    }
  }

  //check whether there are cycles
  /*for(int i=0; i<size; i++) {
    for(int j=0; j<size; j++) {
      if(i!=j && dist[i][j]==1 && dist[j][i]<size) {
	cout<<"Error: cycle "<<nodes[i].getName()<<" "
	    <<nodes[j].getName()<<endl;
      }
    }
    }*/

  FREEP(path, size);
  return dist;
}

void BN::printShortestPath(int ** path) {
  int maxD = 1;
  int d=1;
  while(d<=maxD) {
    for(int i=0; i<size; i++) {
      for(int j=0; j<size; j++) {
	if(path[i][j]<size && path[i][j]>maxD) {
	  maxD = path[i][j];
	}
	if(path[i][j] == d) {
	  cout<<nodes[i].getName()<<" "<<nodes[j].getName()<< " "
	      <<path[i][j] <<endl;
	}
      }
    }
    d++;
  }
}

/*
 * find the Markov blanket of a node
 * Markev blanet consists of direct parents, direct children, and direct parent
 * of direct children
 */
void BN::findMarkovBlanket(int seed) {
  //get direct parent and direct children
  for(int i=0; i<size; i++) {
    if(relations[i][seed] || relations[seed][i]) {
      cout<<nodes[i].getName()<<endl;
    }
  }
  // get direct parent of direct children
  for(int i=0; i<size; i++) {
    if(relations[seed][i]) {
      for(int j=0; j<size; j++) {
	if(relations[j][i] && !relations[j][seed]&& j!=seed) {
	  cout<<nodes[j].getName()<<endl;
	}
      }
    }
  } 
}


/*
 * refine structure.
 * if a pair of linked nodes have commond parents, test whether they are
 * still related given the commond parents (conditional mutual information)
 */
void BN::refine(float threthold) {
  cout<<"refining...\n";
  bool change = true;
  while(change) {
    change = false;
    for(int p=0; p<size; p++) {
      Node* np = getNode(p);
      for(int i=0; i<size; i++) {
	if(relations[p][i]) {
	  Node* ni = getNode(i);
	  for(int j=0; j<size; j++) {
	    /*we want to testing whether there are redundant here for
	     *p-->i, p-->j, and i-->j
	     *we may delete either p-->j or i-->j
	     *if mi(i,j) >= mi(j,p), we test whether ni->cmi(j,p) is import
             *if mi(i,j) < mi(j,p), we test whether np->cmi(j,i) is import
	     */
	    if(relations[p][j] && i!=j && relations[i][j]) {
	      Node* nj = getNode(j);
	      if(false && nj->mutualInfo(ni)>= nj->mutualInfo(np)) {
		cout<<ni->getName()<<" "
		      <<np->getName()<<"->"<<nj->getName()<<" "
		      <<ni->cmi(np, nj)<<endl;
		if(ni->cmi(np, nj) <threthold) {
		  cout<<"remove "<<ni->getName()<<" "
		      <<np->getName()<<"->"<<nj->getName()<<" "
		      <<ni->cmi(np, nj)<<endl;
		  removeRelation(p, j);
		  change= true;
		}
	      }
	      else {
		cout<<np->getName()<<" "
		    <<ni->getName()<<"->"<<nj->getName()<<" "
		    <<np->cmi(ni, nj)<<endl;
		if(np->cmi(ni, nj) <threthold) {
		  cout<<"remove "<<np->getName()<<" "
		    <<ni->getName()<<"->"<<nj->getName()<<" "
		    <<np->cmi(ni, nj)<<endl;
		  removeRelation(i, j);
		  change= true;
		}
	      }
	    }
	  }
	} 
      }
    }
  }

  //remove cycles involving two nodes
  for(int i=0; i<size; i++) {
    for(int j=0; j<size; j++) {
      if(relations[i][j] && relations[j][i]) {
	  cout<<"cycle "<<nodes[i].getName()<<" "
	      <<nodes[j].getName()<<endl;
	  //there is a cycle
	  if(priors[i][j]<priors[j][i]) {
	    cout<<"removing "<<nodes[i].getName()<<"->"
		<<nodes[j].getName()<<" "<<priors[i][j]<<endl;
	    removeRelation(i, j);
	  }
	  else if(priors[i][j]>priors[j][i]){
	    cout<<"removing "<<nodes[j].getName()<<"->"
		<<nodes[i].getName()<<" "<<priors[j][i]<<endl;
	    removeRelation(j, i);
	  }
	  else {
	    cout <<"Warning: "<<nodes[j].getName()<<"<->"
		<<nodes[i].getName()<<" equal weight"<<endl;
	  }
      }
    }
  }

  //there may be cycle in the graph, break it depend on cmi
  change = true;
  cout<<"checking cycles....\n";
  while(change) {
    change = false;
    for(int p=0; p<size; p++) {
      Node* np = getNode(p);
      for(int i=0; i<size; i++) {
	if(relations[p][i]) {
	  Node* ni = getNode(i);
	  for(int j=0; j<size; j++) {
	    if(relations[i][j] && relations[j][p] && relations[p][i]) {
	      Node* nj = getNode(j);
	      float cp = np->cmi(ni, nj);
	      float ci = ni->cmi(np, nj);
	      float cj = nj->cmi(np, ni);
	      if(cp <ci && cp <cj) {
		cout<<"remove "<<np->getName()<<" "
		    <<ni->getName()<<"->"<<nj->getName()<<" "
		    <<cp<<endl;
		removeRelation(i, j);
		change= true;
	      }
	      if(ci <cp && ci <cj) {
		cout<<"remove "<<ni->getName()<<" "
		    <<nj->getName()<<"->"<<np->getName()<<" "
		    <<ci<<endl;
		removeRelation(j, p);
		change= true;
	      }
	      if(cj <cp && cj <ci) {
		cout<<"remove "<<nj->getName()<<" "
		    <<np->getName()<<"->"<<ni->getName()<<" "
		    <<cj<<endl;
		removeRelation(p, i);
		change= true;
	      }
	    }
	  }
	} 
      }
    }
  }
  
  //checking cycle involving more than three nodes
  change = true;

  while(change) {
    change = false;
    cout<<"calculating shortest paths....\n";
    int** path = shortestPath();
    for(int i=0; i<size; i++) {
      for(int j=0; j<size; j++) {
	 if(path[i][j]==1 && path[j][i]<size) {
	  cout<<"cycle "<<nodes[i].getName()<<" "
	      <<nodes[j].getName()<<endl;
	  //there is a cycle
	  Relationship r = weakestRelation(j, i, path);
	  if(priors[i][j]<r.getPrior()) {
	    cout<<"removing "<<nodes[i].getName()<<"->"
		<<nodes[j].getName()<<" "<<priors[i][j]<<endl;
	    removeRelation(i, j);
	  }
	  else {
	    cout<<"removing "<<nodes[r.getParent()].getName()<<"->"
		<<nodes[r.getChild()].getName()<<" "<<r.getPrior()<<endl;
	    removeRelation(r.getParent(), r.getChild());
	  }
	  change = true;

	  path=blockpath(path, i, j);
	 }
      }
    }
    //free old space first
    FREEP(path, size);
  }

  //if a node's in-degree is greater than MAXNUMPARENTS, then remove the 
  // weakest link
  for(int i=0; i<size; i++) {
    while(numParent[i]>MAXNUMPARENTS) {
	  Relationship r = weakestRelation(i);
	  cout<<nodes[i].getName()<<" "<<numParent[i]<<" "
	      <<"removing "<<nodes[r.getParent()].getName()<<"->"
	      <<nodes[r.getChild()].getName()<<endl;
	  removeRelation(r.getParent(), r.getChild());
    }
  }
}

//given a path and two nodes i and j, block anything i or j related.
int** BN::blockpath(int** path, int ni, int nj) {
  map <int, bool> t;

  //find all related
  for (int i=0; i<size; i++) {
    if(path[i][ni]<size || path[ni][i] <size ||
       path[i][nj]<size || path[nj][i] <size) {
      t[i] = true;
    }
  }

  //block all related
  for(map<int, bool>::iterator n= t.begin(); n!=t.end(); n++) {
    for(int i=0; i<size; i++) {
      path[i][n->first]=size; 
      path[n->first][i]=size;
    }
  }

  return path;

}

//give the starting node, ending node, and all-pair shortest path
//find the weakest link (relationship with lowest prior)
Relationship BN::weakestRelation(int ns, int ne, int** d) {
  //recursively find path
  int maxD = d[ns][ne];
  int pid, cid;
  float v=1;

  while(maxD>=2) {
    for(int i=0; i<size; i++) {
	if(d[ns][i]==1 &&
	   d[ns][i]+d[i][ne]<=maxD) {
	  if(priors[ns][i]<v) {
	    pid = ns;
	    cid = i;
	    v = priors[ns][i];
	  }
	  ns = i;
	  maxD--;
	  break;
	}
    }
  }
  //the last link
  if(priors[ns][ne]<v) {
    pid = ns;
    cid = ne;
    v = priors[ns][ne];
  } 

  //create a relationhip
  Relationship r = Relationship(pid, cid);
  r.setPrior(getPrior(pid, cid));

  return r;
}

//give a node, find the weakest link (relationship with lowest prior) to it
Relationship BN::weakestRelation(int c) {
  float v=10000;
  int pid=c; 
  for(int i=0; i<size; i++) {
    if(relations[i][c] && priors[i][c]<v) {
      pid = i;
      v = priors[i][c];
    }
  }

  //create a relationhip
  Relationship r = Relationship(pid, c);
  r.setPrior(getPrior(pid, c));

  return r;
}

/*
 * initialize the structure based on reachable information
 * 
 * If reachable[i][j] >1, and potentialRelations[i][j] is true
 * add the relations i->j.
 */
void BN::initStructure_reachable() {
  for(int i=0; i<size; i++) {
    for(int j=0; j<size; j++) { 
      if(potentialRelations[i][j] && reachable[i][j]>1) {
	addRelation(i,j);
	if(hasCycle(i)) {
	  removeRelation(i,j);
	}
      }
    }
  }
}

/*
 * The learning stratagy is local optimal search
 * 
 * the acceptance criteria is BIC
 */
void BN::learnStructure() {
  int i,j;
  bool** Relations;         // a structure for intermedium result
  NEWP(Relations, size, bool);
  for (i=0; i<size; i++) {
    NEW(Relations[i], size, bool);
  }
  for(i=0; i<size; i++) {
    for(j=0; j<size; j++) {
      Relations[i][j] = false;
    }
  }

  //add fixed relations
  for(i=0; i<size; i++) {
    for(j=0; j<size; j++) {
      if(fixedRelations[i][j]) cout<<"fixedR "<<i<<" "<<j<<endl;
      if(relations[i][j]) cout<<"relation "<<i<<" "<<j<<endl;
    }
  }

  //ancestor matrix
  bool** A;        
  NEWP(A, size, bool);
  for (i=0; i<size; i++) {
    NEW(A[i], size, bool);
  }
  Graph::init_ancestor_matrix(A, relations, size);

  //for intermedium result
  double maxL = BIC()-999999;

  //add relations based on reachable information
  //initStructure_reachable();

  double liken=BIC();
  double likep;
  cout<<"Baseline likelihood: "<<liken<<endl;
  int nUnchange =0;
  int maxIter = size*15;
  //int maxIter = size*7;
  if(maxIter<1500) maxIter = 1500;
  if(maxIter>40000) maxIter = size*10;
  //float Burnin = maxIter-1.5*size;
  float Burnin = 3*size;

  float T=1; // Gibbs sampling
  //float T=-1;  // greedy sampling
  if(T>0) permute();
  for(i=0; i<maxIter && nUnchange <2.1*size; i++) {
    int n=randomi(size*10)%size;
    if (nUnchange>=size) {
      n=i%size;
    }
    
    cout<<"Iteration #"<<i<<" Node #"<<n<<"  "<<numPotentialParent[n]<<endl;
    int   naction=-1;
    liken = BIC(n)-5;
    //test each potential parent
    //for(int pp=0; pp<size; pp++) {
    for (int npp=0; npp<size; npp++) {
      int pp= permIdx[npp];
      //the potential parent is not a child 
      if(potentialRelations[pp][n] && !relations[n][pp]) {
	if(relations[pp][n]) {
	  if(!fixedRelations[pp][n]) {
	    //choice #1 reverse
	    if(numParent[pp]<MAXNUMPARENTS &&
	       potentialRelations[n][pp]) {
	      //make sure reachable nodes in Markov blanket are reachable
	      bool possible = true;
	      /*
	      for(int j=0; possible && j<size; j++) {
		if(reachable[j][n]>1 && 
		   relations[j][pp] && !relations[j][n]) {
		  possible = false;
		}
	      }
	      */
	      likep = BIC(pp);
	      if(possible) {
		reverseRelation(pp, n);		
		//printStructure();
		if(!hasCycle(n)) {
		  //cout<<"reversing parent "<<n<<" "<<pp<<" "
		  //    <<BIC()<<" "<<liken<<endl;
		  if(MC(BIC(n)+BIC(pp),liken+likep, -1)) {
		    //cout<<"reversing parent "<<n<<" "<<pp<<endl;
		    liken = BIC(n)+BIC(pp)-likep;
		    naction = pp*3;
		  }
		}
		else {
		  //cout<<"has cycle.\n";
		}
		//recover
		reverseRelation(n, pp);
	      }
	    }
	    //choice #2
	    //removing parent
	    removeRelation(pp, n);
	    //cout<<"removing parent "<<pp<<" from "<<n<<" "
	    //	<<BIC() <<" "<<liken<<endl;
	    //printStructure();
	    if(MC(BIC(n),liken, -1)) {
	      liken = BIC(n);
	      naction = pp*3+1;
	      //cout<<"removing parent "<<pp<<" from "<<n <<endl;
	    }
	  
	    //recover
	    addRelation(pp, n);	  
	  }
	}
	else {//add it to parant list if it is not on the list now
	    //choice #3
	  if(numParent[n] < MAXNUMPARENTS) {
	    //make sure reachable nodes in Markov blanket are reachable
	    bool possible = true;
	    /*
	    for(int j=0; possible && j<size; j++) {
	      if(reachable[j][pp]>1 && relations[j][n] &&
		 !relations[j][pp]) {
		possible = false;
	      }
	    }
	    */
	    if(possible) {
	        addRelation(pp,n);
		/*if(i>0) {
		  if (hasCycle(pp)) {
		    if(!A[n][pp]) {
		      cout<<"cycle "<<pp<<" "<<n<<endl;
		      cout<<"A "<<A[n][pp]<<endl;
		      exit(0);
		    }
		  }
		  else if(A[n][pp]) {
		    cout<<"no cycle "<<pp<<" "<<n<<endl;
		    cout<<"A "<<A[n][pp]<<endl;
		    exit(0);
		  }
		  
		}*/
		//printStructure();
		if(!A[n][pp]) {
		    //cout<<"adding parent "<<pp<<" to "<<n<<" "
		    //  	<<BIC(n)<<" "<<liken<<endl;
		    //cout<<getNode(pp)->getName()<<" "<<BIC(n)<<endl;
		    if(MC(BIC(n), liken, -1)) {
		        //cout<<"adding parent "<<pp<<" to "<<n<<endl;
		        liken = BIC(n);
			naction = pp*3+2;
		    }
		}
		//cout<<"recovering...\n";
		//recover
		removeRelation(pp, n);
	    }
	  }
	}
      }
    }
    //update structure
    if(naction>=0 && MC(liken, BIC(n), T)) { //likelihood improves
      int pp = naction/3;
      cout<<"action: "<<pp<<" "<<n<<" "<<naction%3<<endl;
      switch (naction%3) { 
      case 0:
	  //reverseRelation(pp, n);
	  removeRelation(pp,n);
	  Graph::do_removal(A, pp, n, relations, size);
	  addRelation(n,pp);
	  Graph::do_addition(A, n, pp, size);
	  break;
      case 1: 
	  removeRelation(pp, n);
	  Graph::do_removal(A, pp, n, relations, size);
	  break;
      case 2:
	  addRelation(pp, n);
	  Graph::do_addition(A, pp, n, size);
	  break;
      }
      liken=BIC();
      //save the maximum structure
      if(i>Burnin & liken>maxL) {
	maxL=liken;
	//copy structure
	for(int ii=0; ii<size; ii++) {
	  for (j=0; j<size; j++) {
	    Relations[ii][j] = relations[ii][j];
	  }
	}
	//printStructure();
	//likelihood();
      }
      nUnchange = 0;
      cout<<"Likelihood="<<liken<<endl;
    }
    nUnchange++;
    
    //if(i==4489) break;
  }

  // in case exist too quickly
  //save the maximum structure
  if(i<Burnin & liken>maxL) {
    maxL=liken;
    //copy structure
    for(int ii=0; ii<size; ii++) {
      for (j=0; j<size; j++) {
	Relations[ii][j] = relations[ii][j];
      }
    }
  }

  //clear the structure
  for(i=0; i<size; i++) {
    for (j=0; j<size; j++) {
      if(relations[i][j]) removeRelation(i,j);
    }
  }

  //set the structure to the best structure
  for(i=0; i<size; i++) {
    for (j=0; j<size; j++) {
      if(Relations[i][j]) addRelation(i,j);
    }
  }
  
  //free space
  FREEP(Relations, size);
  FREEP(A, size);

  cout<<"nTrylist "<<trylist.size()<<endl;
}

/*
 * enumerate possible configurations
 * then calculate scores
 */
void BN::enumerateStructures() {
  //no parent
  for(int i=0; i<size; i++) {
    BIC(i);
  }

  //one parent
  for(int i=0; i<size; i++) { 
    for (int p=0; p<size; p++) {
      if(potentialRelations[p][i]) {
	addRelation(p,i);
	BIC(i);
	removeRelation(p, i);
      }
    }
  }

  cout<<"nTrylist "<<trylist.size()<<endl;
}

/*
 * enumerate possible configurations
 * then calculate scores
 */
void BN::enumerateStructures(string enumerationFile) {

  ifstream tfile(enumerationFile.c_str());
  if(!tfile) {
    cout <<"Can not open file "<<enumerationFile<<"\n";
    exit(-1);
  }
  else {
    int maxlen=1024;
    char line[maxlen];
    int  pa[MAXNUMPARENTS];
    int  np=0;
    int  cid=0;
    //read in trylist
    while(tfile.getline(line,maxlen,'\n') ) {
      //child node
      sscanf(line, "%d", &cid);

      //count number of parents
      np=0;
      for (int i=0; i<strlen(line); i++) {
	if (line[i] == ':') {
	  sscanf(&line[i+1], "%d",&pa[np]);
	  addRelation(pa[np], cid);
	  np++;
	}
      }
      BIC(cid);

      //reset structure
      for (int i=0; i<np; i++) {
	removeRelation(pa[i], cid);
      }
    }
    tfile.close();
  }
}

void BN::printStructure() {
  cout<<"Structure:"<<endl;
  for(int i=0; i<nodes.size(); i++) {
    for(int j=0; j<nodes.size(); j++) {
      if(relations[i][j]) {
	cout<<nodes[i].getName()<<"->"
	    <<nodes[j].getName()<<endl;
      }
    }
  }
  cout<<"Likelihood="<<BIC()<<endl;
}

void BN::toDotInput(string file) {
  ofstream out(file.c_str(),ios::out);
  //header
  out<<"digraph G {"<<endl;
  for(int i=0; i<nodes.size(); i++) {
    for(int j=0; j<nodes.size(); j++) {
      if(relations[i][j]) {
	out<<nodes[i].getName()<<"->"
	    <<nodes[j].getName()<<endl;
      }
    }
  }

  for(int i=0; i<size; i++) {
      if(numParent[i] ==0 && numChildren[i] ==0) {
	  out<<nodes[i].getName()<<";"<<endl;
      }
  }
  out<<"}"<<endl;
  out.close();
}


/*check whether there are cycle for a node */
bool BN::hasCycle(int node) {
  bool cycleF=false;
  map<int, int> waitList;
  map<int, int> checkedList;


  for(int i=0; i<size; i++) {
    if(relations[node][i]) {
      waitList[i]=1;
    }
  }

  while(!cycleF && waitList.size()>0) {
    int nodetmps = waitList.begin()->first;
    waitList.erase(waitList.begin());
    checkedList[nodetmps]=1;

    for(int i=0; i<size; i++) {
      if(relations[nodetmps][i]) {
	if(i==node) {
	  cycleF=true;
	}
	else {
	  if(checkedList.find(i)==checkedList.end() ) {
	    //add to waiting list
	    waitList[i]=1;
	  }    
	}
      }
    }
  }
  return cycleF;
}

/* simulated annealing 
boolean BN:SA(float E1, float E2, int T) {
  if(E1>E2) {
    return true;
  }
  else {
  }
}
*/

void BN::toXMLBIF(string xmlfile) {
  ofstream out(xmlfile.c_str(),ios::out);
  //write header
  out<<"<?xml version=\"1.0\"?>"<<endl;
  out<<endl;
  out<<"<!-- DTD for the BIF format -->"<<endl;
  out<<"<!DOCTYPE BIF ["<<endl;
  out<<"      <!ELEMENT BIF ( NETWORK )*>"<<endl;
  out<<"      <!ELEMENT PROPERTY (#PCDATA)>"<<endl;
  out<<"      <!ELEMENT TYPE (#PCDATA)>"<<endl;
  out<<"      <!ELEMENT VALUE (#PCDATA)>"<<endl;
  out<<"      <!ELEMENT NAME (#PCDATA)>"<<endl;
  out<<"      <!ELEMENT NETWORK"<<endl;
  out<<"          ( NAME, ( PROPERTY | VARIABLE | PROBABILITY | LIKELIHOOD)* )>"<<endl;
  out<<"      <!ELEMENT VARIABLE ( NAME, TYPE, ( VALUE |  PROPERTY )* ) >"<<endl;
  out<<"      <!ELEMENT PRIOR ( FOR | GIVEN | PROPERTY )* >"<<endl;
  out<<"      <!ELEMENT PROBABILITY"<<endl;
  out<<"          ( FOR | GIVEN | TABLE | ENTRY | DEFAULT | PROPERTY )* >"<<endl;
  out<<"      <!ELEMENT FOR (#PCDATA)>"<<endl;
  out<<"      <!ELEMENT GIVEN (#PCDATA)>"<<endl;
  out<<"      <!ELEMENT TABLE (#PCDATA)>"<<endl;
  out<<"      <!ELEMENT DEFAULT (TABLE)>"<<endl;
  out<<"      <!ELEMENT ENTRY ( VALUE* , TABLE )>"<<endl;
  out<<"      <!ELEMENT LIKELIHOOD (#PCDATA)>"<<endl;
  out<<"]>"<<endl;

  //network
  out<<"<BIF>"<<endl;
  out<<"<NETWORK>"<<endl;
  out<<"<NAME>"<<name<<"</NAME>"<<endl;
  out<<"<LIKELIHOOD>"<<BIC()<<"</LIKELIHOOD>"<<endl;
  //variables
  out<<"<!-- Variables -->"<<endl;

  for(int i=0; i<size; i++) {
    out<<"<VARIABLE>"<<endl;
    out<<"    <NAME>"<<nodes[i].getName()<<"</NAME>"<<endl;
    out<<"    <TYPE>discrete</TYPE>"<<endl;
    for(int j=0; j<nDiscrete[i]; j++) {
      out<<"    <VALUE>"<<nodes[i].getNotion(j)<<"</VALUE>"<<endl;
    }

    out<<"</VARIABLE>"<<endl;    
  }

  //Probability distributions
  out<<"<!-- Probability distributions -->"<<endl;
  for(int i=0; i<size; i++) {
    out<<"<PROBABILITY>"<<endl;
    out<<"    <FOR>"<<nodes[i].getName()<<"</FOR>"<<endl;
    for(int j=0; j<size; j++) {
        if(relations[j][i])
	  out<<"    <GIVEN>"<<nodes[j].getName()<<"</GIVEN>"<<endl;
    }
    //output table
    if(numParent[i]>0) {
      //find out table dimension
      int np3=1;
      for(int j=0; j<size; j++) {
	if(relations[j][i]) 
        np3=np3*nDiscrete[j];
      }      
      for(int np=0; np<np3; np++) {
	out<<"    <TABLE>";
	for(int n=0; n<nDiscrete[i]; n++) {
	  out<<nodes[i].getTable(np, n)<<" ";
	}
	out<<"</TABLE>"<<endl;
      }
    }
    else {
      out<<"    <TABLE>";
      for(int n=0; n<nDiscrete[i]; n++) {
	out<<nodes[i].getPrior(n)<<" ";
      }
      out<<"</TABLE>"<<endl;
    }
    out<<"</PROBABILITY>"<<endl;    
  }  
  /*  
  //Prior distributions
  out<<"<!-- Prior distributions -->"<<endl;
  for(int i=0; i<nodes.size(); i++) {
    for(int j=0; j<nodes.size(); j++) {
      if(i!=j && priors[i][j] != -1) {
	out<<"<PRIOR>"<<endl;
	out<<"    <FOR>"<<nodes[j].getName()<<"</FOR>"<<endl;
	out<<"    <GIVEN>"<<nodes[i].getName()<<"</GIVEN>"<<endl;
	out<<"    <PROPERTY>"<<priors[i][j]<<"</PROPERTY>"<<endl;
	out<<"</PRIOR>"<<endl;
      }
    }
  }
  */
  //network
  out<<"</NETWORK>"<<endl;
  out<<"</BIF>"<<endl;
  out.close();
}

void BN::toBIF(string BIFfile) {
  cout<<"BIF output\n";
  ofstream out(BIFfile.c_str(),ios::out);
  //write header

  //network
  out<<"network "<<name<<" {"<<endl;
  out<<"}"<<endl;

  //variables
  for(int i=0; i<size; i++) {
    out<<"variable "<<nodes[i].getName()<<"{"<<endl;
    out<<"    type discrete["<<nDiscrete[i]<<"] {" ;
    for(int j=0; j<nDiscrete[i]; j++) {
      out<<nodes[i].getNotion(j)<<" ";
    }
    out<<"};\n";
    out<<"}\n";
  }

  //Probability distributions
  for(int i=0; i<size; i++) {
    out<<"probability ( "<<nodes[i].getName();
    if(numParent[i]>0) {
      out<<"| ";
    }
    for(int j=0; j<size; j++) {
        if(relations[j][i])
	  out<<"  "<<nodes[j].getName();
    }
    out<<") {\n";
    //output table
    out<<"  table ";
    if(numParent[i]>0) {
      //find out table dimension
      int np3=1;
      for(int j=0; j<size; j++) {
	if(relations[j][i]) 
        np3=np3*nDiscrete[j];
      }      
      for(int n=0; n<nDiscrete[i]; n++) {
	for(int np=0; np<np3; np++) {
	  out<<nodes[i].getTable(np, n)<<" ";
	}
      }
      out<<";"<<endl;
    }
    else {
      for(int n=0; n<nDiscrete[i]; n++) {
	out<<" "<<nodes[i].getPrior(n);
      }
      out<<";"<<endl;
    }
    out<<"}"<<endl;    
  }  
  out.close();
}

/* create script to create a network in matlab BNT package */ 
void BN::toBNT(string BNTfile) {
  //ofstream out(BNTfile.c_str(),ios::out,0644);
  ofstream out(BNTfile.c_str(),ios::out);
  //write header
  
  //variables
  out<<"node= struct (";
  for(int i=0; i<size-1; i++) {
    out<<"'"<<nodes[i].getName()<<"',"<<i+1<< ", ..."<<endl;
  }
  out<<"'"<<nodes[size-1].getName()<<"',"<<size<<");"<<endl;

  out<<endl;
  // adjacency table
  out<<"adjacency = zeros("<<size<<");"<<endl;
  for(int i=0; i<size; i++) {
    if(numParent[i]>0) {
      out<<"adjacency([";
      for(int j=0; j<size; j++) {
        if(relations[j][i])
	  out<<j+1<<" ";
      }
      out<<"]," <<i+1<<")=1;"<<endl;
    }
  }
  out<<endl;

 
  out<<"nDiscrete = [";
  for(int i=0; i<size; i++) {
    out<<nDiscrete[i]<<" ";
  }
  out<<"];"<<endl;

  out<<"bnet=mk_bnet(adjacency, nDiscrete);"<<endl;

  //Probability distributions
  for(int i=0; i<size; i++) {

    out<<"bnet.CPD{"<<i+1<<"}= tabular_CPD(bnet, "<<i+1<<",[";
    //output table
    if(numParent[i]>0) {
      //find out table dimension
      int np3=1;
      for(int j=0; j<size; j++) {
	if(relations[j][i]) 
        np3=np3*nDiscrete[j];
      }      
      for(int n=0; n<nDiscrete[i]; n++) {
	for(int np=0; np<np3; np++) {
	  out<<nodes[i].getTable(np, n)<<" ";
	}
      }
    }
    else {
      for(int n=0; n<nDiscrete[i]; n++) {
	out<<" "<<nodes[i].getPrior(n);
      }
    }
    out<<"]);"<<endl;    
  }  
  out.close();
}
