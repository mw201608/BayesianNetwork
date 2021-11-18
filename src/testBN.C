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

/**
 *  construct a Bayesian network given a set of data 
 */

#include <iostream> 
#include <fstream>
#include <strstream>
#include <math.h>
#include <vector>
#ifdef LINUX
    #include <getopt.h>
#endif
#include "common.h"
#include "BN.h"
using namespace std;
/**
 *$Revision: 1.9 $
 *$Date: 2004/04/12 15:49:53 $
 *@Author:  Jun Zhu
 *@Creation Date: 3/23/01
 */

void usagemessage (char *name) {
    cout << "Usage: "<< name << " -s seed -b BIF file -d datafile [-e eQTL file] -t threthold_i -T threthold_e -D data dimension -o outfile [-x xmlfile] [-T timeSeriesFlag] [-c  convert file] [-p 0/1] [-i inference data -I inference dimension] [-r rescale_ratio] [-a alpha] [-C change node -S state -N times] [-P priorfile]\n";
    exit(0);
}

int main(int argc, char **argv) {
    bool   pseudoFlag = 1;
    string bannedListFile;  
    string biffile;  //input
    string BIFfile;  //output
    string BNTfile;  //BNT output
    string datafile;
    string eqtlfile;
    string priorfile;
    string outfile;
    string xmlfile;
    string inferencefile;
    string convertfile;
    string corrfile;
    string trylistFile;
    string updateTryListFile;
    string enumerationFile;
    float  threthold_i, threthold_e;
    int    timeSeriesFlag =0;
    int    listF=0;
    int    fitF=0;
    int    dim=0;
    int    iDim=0;
    int    seed=0;
    int    updateNodeId =-1;
    int    maxTrylist =-1;
    string initP="d"; //o: overlap; c:causal; n: number of QTL
    string pathfile ;
    string MarkovBlanketSeed ;
    string nodeSeed;
    vector <string> changeNodes;
    vector <int>    changeStates;
    int    nPrediction = 0;         
    float  ratio=1;                 //the scaling factor for the prior
    float  alpha =0.5 ;             //factor for BIC
    float  acutoff = 0.01;          //association cutoff
    float  qcutoff = 2.0;           //QTL cutoff
    float  qratio_cutoff =0.3;      //QTL ratio cutoff, used in direction determination
    bool**  bannedList;              //prohibited relationships

    int option_char;
    while ((option_char = getopt (argc, argv, "a:b:c:d:e:f:g:i:l:m:n:o:p:q:r:s:t:u:x:A:B:C:D:E:F:I:L:M:N:P:Q:R:S:T:U:X:")) != -1) {
        switch (option_char) {
	case 'a': alpha = atof(optarg);break;
	case 'b': biffile = optarg; break;
	case 'c': convertfile = optarg; break;
	case 'd': datafile = optarg; break;
	case 'e': eqtlfile = optarg; break;
	case 'f': pseudoFlag = atoi(optarg); break;
	case 'g': bannedListFile = optarg; break;
	case 'i': inferencefile = optarg; break;
        case 'l': trylistFile = optarg; break;
	case 'm': MarkovBlanketSeed = optarg; break;
	case 'n': nodeSeed = optarg; break;
	case 'o': outfile = optarg; break;
	case 'p': pathfile = optarg; break;
	case 'q': qratio_cutoff = atof(optarg); break;
	case 'r': ratio = atof(optarg);break;
	case 's': seed = atoi(optarg); break;
	case 't': threthold_i = atof(optarg); break;
        case 'u': updateNodeId = atoi(optarg); break;
	case 'x': xmlfile = optarg; break;
	case 'A': acutoff = atof(optarg); break;
	case 'B': BNTfile = optarg; break;
        case 'C': changeNodes.push_back(optarg); break;
	case 'D': dim = atoi(optarg); break;
	case 'E': enumerationFile = optarg; break;
	case 'F': fitF =atoi(optarg); break;
	case 'I': iDim = atoi(optarg); break;
	case 'L': listF = atoi(optarg);break; 
	case 'M': maxTrylist = atoi(optarg); break;
	case 'N': nPrediction = atoi(optarg);break;
	case 'P': priorfile = optarg;break; 
	case 'Q': qcutoff = atof(optarg); break;
	case 'R': corrfile = optarg; break;
        case 'S': changeStates.push_back(atoi(optarg)); break;
	case 'T': threthold_e = atof(optarg); break;
	case 'U': updateTryListFile = optarg; break;
	case 'X': initP = optarg; break;
	default: usagemessage(argv[0]);
	}
    }
    cout<<"X "<<initP<<endl;

    /* data dimension is not automatically checked, so we need to
     * know this ahead of time.
     */
    if(dim==0) {
      cout<<"Data dimension is required.\n";
      usagemessage(argv[0]);
    }
    
    if(biffile.length()==0) {
      cout<<"BIF file is required.\n";
      usagemessage(argv[0]);
    }

    /* at least one output file is required. */
    if(outfile.length()==0 && xmlfile.length() ==0) {
      cout<<"Normal outfile or xmlfile output is required.\n";
      usagemessage(argv[0]);
    }
        
    //initialize random seed
    cerr<<"seed="<<seed<<endl;
    if(seed==0) {
      initrandom();
    }
    else {
      srand(seed);
    }

    //create an empty Bayesian network;
    BN bn = BN();

    //set alpha for BIC
    bn.setAlpha(alpha);

    //set data dimension
    bn.setDataDim(dim);

    //set pseudoFlag
    bn.setPseudoFlag(pseudoFlag);

    //set acutoff
    bn.setACutoff(acutoff);

    //set qcutoff
    bn.setQCutoff(qcutoff);

    if(maxTrylist>0) {
      bn.setMaxTryList(maxTrylist);
    }

    //read in initial Bayesian network if exists
    if(biffile.length()>0) {
      cerr<<"reading "<<biffile<<endl;
        if(bn.init(biffile)!=0) exit(0);
    }

    //read correlation file
    if(corrfile.length()>0) {
      bn.readCorr(corrfile);
      cerr<<"correlation file "<<corrfile <<" is loaded."<<endl;
    }

    int maxlen = 819200;
    char* line;
    NEW(line, maxlen, char);

    //read in data
    ifstream dfile(datafile.c_str()); 
    if(!dfile) {
        cerr <<"Can not open data file "<<datafile<<"\n";
	exit(-1);
    }
    else {
	//read in data
	while(dfile.getline(line,maxlen,'\n') ) {
	    // comment line starting with #, so skip it
	    if(line[0] != '#') {
	        istrstream istr(line);
		string stmp;
		istr>>stmp;
		if(bn.getNodeByName(stmp) ==NULL) {
		    Node n = Node(stmp);
		    bn.addNode(n);
		}
		Node *nd= bn.getNodeByName(stmp);
		for(int i=0; i<dim; i++) {
		    float datav;
		    istr>>datav;
		    nd->addData(datav);
		}
		if(nd->isDiscreteNode()) {
		  nd->calculatePrior();
		}
	    }
	}
	dfile.close(); 
	//put all data into data matrix
	bn.readdata();
    }
    
    //convert grahp to directed graph 
    if(convertfile.length()>0) {
      cout <<"reading  "<<convertfile<<"....\n";
      ifstream cfile(convertfile.c_str()); 
      if(!cfile) {
        cout <<"Can not open convert file "<<convertfile<<"\n";
	exit(-1);
      }
      else {
	//read in data
	string n1, n2, label;
	string pat;
	while(cfile.getline(line, maxlen,'\n') ) {
	    string sl = string(line);
	    // skip header line
	    if(sl.find("--") != -1 || sl.find("->") != -1 ) {
	      if(sl.find("--") != -1 ) pat = "--";
	      else if (sl.find("->") != -1 ) pat = "->";
	      n1 = sl.substr(0, sl.find(pat));
	      if(sl.find("--") !=-1) {
		if(sl.find(' ')!=-1) {
		  n2 = sl.substr(sl.find("--")+2, 
				 sl.find(' ')-sl.find("--")-2);
		  label = sl.substr(sl.find(' '), sl.length());
		}
		else {
		  n2 = sl.substr(sl.find("--")+2, sl.length());
		  label = "";
		}
		//the relationship is two-way
		bn.addRelation(bn.getNodeIdByName(n2),
			       bn.getNodeIdByName(n1));
	      }
	      else {
		if(sl.find(" ")!=-1) {
		  n2 = sl.substr(sl.find("->")+2, 
				      sl.find(' ')-sl.find("->")-2);
		  label = sl.substr(sl.find('=')+1, sl.find(']')-sl.find('=')-1);
		}
		else {
		  n2 = sl.substr(sl.find("->")+2, sl.length());
		  label = "";
		}
	      }
	      #ifdef DEBUG
	         cout<<sl<<endl;
	         cout<<"n1="<<n1<<endl;
		 cout<<"n2="<<n2<<endl;
		 cout<<"label="<<label<<endl;
	      #endif
	      bn.addRelation(bn.getNodeIdByName(n1),
			     bn.getNodeIdByName(n2));
	      if(label.length()>0) {
		Relationship r = Relationship(bn.getNodeIdByName(n1),
					      bn.getNodeIdByName(n2));
		r.setPrior(atof(label.data()));
		bn.setPrior(&r);
	      }
	    }
	}
	cfile.close();   
      }
      

      //find the Markov blanket of a seed node
      if(MarkovBlanketSeed.length() >0) {
	bn.findMarkovBlanket(bn.getNodeIdByName(MarkovBlanketSeed));
	exit(0);
      }

      if(nodeSeed.length()>0) {
	bn.likelihood();
	bn.updateNodeTable(bn.getNodeIdByName(nodeSeed), true);
	bn.printNodeTable(bn.getNodeIdByName(nodeSeed));
	exit(0);
      }

      if(changeNodes.size()>0) {
	bn.likelihood();
	/*for(int i=0; i<nPrediction; i++) {
	  bn.predictChange(changeNodes, changeStates);
	  }*/
	bn.predictOutcome(changeNodes, changeStates, pathfile, nPrediction);
	exit(0);
      }

      if(changeNodes.size()==0 & nPrediction >0) {
	bn.likelihood();
	bn.predictOutcome(nPrediction);
	exit(0);
      }
      

      //calculate the probability of a specific state
      if(inferencefile.length()>0) {
	cout<<"calculate sample probability ...\n";
	bn.likelihood();
	bn.inference(inferencefile, iDim, outfile);
	exit(0);
      }

      //calculate shortest path
      if(pathfile.length()>0) {
	bn.printShortestPath(bn.shortestPath()); exit(0);
      }

      if(!fitF && eqtlfile.length()==0) {
        //refine graph
        bn.refine(threthold_i);
        //bn.printStructure();
        if(outfile.length()>0)  bn.toDotInput(outfile);
        exit(0);
      }
    }

    //read in eQTL information
    if(eqtlfile.length()>0) {
      ifstream efile(eqtlfile.c_str()); 
      if(!efile) {
        cout <<"Can not open eQTL file "<<eqtlfile<<"\n";
	exit(-1);
      }
      else {
	/* read in eqtl
	 * each line in the file is the following format
	 * Node_name   Locus    LODS
	 */
	while(efile.getline(line,maxlen,'\n') ) {
	    // comment line starting with #, so skip it
	    if(line[0] != '#') {
	        istrstream istr(line);
		string stmp;
		float loc, lod, R2;
		istr>>stmp>>loc>>lod>>R2;
		Node *nd= bn.getNodeByName(stmp);
		nd->getEQTL()->add(loc, lod, R2);
	    }
	}
	efile.close();   
	cout<<"eQTL file "<<eqtlfile<<" is loaded."<<endl;
      }
    }

    if(eqtlfile.length()>0) {
      if(initP =="c"){
	cerr<<"Prior by QTL causal\n"<<endl; 
	bn.initPrior_causal();
      }
      else if(initP =="n") {
	cerr<<"Prior by QTL complexity\n"<<endl; 
	bn.initPrior_complexity();
      }
      else if(initP =="o") {
	cerr<<"Prior by QTL match\n"<<endl; 
	bn.initPrior_qtlmatch(); 
      }
      else if(initP =="p"){
	cerr<<"Prior by QTL peak\n"<<endl; 
	bn.initPrior_peak();
      }
      else {
	cerr<<"Prior by QTLs\n"<<endl; 
	bn.initPrior(threthold_e);
      }
    }
    else if (priorfile.length()>0) {
      ifstream pfile(priorfile.c_str()); 
      if(!pfile) {
        cout <<"Can not open prior file "<<priorfile<<"\n";
	exit(-1);
      }
      else {
	/* read in prior infomation
	 * each line in the file is the following format
	 * parent_name   child_name   log(prior)    mutualinformation
	 */
	while(pfile.getline(line,maxlen,'\n') ) {
	    // comment line starting with #, so skip it
	    if(line[0] != '#') {
	        istrstream istr(line);
		string stmp,np, nc;
		float lp;
		istr>>np>>stmp>>nc>>lp;
		//cout<<np<<" "<<nc<<" "<<lp<<endl;
		Relationship r = Relationship(bn.getNodeIdByName(np),
					      bn.getNodeIdByName(nc));
		r.setPrior(lp);
		bn.setPrior(&r);
	    }  
	}
      }
    }
    else {//non-informative priors
      bn.initPrior();
      cout<<"non-informative priors\n";
    }

    //prohibited relatioships
    NEWP(bannedList, bn.getSize(), bool);
    for (int i=0; i<bn.getSize(); i++) {
        NEW(bannedList[i], bn.getSize(), bool);
	for (int j=0; j<bn.getSize(); j++) {
	  bannedList[i][j] = 0;
	}
    }
    //read banned list
    if(bannedListFile.length()>0) {
      //read in data
      ifstream dfile(bannedListFile.c_str()); 
      if(!dfile) {
        cout <<"Can not open banned list file "<<bannedListFile<<"\n";
	exit(-1);
      }
      else {
	int i=0;
	while(dfile.getline(line,maxlen,'\n') ) {
	  istrstream istr(line);
	  for(int j=0; j<bn.getSize(); j++) {
	    istr>>bannedList[i][j];
	  }
	  i++;
	}
	dfile.close();
      }
    }  
    /* A noive way to decide who is whose parent for time series
     * The rules are: 
     *  (1) if two genes have the same number of transition points, then 
     *      the gene with earlier transition point will be parent. If the 
     *      transition points coincide with each other, then check the next
     *      ones. If they have all the same transition points, they can be 
     *      parent of each other.
     *  (2) if two genes have different number of transition points, then
     *      the gene with first early transition point will be the parent. If
     *      the first transition points coincide, then check the next one. 
     *      If they have all the same transition points, they can be parent 
     *      of each other.
     */

    //calculate probability of node I being parent of node J
    #ifdef DEBUG
    for(int i=0; i<bn.getSize(); i++) {
      Node* ni = bn.getNode(i);
      cout<<"Node #"<<i<<": "<<ni->getName()<<endl;
      cout<<ni->getPrior(0)<<endl;
    }
    #endif

    //float threthold = log(threthold_e);
    //float threthold = log(1/float(bn.getSize()));
    ratio = log(ratio);
    qratio_cutoff = log(qratio_cutoff);
    for(int i=0; i<bn.getSize(); i++) {
      Node* ni=bn.getNode(i);
      for(int j=0; j<bn.getSize(); j++) {
	if(i!=j & bannedList[i][j] ==0) {
	  Node* nj=bn.getNode(j);
	  float pi=bn.getPrior(i,j);
	  if(listF && bannedList[i][j]==0) {
	      cout<<ni->getName()<<" -> "<<nj->getName()<<"   "<<pi+ratio<<" "
		  <<ni->mutualInfo(nj)<<endl;
	      cout<<"prior "<<ni->getName()<<" "<<nj->getName()<<" "
		  <<bn.getPrior(i,j)+ratio
		  <<" "<<bn.getNode(i)->getEQTL()->getComplexity(qcutoff)
		  <<" "<<bn.getNode(j)->getEQTL()->getComplexity(qcutoff)<<endl;
	  }

	  if(pi>qratio_cutoff && bannedList[i][j] ==0 && 
             (threthold_i <0 || ni->mutualInfo(nj)>threthold_i)) {
	    //re-scale the priors
	    Relationship r = Relationship(i,j);
	    
	    //if(pi>qratio_cutoff && bannedList[i][j] ==0) {
	    r.setPrior(pi+ratio);
	    bn.setPrior(&r);
	        #ifdef DEBUG
		cout<<ni->getName()<<" -> "<<nj->getName()<<"   "<<pi<<" "
		    <<ni->mutualInfo(nj)<<endl;
                #endif
		bn.addPotentialRelation(i, j);
	    //}
	  }  
	}
      }
    }

    //free banned list
    FREEP(bannedList, bn.getSize());
    
    //free line for text file reading
    free(line);

    //list priors
    if(listF) {
      exit(0);
    }

    //load trylist
    if(trylistFile.length()>0) {
      bn.loadTryList(trylistFile);
    }
    if(fitF) {
      cout<<"LIKELIHOOD "<<bn.BIC()<<endl;
      if(xmlfile.length()>0) bn.toXMLBIF(xmlfile);
      if(BIFfile.length()>0)  bn.toBIF(BIFfile);
      if(BNTfile.length()>0)  bn.toBNT(BNTfile);
      exit(0);
    }

    //update trylist
    if(updateNodeId>=0) {
        bn.updateTryList(updateNodeId);
	bn.outputTryList(updateTryListFile);
        exit(0);
    }

    //enumerate structures
    if(enumerationFile.length()>0) {
        bn.enumerateStructures(enumerationFile);
	bn.outputTryList(updateTryListFile);
        exit(0);
    }
    
    bn.learnStructure();
    //bn.printStructure(); 
    cout<<"LIKELIHOOD "<<bn.BIC()<<endl; 
    if(updateTryListFile.length()>0) bn.outputTryList(updateTryListFile);
    if(outfile.length()>0)  bn.toDotInput(outfile);
    if(xmlfile.length()>0)  bn.toXMLBIF(xmlfile);
    if(BIFfile.length()>0)  bn.toBIF(BIFfile);
}
