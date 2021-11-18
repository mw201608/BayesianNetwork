/**
 *  test CMI (conditional mutual information) on a set of data 
 */

#include <iostream.h>
#include <fstream.h> 
#include <strstream.h>
#include <math.h>
#include <vector.h>
#ifdef LINUX
    #include <getopt.h>
#endif
#include "common.h"
#include "BN.h"

/**
 *$Revision: $
 *$Date: $
 *@Author:  Jun Zhu
 *@Creation Date: 5/6/01
 */

void usagemessage (char *name) {
    cout << "Usage: "<< name << " -s seed -b BIF file -d datafile [-e eQTL file] -t threthold_i -T threthold_e -D data dimension -1 n1 -2 n2 [-3 n3]\n";
    exit(0);
}

int main(int argc, char **argv) {
    string biffile;
    string datafile;
    string eqtlfile;
    float  threthold_i=0;
    float  threthold_e=0;
    int    dim=0;
    string n1, n2, n3;

    int option_char;
    while ((option_char = getopt (argc, argv, "b:e:d:t:D:1:2:3:")) != -1) {
        switch (option_char) {
	case 'b': biffile = optarg; break;
	case 'd': datafile=optarg; break;
	case 'e': eqtlfile = optarg; break;
	case 't': threthold_i=atof(optarg); break;
	case 'T': threthold_e=atof(optarg); break;
	case 'D': dim=atoi(optarg); break;
	case '1': n1 = optarg; break;
	case '2': n2 = optarg; break;
	case '3': n3 = optarg; break;
	default: usagemessage(argv[0]);
	}
    }

    /* data dimension is not automatically checked, so we need to
     * know this ahead of time.
     */
    if(dim==0) {
      cout<<"Data dimension is required.\n";
      usagemessage(argv[0]);
    }
    
    //create an empty Bayesian network;
    BN bn = BN();

    //read in initial Bayesian network if exists
    if(biffile.length()>0) {
        if(bn.init(biffile)!=0) exit(0);
    }

    int maxlen = 81920;
    char line[maxlen];

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
		float loc, lod;
		istr>>stmp>>loc>>lod;
		Node *nd= bn.getNodeByName(stmp);
		nd->getEQTL()->add(loc, exp(lod));
	    }
	}
	efile.close();   
      }
    }

    //read in data
    ifstream dfile(datafile.c_str()); 
    if(!dfile) {
        cout <<"Can not open data file "<<datafile<<"\n";
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
		    int datav;
		    istr>>datav;
		    nd->addData(datav);
		}
		nd->calculatePrior();
		nd->calculateLikelihood();  
	    }
	}
	dfile.close();   
    }

    Node* nd1 = bn.getNodeByName(n1);
    Node* nd2 = bn.getNodeByName(n2);
    if(n3.length()>0) {
      //calculate conditional mutual information
      cout<<"condition on "<<n1<<" "<<n2<<"--"<<n3<<" "
	  <<nd1->cmi(nd2, bn.getNodeByName(n3))<<endl;
      cout<<"mutual info "<<n2<<"--"<<n3<<" "
	  <<nd2->mutualInfo(bn.getNodeByName(n3))<<endl;
    }
    else {
      bn.initPrior(threthold_i, threthold_e);
      //caluclate mutual information
      int id1 = bn.getNodeIdByName(n1);
      int id2 = bn.getNodeIdByName(n2);
      cout<<n1<<" ";
      nd1->getEQTL()->printMode();
      cout<<n2<<" ";
      nd2->getEQTL()->printMode();
      cout<<n1<<" "<<n2<<" "<<nd1->mutualInfo(nd2)<<" "
	  <<bn.matcheQTL(nd1->getEQTL(), nd2->getEQTL())<<" "
          <<bn.getPrior(id1, id2)<<endl;
      cout<<n2<<" "<<n1<<" "<<nd2->mutualInfo(nd1)<<" "
	  <<bn.matcheQTL(nd2->getEQTL(), nd1->getEQTL())<<" "
          <<bn.getPrior(id2, id1)<<endl;
    }
}
