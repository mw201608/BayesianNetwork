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

/*
 *  Bayesian network SAX parser handler
 */

#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <iostream>
#include "BNHandler.h"

BNHandler::BNHandler(BN* bn) {
  this->bn = bn;
  NEW(currentPrior, 1, Relationship);
  NEW(currentReachable, 1, Relationship);
}

// ---------------------------------------------------------
//  BNHandlers: Overrides of the SAX ErrorHandler interface
// ---------------------------------------------------------
void BNHandler::error(const SAXParseException& e)
{
    cerr << "\nError at file " << e.getSystemId()
	 << ", line " << e.getLineNumber()
	 << ", char " << e.getColumnNumber()
         << "\n  Message: " << e.getMessage() << endl;
}

void BNHandler::fatalError(const SAXParseException& e)
{
    cerr << "\nFatal Error at file " << e.getSystemId()
	 << ", line " << e.getLineNumber()
	 << ", char " << e.getColumnNumber()
         << "\n  Message: " << e.getMessage() << endl;
}

void BNHandler::warning(const SAXParseException& e)
{
    cerr << "\nWarning at file " << e.getSystemId()
	 << ", line " << e.getLineNumber()
	 << ", char " << e.getColumnNumber()
         << "\n  Message: " << e.getMessage() << endl;
}

void BNHandler::endDocument()
{
}

void BNHandler::endElement(const XMLCh* const uri,
			   const XMLCh* const localname,
			   const XMLCh* const qname)
{
    //cout<<"parent.size()="<<parent.size()<<endl;
    parent.pop();
    if(parent.size()>0) {
        string top = parent.top();
	string name = XMLString::transcode(localname);
	//cout<<"name="<<name<<", top="<<top<<", value="<<value<<endl;
	if(name.compare("NAME")==0) {
	    if(top.compare("NETWORK")==0) {
	        bn->setName(value);
	    }
	    else if(top.compare("VARIABLE")==0) {
	        // update node name
	        currentNode->setName(value);
	    }	 
	}
	else if(name.compare("TYPE") ==0) {
	    if(top.compare("VARIABLE") ==0) {
	      currentNode->setType(value);
	    }
	}
	else if(name.compare("VALUE") ==0) {
	    if(top.compare("VARIABLE") ==0) {
	      currentNode->addNotion(value);
	    }
	}
	else if(name.compare("FOR") ==0) {
	  if(top.compare("PRIOR") ==0) {
	    currentPrior->setChild(bn->getNodeIdByName(value));
	  }
	}
	else if(name.compare("GIVEN") ==0) {
	  if(top.compare("PRIOR")==0) {
	    currentPrior->setParent(bn->getNodeIdByName(value));
	  }
	}
	else if(name.compare("TO") ==0) {
	  if(top.compare("REACHABLE") ==0) {
	    currentReachable->setChild(bn->getNodeIdByName(value));
	  }
	}
	else if(name.compare("FROM") ==0) {
	  if(top.compare("REACHABLE")==0) {
	    currentReachable->setParent(bn->getNodeIdByName(value));
	    bn->addFromNode(bn->getNodeIdByName(value));
	  }
	}
	else if(name.compare("PROPERTY") ==0) {
	  if(top.compare("PRIOR")==0) {
	    currentPrior->setPrior(atof(value));
	    bn->setPrior(currentPrior);
	    //set fixed relationship if prior is 1
	    if(currentPrior->getPrior()==1) {
	        //cout<<"fixed relationship "<<currentPrior->getParent()<<" "
                //      <<currentPrior->getChild()<<endl;
                bn->addFixedRelation(currentPrior->getParent(), 
				     currentPrior->getChild());
            }
	  }
	  else if(top.compare("REACHABLE")==0) {
	    currentReachable->setPrior(atof(value));
	    bn->setReachable(currentReachable);
	  }
	}
	else if(name.compare("NETWORK") ==0) {
	    //when we see network, all nodes' info should be input
	    //note info may be already updated if there are priors 
	    bn->updateNodesInfo();
	}
    }
}

void BNHandler::characters(const XMLCh* const chars, const unsigned int length) {
    if(length>0) {
        value = XMLString::transcode(chars);
        #ifdef DEBUG
	cout<<"value="<<value<<endl;
        #endif
    }
}

void BNHandler::startDocument()
{
}

void BNHandler::startElement(const   XMLCh* const    uri,
			     const   XMLCh* const    localname,
			     const   XMLCh* const    qname,
			     const   Attributes&     attributes)
{
    // don't care naming space (uri) much
    #ifdef DEBUG
    cout<< "localname="<<XMLString::transcode(localname)<<endl;
    cout<< "qname="<<XMLString::transcode(qname)<<endl;
    #endif
        
    string name = XMLString::transcode(localname);
    parent.push(name);
    if(string("NETWORK").compare(name)==0) {
        unsigned int len = attributes.getLength();
	for (unsigned int index = 0; index < len; index++) {
          #ifdef DEBUG
	  cout<<"Attribute # "<<index <<endl;
	  cout<<"localname="
	      << XMLString::transcode(attributes.getLocalName(index))
	      <<endl;
	  cout<<"qname="<< XMLString::transcode(attributes.getQName(index))
	      <<endl   ;
	  cout<<"value="<< XMLString::transcode(attributes.getValue(index))
	      <<endl;
          #endif
	  if(string("size").
	     compare(XMLString::transcode(attributes.getQName(index)))==0) {
	    bn->setSize(atoi(XMLString::transcode(attributes.getValue(index))));
	  }
	}
        bn->initSpace();
    }
    else if(string("VARIABLE").compare(name)==0) {
        //create a new node 
        Node node = Node("TMP");
	bn->addNode(node);
	currentNode = bn->getNodeByName("TMP");
    }
    else if(string("PRIOR").compare(name)==0) {
        //when we see PRIOR or REACHABLE, all nodes' info should be input
        bn->updateNodesInfo();
        //create a new relationship
        Relationship r = Relationship();
	*currentPrior = r;
    }
    else if(string("REACHABLE").compare(name)==0) {
        //when we see PRIOR or REACHABLE, all nodes' info should be input
        bn->updateNodesInfo();
        //create a new relationship
        Relationship r = Relationship();
	*currentReachable = r;
    }
}
