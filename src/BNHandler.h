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
 *  XML handler (SAX2 based) for Bayesian network.
 *
 */

#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/framework/XMLFormatter.hpp>
#include <stack>
#include "BN.h"

XERCES_CPP_NAMESPACE_USE

class BNHandler : public DefaultHandler {
public:
    BNHandler(BN* bn); 
    ~BNHandler(){}


    // -------------------------------------------------------
    //  Implementations of the SAX DocumentHandler interface
    // -------------------------------------------------------
    void endDocument();

    void endElement( const XMLCh* const uri, 
		     const XMLCh* const localname, 
		     const XMLCh* const qname);

    void characters(const XMLCh* const chars, const unsigned int length);
    void startDocument();
    void startElement(	const   XMLCh* const    uri,
			const   XMLCh* const    localname,
			const   XMLCh* const    qname,
			const   Attributes&	attributes);

    // -----------------------------------------------------------------------
    //  Implementations of the SAX ErrorHandler interface
    // -----------------------------------------------------------------------
    void warning(const SAXParseException& exception);
    void error(const SAXParseException& exception);
    void fatalError(const SAXParseException& exception);


private:
    BN* bn;
    Node* currentNode;
    Relationship* currentPrior;
    Relationship* currentReachable;
    char* value;
    stack <string> parent;
};
