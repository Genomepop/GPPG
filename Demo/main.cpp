#include <iostream>
using namespace std;

//#include "GreedyLoad.h"
#include "SequenceOperation.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace GPPG;
using namespace GPPG::Model;
using namespace boost::numeric;

int main (int argc, char * const argv[])
{
	//GreedyLoad gl(3);
	
	ublas::vector<double> distr = ublas::vector<double>(4);
	for (int i=0; i<distr.size(); i++) {
		distr(i) = 1.0/distr.size();
	}
	
	cout << distr << endl;
	SequenceFactory factory(20, distr);
	
	SequenceRoot* sr = factory.random();
	SequenceData& sd = *sr->data();
	cout << sd << endl; 

	ublas::matrix<double> T(4,4);
    for (unsigned i = 0; i < T.size1(); ++ i)
        for (unsigned j = 0; j < T.size2(); ++ j)
            T (i, j) = 0.25;
	
	SequencePointMutator spm(0.00001, T);
	
	cout << spm << endl;
	//SequencePointChange( dop, 
	//SequenceData* sr = dop.evaluate();
	
	cout << "Hello World2!\n";
	return 0;
}