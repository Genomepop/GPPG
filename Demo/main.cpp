#include <iostream>
using namespace std;

//#include "GreedyLoad.h"
#include <Model/Sequence/Operation.h>
#include <Model/Sequence/IO.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <Base/Simulator.h>

using namespace GPPG;
using namespace GPPG::Model;
using namespace boost::numeric;

void testSequence() {
	//GreedyLoad gl(3);
	
	ublas::vector<double> distr = ublas::vector<double>(4);
	for (int i=0; i<distr.size(); i++) {
		distr(i) = 1.0/distr.size();
	}
	
	cout << distr << endl;
	SequenceRootFactory factory(10, distr);
	
	SequenceRoot* sr = factory.random();
	SequenceData& sd = *sr->data();
	cout << sd << endl; 
	
	ublas::matrix<double> T(4,4);
    for (unsigned i = 0; i < T.size1(); ++ i)
        for (unsigned j = 0; j < T.size2(); ++ j)
            T (i, j) = 0.25;
	
	SequencePointMutator spm(0.1, T);
	
	cout << spm << endl;
	
	cout << "About to mutate.\n";
	OpSequence* m1 = spm.mutate( *sr );
	cout << "Finished mutate.\n";
	cout << *m1 << endl;
	cout << m1->length() << endl;
	for (int i=0; i<m1->length(); i++) {
		cout << m1->get(i);
	}
	cout << endl;
	cout << "Hello World2!\n";
}

void testSimulator() {
	PopulationSimulator psim(0);
	
	ublas::vector<double> distr = ublas::vector<double>(4);
	for (int i=0; i<distr.size(); i++) {
		distr(i) = 1.0/distr.size();
	}
	SequenceRootFactory factory(10, distr);
	SequenceRoot* sr = factory.random();
	
	ublas::matrix<double> T(4,4);
    for (unsigned i = 0; i < T.size1(); ++ i)
        for (unsigned j = 0; j < T.size2(); ++ j)
            T (i, j) = 0.25;
	
	SequencePointMutator *spm = new SequencePointMutator(0.01, T);
	
	psim.addGenotype( sr, 1.0 );
	psim.addMutator( (IMutator*)spm );
	
	psim.evolve( 20, 10 );
	
	cout << "Generation: " << psim.clock() << endl;
	for (int i=0; i<psim.activeCount(); i++) {
		IGenotype& g = psim.genotype(i);
		cout << "Genotype " << i << "/" << g.order() << ": " << g.frequency() << endl;
	}
	
}

int main (int argc, char * const argv[])
{
	//testSequence();
	testSimulator();
	return 0;
}

