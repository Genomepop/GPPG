#include <iostream>
using namespace std;

#include "Operation/GreedyLoad.h"
#include <Model/Sequence/Operation.h>
#include <Model/Sequence/IO.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <Base/Simulator.h>
#include <Operation/OperationHeap.h>
#include <Operation/BaseCompressionPolicy.h>
#include <Operation/UbiOperationGraph.h>
#include <Simulator/EvoSimulator.h>

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
	
	//PopulationSimulator psim( new OperationGraph(new BaseCompressionPolicy(STORE_ACTIVE)) );
	EvoSimulator psim( new OperationGraph(new BaseCompressionPolicy(STORE_ROOT)) );
	//EvoSimulator psim( new OperationGraph(new GreedyLoad(20, 10) ));
	
	ublas::vector<double> distr = ublas::vector<double>(4);
	for (int i=0; i<distr.size(); i++) {
		distr(i) = 1.0/distr.size();
	}
	double scaling = 1e2;
	long N = 1e4;
	long L = 1e6;
	long G = 2e7;
	double u = 1e-9;
	double ud = 1e-9;
	double ui = 1e-9;
	double ur = 1e-9;
	cout << "N = " << N/scaling << endl;
	cout << "G = " << G/scaling << endl;
	cout << "u = " << u*scaling << endl;
	
	SequenceRootFactory factory(L, distr);
	SequenceRoot* sr = factory.random();
	
	ublas::matrix<double> T(4,4);
    for (unsigned i = 0; i < T.size1(); ++ i)
        for (unsigned j = 0; j < T.size2(); ++ j)
            T (i, j) = 0.25;
	
	SequencePointMutator *spm = new SequencePointMutator( u*scaling, T);
	SequenceDeletionMutator *sdm = new SequenceDeletionMutator( ud*scaling , 10, 20);
	SequenceInsertionMutator *sim = new SequenceInsertionMutator( ui*scaling , 10, 20, distr);

	psim.addGenotype( sr, 1.0 );
	psim.addMutator( spm );
	//psim.addMutator(sdm);
	//psim.addMutator(sim);

	
	psim.addRecombinator(new SequenceRecombinator(ur*scaling) );
	int steps = 1000;
	for (int i=0; i<steps; i++) {
		psim.evolve( N/scaling, (G/scaling)/steps );	
		//usleep(50000);
		cout << "Done with " << i << endl;
	}
	//psim.evolve( 500, 10000 );
	
	cout << "Generation: " << psim.clock() << endl;
	set<IGenotype*>::iterator git;
	const set<IGenotype*>& active = psim.activeGenotypes();
	int i=0;
	for (git=active.begin(); git!=active.end(); git++) {
		IGenotype* g = *git;
		cout << "Genotype " << i << "/" << g->order() << ": " << g->frequency() << endl;
		i++;
	}
	
}

int main (int argc, char * const argv[])
{
	//testSequence();
	testSimulator();
	return 0;
}

