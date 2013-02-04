#include <Model/Pathway/Fitness.h>
#include <Model/Pathway/Operation.h>

using namespace GPPG::Model::TransReg;

double ConnectedFitness::calculate(OpPathway* g) {
	BindingSiteChange* bsc = dynamic_cast<BindingSiteChange*>(g);
	if( bsc ) {
		// Check if any delta sites are 0
		if( bsc->areAllGenesRegulated() ) return 1.0;
		return 0.0;
	} 
	
	// If we can't use our shortcut, then we have to check all the genes!
	for(int i=0; i<g->numGenes(); i++) {
		if( g->numSitesForGene(i) == 0 )
			return 0.0;
	}
	return 1.0;
}