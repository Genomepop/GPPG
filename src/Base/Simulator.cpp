/*
 *  Simulator.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 8/17/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "Simulator.h"
#include "Base/Genotype.h"
#include "Base/Mutator.h"

using namespace GPPG;

GenotypeSimulator::GenotypeSimulator() {}

void GenotypeSimulator::addMutator(IMutator* mutator) {
	_mutators.push_back( mutator );
}

void GenotypeSimulator::addGenotype(IGenotype* g) {
	configureGenotype( g );
}

void GenotypeSimulator::configureGenotype(IGenotype *g) {

	g->configure();
	
	//TODO: Evaluate phenotype function
	
	//TODO: Evaluate fitness function
	
	g->setTotal(0);
	
	// Add genotype to lookup
	_genotypes.push_back(g);
	
}

IGenotype* GenotypeSimulator::handleGenotype(IGenotype *g) {
	addGenotype(g);
	return g;
}

void GenotypeSimulator::removeGenotype(IGenotype *g) {
	
}

PopulationSimulator::PopulationSimulator(): _curr_gen(0) {}

void PopulationSimulator::addGenotype(IGenotype* g, double freq) {
	GenotypeSimulator::addGenotype(g);
	g->setIndex(-1);
	if (freq > 0.0) {
		activateGenotype(g, freq);
	}
}


IGenotype& PopulationSimulator::genotype(int i) {
	return *_active[i];
}


double PopulationSimulator::frequency(int i) {
	return _freqs[i];
}


int PopulationSimulator::activeCount() const {
	return _active.size();
}


int PopulationSimulator::clock() const {
	return _curr_gen;
}


void PopulationSimulator::evolve(long N, long G) {
	
	int num_active;
	double one_individual = 1.0/N;
	
	while (_curr_gen < G) {
		// Clear the frequency arrays, reset status of genotypes
		for (int i=0; i<_freqs.size(); i++) {
			_freqs_m[i] = _freqs[i];
		}
		num_active = _freqs.size();
		
		for (int i=0; i<_mutators.size(); i++) {
			IMutator* mutator = _mutators[i];
			for (int gi=0; gi<num_active; gi++) {
				IGenotype* g1 = _active[i];
				int num_mutants = mutator->numMutants(*g1, N, _freqs[gi]);
				for (int muti; muti<num_mutants; muti++) {
					IGenotype* g2 = handleGenotype( mutator->mutate( *g1 ) );
					if (g2 != NULL)
						_freqs_m[ g2->index() ] += one_individual;
				}
				_freqs_m [ g1->index() ] -= num_mutants*one_individual;
			}
		}
		_curr_gen ++;
	}
}

IGenotype* PopulationSimulator::handleGenotype(IGenotype* g) {
	if (g == NULL) return NULL;
	GenotypeSimulator::handleGenotype(g);
	
	return activateGenotype(g, 0.0);
}

void PopulationSimulator::compactActive(long N) {
	int nactive = 0; // Keeps track of the number of active genotypes
	
	for (int i=0; i<_active.size(); i++) {
		if (_freqs[i]*N >= 1) {
			if (i > nactive) {
				_active[nactive] = _active[i];
				_freqs[nactive] = _freqs[i];
				_active[nactive]->setIndex( nactive );
			}
			nactive++;
		} else {
			// Retire the genotype 
			IGenotype* g = _active[i];
			retireGenotype( g );
			removeGenotype( g );
		}
	}
	
	_active.resize( nactive );
	_freqs.resize( nactive );
	_freqs_m.resize( nactive );
	_freqs_r.resize( nactive );
}

IGenotype* PopulationSimulator::activateGenotype(IGenotype* g, double freq) {
	if (g->index() >= 0) {
		return g;
	}
	
	// Put g at the end of the list
	int idx = _active.size();
	g->setIndex(idx);
	g->setOrder( clock() );
	
	_active.push_back( g );
	_freqs.push_back( freq );
	_freqs_m.push_back( freq );
	_freqs_r.push_back( freq );
	
	return g;
}

void PopulationSimulator::retireGenotype(IGenotype* g) {
	g->setFrequency(0);
	g->setIndex(-1);
}

void PopulationSimulator::finishGeneration() {

}