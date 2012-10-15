/*
 *  Simulator.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 8/17/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "GPPG.h"
#include "Simulator.h"
#include "Base/Genotype.h"
#include "Base/Mutator.h"
#include "Base/GenotypeHeap.h"
#include "Base/Recombinator.h"

#include "Util/Random.h"
#include <iostream>

using namespace GPPG;

void printArray(std::vector<double>& arr) {
	std::cout << "[";
	for (int i=0; i<arr.size(); i++) {
		std::cout << arr[i] << " ";
	}
	std::cout << "]";
}

inline void normalizeArray( std::vector<double>& arr ) {
	double csum = 0.0;
	for (int i=0; i<arr.size(); i++) csum += arr[i];
	for (int i=0; i<arr.size(); i++) arr[i] /= csum;
}




void samplePopulation( std::vector<double>& F, long N, std::vector<double>& res ) {
	int ig = 0;
	long n = N;
	int draw = -1;
	double sump, pp, tot;
	sump = 0.0;
	tot = 0.0;
	int k = F.size();
	
	while (ig < k) {
		if (n > 0 && F[ig] > 0) {
			pp = F[ig]/(1.0-sump);
			if (pp >= 1.0)
				draw = n;
			else 
				draw = binomial(n, pp);

			res[ig] = (1.0*draw)/N;
			tot == res[ig];
			n = n-draw;
			sump += F[ig];
		} else {
			res[ig] = 0;
		}
		ig += 1;
	}

}

GenotypeSimulator::GenotypeSimulator(IGenotypeHeap* h): _heap(h), _gcount(0) {}

IGenotypeHeap* GenotypeSimulator::heap() { return _heap; }

void GenotypeSimulator::addMutator(IMutator* mutator) {
	_mutators.insert( mutator );
}

void GenotypeSimulator::addRecombinator(IRecombinator* r) {
	_recombinators.insert( r );
}

void GenotypeSimulator::addGenotype(IGenotype* g) {
	configureGenotype( g );
	_heap->addGenotype( g );	
}

void GenotypeSimulator::configureGenotype(IGenotype *g) {
	// Give the genotype a unique key
	g->setKey(_gcount++);
	g->configure();
	
	//TODO: Evaluate phenotype function
	
	//TODO: Evaluate fitness function
	
	g->setTotal(0);
	
	// Add genotype to lookup
	//_genotypes.insert(g);
	
}

IGenotype* GenotypeSimulator::handleGenotype(IGenotype *g) {
	addGenotype(g);
	return g;
}

void GenotypeSimulator::removeGenotype(IGenotype *g) {
	_heap->removeGenotype( g );
}

PopulationSimulator::PopulationSimulator(IGenotypeHeap* h): GenotypeSimulator(h), _curr_gen(0) {
	initRandom();
}

void PopulationSimulator::addGenotype(IGenotype* g) {
	addGenotype(g, 0.0);
}

void PopulationSimulator::addGenotype(IGenotype* g, double freq) {
	GenotypeSimulator::addGenotype(g);
#ifdef DEBUG_0
	std::cout << "PopulationSimulator::addGenotype " << g << ": " << freq;
#endif
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
	
	normalizeArray( _freqs );
	long Gtot = _curr_gen + G;
	
	while (_curr_gen < Gtot) {
		// Clear the frequency arrays, reset status of genotypes
		for (int i=0; i<_freqs.size(); i++) {
			_freqs_m[i] = _freqs[i];
		}
		num_active = _freqs.size();
		
		for (std::set<IMutator*>::iterator it = _mutators.begin(); it!=_mutators.end(); it++) {
		//for (int i=0; i<_mutators.size(); i++) {
			IMutator* mutator = *it; //_mutators[i];
			for (int gi=0; gi<num_active; gi++) {
				IGenotype* g1 = _active[gi];
				int num_mutants = mutator->numMutants(*g1, N, _freqs[gi]);
#ifdef DEBUG_0
				std::cout << num_mutants << " Mutants\n";
#endif
				for (int muti=0; muti<num_mutants; muti++) {
#ifdef DEBUG_0
					std::cout << "Mutating " << muti << std::endl;
#endif
					addGenotype( mutator->mutate( *g1 ), one_individual );
					//if (g2 != NULL)
					//	_freqs_m[ g2->index() ] += one_individual;
				}
				_freqs_m [ g1->index() ] -= num_mutants*one_individual;
			}
		}
		
		// Do recombination
		// TODO: recombination
		
		// Apply fitness
		for (int i=0; i<_active.size(); i++) {
			_freqs_m[i] *= _active[i]->fitness();
		}
		
		// Normalize array
		normalizeArray( _freqs_m );
		
		// Genetic Drift
		//std::cout << "DRIFT\n";
		//printArray( _freqs_m ); std::cout<<std::endl;
		samplePopulation( _freqs_m, N, _freqs);
		//printArray( _freqs); std::cout<<std::endl;
		

		compactActive(N);
		
		//std::cout << "COMPACT ACTIVE\n";
		//printArray( _freqs); std::cout<<std::endl;

		normalizeArray( _freqs );
		
		//std::cout << "NORMED\n";
		//printArray( _freqs); std::cout<<std::endl;
		
		finishGeneration();
		
		
		_curr_gen ++;
		
		
		
		// TODO: Record information
		//for (int i=0; i<_active.size(); i++) {
		//	_active[i]->setFrequency( _freqs[i] );
		//}
	}
	
	for (int i=0; i<_active.size(); i++) {
		_active[i]->setFrequency( _freqs[i] );
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
	_heap->generationFinished( _active );
}