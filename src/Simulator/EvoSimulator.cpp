/*
 *  EvoSimulator.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 8/25/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "EvoSimulator.h"
#include "GPPG.h"
#include "Base/Genotype.h"
#include "Base/Mutator.h"
#include "Base/GenotypeHeap.h"
#include "Base/Recombinator.h"

#include <iostream>
#include <sstream>

//#include "Util/Random.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/binomial_distribution.hpp>

#ifdef UBIGRAPH
extern "C" {
#include <Util/Ubigraph/ubiclient.h>
}
#endif

using namespace GPPG;
using std::set;
using std::cout;
using std::endl;

extern boost::mt19937 gen;

typedef set<IGenotype*>::iterator GIter;

template <class T> std::string TToStr( const T &t )
{
    std::ostringstream oss;
    oss << t;
    return std::string (oss.str());
}

inline void normalizeArray( std::set<IGenotype*>& genos ) {
	double csum = 0.0;
	for (GIter git=genos.begin(); git!=genos.end(); git++) csum += (*git)->frequency();
	for (GIter git=genos.begin(); git!=genos.end(); git++) (*git)->setFrequency((*git)->frequency()/csum);	
}

inline int binomial(int n, double r) {
	boost::random::binomial_distribution<> dist( n, r );
	return dist(gen);
}

void printGenos( set<IGenotype*>& genos ) {
	cout << "[";
	for (GIter git=genos.begin(); git!=genos.end(); git++) cout << (*git)->frequency() << ",";
	cout << "]" << endl;
}

void samplePopulation( set<IGenotype*>& genos, long N) {
	int ig = 0;
	long n = N;
	int draw = -1;
	double sump, pp, tot, res, f;
	sump = 0.0;
	tot = 0.0;
	int k = genos.size();
	for (GIter git=genos.begin(); git!=genos.end(); git++) {
		IGenotype* g = *git;
		f = g->frequency();
		res = 0;
		if (ig == k-1) {
			if (f > 0) res = 1.0-tot;
		} else if (n > 0 && f > 0) {
			pp = f/(1.0-sump);
			if (pp >= 1.0) {
				draw = n;
			} else {
				draw = binomial(n, pp);
			}
			res = (1.0*draw)/N;
			tot += res;
			n = n-draw;
			sump += f;
		} 
		
		g->setFrequency( res );
		ig++;
	}
	
}

EvoSimulator::EvoSimulator(IGenotypeHeap* h): GenotypeSimulator(h), _curr_gen(0) {
	gen.seed((unsigned int)time(0));
}

void EvoSimulator::addGenotype(IGenotype* g) {
	addGenotype(g, 0.0);
}

void EvoSimulator::addGenotype(IGenotype* g, double freq) {
	GenotypeSimulator::addGenotype(g);
	
	g->setIndex(-1);
	g->setState(-1);
	
	if (freq > 0.0) {
		activateGenotype(g, freq);
	}
}


int EvoSimulator::activeCount() const {
	return _active.size();
}


int EvoSimulator::clock() const {
	return _curr_gen;
}


void EvoSimulator::evolve(long N, long G) {
	
	double one_individual = 1.0/N;
	
	normalizeArray( _active );
	long Gtot = _curr_gen + G;
	
#ifdef DEBUG_0
	cout <<	"Starting Evo with " << _active.size() << " genotypes.\n";
#endif
	
	while (_curr_gen < Gtot) {
		
		for (std::set<IMutator*>::iterator it = _mutators.begin(); it!=_mutators.end(); it++) {
			IMutator* mutator = *it;
#ifdef DEBUG_0
			cout <<	"Using Mutator: " << mutator << endl;
#endif
			for (GIter git = _active.begin(); git != _active.end(); git++) {
				IGenotype* g1 = *git;
#ifdef DEBUG_0
				cout <<	"Geno: (" << g1->key() << ", " << g1->order() << ")" << endl;
#endif
				if (g1->order() == _curr_gen) continue;
				
#ifdef DEBUG_0
				cout <<	"Mutating: (" << g1->key() << ", " << g1->frequency() << ")" << endl;
#endif
				
				double freq = g1->frequency();
				int num_mutants = mutator->numMutants(*g1, N, freq);
#ifdef DEBUG_0
				cout <<	"Num Mutants: " << num_mutants << endl;
#endif
				for (int muti=0; muti<num_mutants; muti++) {
					
					addGenotype( mutator->mutate( *g1 ), one_individual );
				}
				g1->setFrequency( freq - num_mutants*one_individual );
			}
		}
		
		// Do recombination
		for (std::set<IRecombinator*>::iterator it = _recombinators.begin(); it!=_recombinators.end(); it++) {
			IRecombinator* recombinator = *it;
#ifdef DEBUG_0
			cout <<	"Using Recombinator: " << recombinator << endl;
#endif
			for (GIter git = _active.begin(); git != _active.end(); git++) {
				IGenotype* g1 = *git;
				if (g1->frequency() == 0) continue;
#ifdef DEBUG_0
				cout <<	"Geno: (" << g1->key() << ", " << g1->order() << ")" << endl;
#endif
				//if (g1->order() == _curr_gen) continue;
				for (GIter git2 = git; git2 != _active.end(); git2++) {
					if (*git2==*git) continue;
					IGenotype* g2 = *git2;
					if (g2->frequency() == 0) continue;
				
#ifdef DEBUG_0
					cout <<	"Recombining: (" << g1->key() << ", " << g2->key() << ") " << g1->frequency() << ": " << g2->frequency()<< endl;
#endif
				

					int num_mutants = recombinator->numMutants(*g1, *g2, N);
#ifdef DEBUG_0
					cout <<	"Num Genos: " << num_mutants << endl;
#endif
					for (int muti=0; muti<num_mutants; muti++) {
					
						addGenotype( recombinator->recombine( *g1, *g2 ), one_individual );
					}
					
					g1->setFrequency( g1->frequency() - num_mutants*one_individual/2 );
					g2->setFrequency( g2->frequency() - num_mutants*one_individual/2 );
				}
				
			}
		}
		
		// Apply fitness
		for (GIter git = _active.begin(); git != _active.end(); git++) {
			IGenotype* g1 = *git;
			g1->setFrequency( g1->fitness()*g1->frequency() );
		}
		
		// Normalize array
		normalizeArray( _active );
		
		// Genetic Drift
		samplePopulation( _active, N );
		
		compactActive(N);
		// This is probably not necessary!
		//normalizeArray( _active );
		
		finishGeneration();
		
		
		_curr_gen ++;
	}
	
}

IGenotype* EvoSimulator::handleGenotype(IGenotype* g) {
	if (g == NULL) return NULL;
	GenotypeSimulator::handleGenotype(g);
	
	return activateGenotype(g, 0.0);
}

void EvoSimulator::compactActive(long N) {
	//int nactive = 0; // Keeps track of the number of active genotypes
	GIter git = _active.begin();
	while (git != _active.end()) {
		IGenotype* g = *git;
		git++;
		if (g->frequency() <= 0) {
			retireGenotype( g );
			removeGenotype( g );
		}
	}
	
}

const set<IGenotype*>& EvoSimulator::activeGenotypes() const { return _active; }

IGenotype* EvoSimulator::activateGenotype(IGenotype* g, double freq) {
	//if (g->order() <= _curr_gen) {
	//	return g;
	//}
	
	// Put g at the end of the list
	g->setIndex(1);
	g->setState(1);
	g->setOrder( clock() );
	g->setFrequency( freq );
	_active.insert( g );
	
	
	return g;
}

void EvoSimulator::retireGenotype(IGenotype* g) {
	g->setFrequency(0);
	g->setIndex(-1);
	g->setState(-1);
	_active.erase( g );
}

void EvoSimulator::finishGeneration() {
#ifdef UBIGRAPH_SIM
	for (GIter git = _active.begin(); git != _active.end(); git++) {
		IGenotype* g1 = *git;
		ubigraph_set_vertex_attribute( g1->key(), "size", TToStr<double>(10*g1->frequency() ));
	}
#endif
	_heap->generationFinished( _active );
}
