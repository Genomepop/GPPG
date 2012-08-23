/*
 *  Simulator.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/17/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef	SIMULATOR_
#define SIMULATOR_

#include <vector>

namespace GPPG {
	class IMutator;
	class IGenotype;
	class IGenotypeHeap;
	
	class GenotypeSimulator {
	public:
		GenotypeSimulator(IGenotypeHeap* heap);
			
		void addMutator(IMutator* mutator);
		
		virtual void addGenotype(IGenotype* g);
		
		IGenotypeHeap* heap();
		
	protected:
		virtual void configureGenotype(IGenotype *g);
		virtual void removeGenotype(IGenotype *g);
		virtual IGenotype* handleGenotype(IGenotype *g);
		
		std::vector<IMutator*> _mutators;
		std::vector<IGenotype*> _genotypes;
		
	private:
		IGenotypeHeap* _heap;
	};
	
	class PopulationSimulator : public GenotypeSimulator {
	public:
		PopulationSimulator(IGenotypeHeap* heap);
		
		void addGenotype(IGenotype* g);
		
		void addGenotype(IGenotype* g, double freq);
		
		/** Retrieves the \param i'th genotype.
		 */
		IGenotype& genotype(int i);
		
		/** Retrieves the frequency of the \param i'th genotype
		 */
		double frequency(int i);
		
		/** Returns the number of active genotypes.
		 */
		int activeCount() const;
		
		/** Returns the current generation.
		 */
		int clock() const;
		
		/** Evolve the population of size \param N for \param G generations.
		 */
		void evolve(long N, long G);
		
		
	protected:
		IGenotype* activateGenotype(IGenotype* g, double freq);
		
		IGenotype* handleGenotype(IGenotype* g);
		
		void retireGenotype(IGenotype* g);
		
		void finishGeneration();
		
		void compactActive(long N);
		
	private:
		int _curr_gen;
		std::vector<IGenotype*> _active;
		std::vector<double> _freqs, _freqs_m, _freqs_r;
	};
}

#endif