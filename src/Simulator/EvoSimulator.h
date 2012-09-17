/*
 *  EvoSimulator.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/25/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef EVO_SIMULATOR_
#define EVO_SIMULATOR_

#include <Base/Simulator.h>

namespace GPPG {
	
	class EvoSimulator : public GenotypeSimulator {
	public:
		EvoSimulator(IGenotypeHeap* heap);
		
		void addGenotype(IGenotype* g);
		
		void addGenotype(IGenotype* g, double freq);
		
		const std::set<IGenotype*>& activeGenotypes() const;
		
		/** Returns the number of active genotypes.
		 */
		int activeCount() const;
		
		/** Returns the current generation.
		 */
		int clock() const;
		
		/** Evolve the population of size \param N for \param G generations.
		 */
		void evolve2(long N, long G);
		void evolve(long N, long G);
		
	protected:
		IGenotype* activateGenotype(IGenotype* g, double freq);
		
		IGenotype* handleGenotype(IGenotype* g);
		
		void retireGenotype(IGenotype* g);
		
		void finishGeneration();
		
		void compactActive(long N);
		
		/** Generates a random parent from the individual array
		 */
		int randomParent();
		
	private:
		void checkIndividuals(long N);
		
		int _curr_gen;
		std::set<IGenotype*> _active;
		std::vector<IGenotype*> _ind1, _ind2;
		std::vector<IGenotype*> *_indIn, *_indOut;
		bool _indDirty;
	};
	
}
#endif