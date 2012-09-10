/*
 *  GreedyLoad.h
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef OPERATION_GREEDY_LOAD_
#define OPERATION_GREEDY_LOAD_

#include "Operation/CompressionPolicy.h"
#include <set>

namespace GPPG {
	class GreedyLoad : public CompressionPolicy {
	public:
		/** Create a GreedyLoad policy that uses at most \param maxExplicit uncompressed genotypes and is applied every \param numGens generations.
		 */
		GreedyLoad(int maxExplicit, int numGens);
		
		
		void decompressionReleased( IOperation* op );
		
		
		void operationAdded(IOperation* op);
		
		
		void generationFinished( const std::set<IOperation*>& active );
		
		/** Force an update by the policy.
		 * This resets the count of elapsed generations.
		 */
		void apply( const std::set<IOperation*>& active );
		
		/** Retrieve the maximum number of uncompressed genotypes.
		 * This is the 'k' parameter in the paper.
		 */
		int maxUncompressed() const;
		
		/** Retrieve the number of generations to elapse.
		 * This is the 't' parameter in the paper.
		 */
		int numGenerations() const;
		
	private:
		GreedyLoad(GreedyLoad const&);
		GreedyLoad const& operator=(GreedyLoad const&);
		
		
		
		void advance(IOperation* op);
		
		// Methods for managing the compression sets
		void add(IOperation* op);
		void remove(IOperation* op, bool cache, bool doRecurse);
		void move(IOperation* a, IOperation* b);
		void update();
		void clearCache();
		void split( IOperation* op, int& s1, int& s2, IOperation*& g1, IOperation*& g2 );
		IOperation* getMaxItem( const std::set<IOperation*>& items, bool compare );
		
		std::set<IOperation*> _U;
		IOperation* _root;
		int _maxExplicit, _elapsedGens, _numExplicit, _waitGens;
	};
}
#endif
