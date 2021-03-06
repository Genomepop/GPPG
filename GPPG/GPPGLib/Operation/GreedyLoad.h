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
#include <map>

namespace GPPG {

	
	class GreedyLoad : public CompressionPolicy {
	public:
		/** Create a GreedyLoad policy that uses at most \param maxExplicit uncompressed genotypes and is applied every \param numGens generations.
		 */
		GreedyLoad(int maxExplicit, int numGens);
		
		void decompressionReleased( IOperation* op );
		
		void operationAdded(IOperation* op);
		
		void generationFinished( OperationGraph* heap, const std::set<IOperation*>& active );
		
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
		void remove(IOperation* op);
		bool isCompressed(IOperation* op);
		void move(IOperation* a, IOperation* b);
		void split( IOperation* op, int& s1, int& s2, IOperation*& g1, IOperation*& g2 );
		
		// Methods for annotating load
		void setLoad(IOperation* op, int c);
		void incrLoad(IOperation* op, int c);
		void decrLoad(IOperation* op, int c);
		int load(IOperation* op);
		void clearLoadMap();
		void annotate(const std::set<IOperation*>&);
		void reset(IOperation* op);
		IOperation* findMaxAdvance(IOperation* op, bool doReset);
		IOperation* uncoveredChild(IOperation* op);
		void reverseAnnotate(IOperation* op, int c);
		void resetAnnotation(IOperation* op);
		
		
		IOperation* getMaxItem( const std::set<IOperation*>& items, bool compare );
		
		std::set<IOperation*> _U, _V;
		IOperation* _root;
		int _maxExplicit, _elapsedGens, _numExplicit, _waitGens, _runs;
	};
}
#endif
