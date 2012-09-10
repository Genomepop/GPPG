/*
 *  OperationHeap.h
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */


#ifndef CORE_OPERATION_HEAP_
#define CORE_OPERATION_HEAP_

//#include "Operation/CompressionPolicy.h"
//#include "Operation/Operation.h"
#include "Base/GenotypeHeap.h"

namespace GPPG {
	class ICompressionPolicy;
	class IOperation;
	
	
	class OperationGraph : public IGenotypeHeap {
	public:
		OperationGraph(ICompressionPolicy* policy);
		
		~OperationGraph();
		
		/**
		 Retrieves the policy
		 */
		ICompressionPolicy& compressionPolicy();
		
		void setHeapRemovalPolicy(HeapRemovalPolicy policy);
		
		HeapRemovalPolicy heapRemovalPolicy() const;
		
		/** Casts the genotype to an IOperation
		 */
		void addGenotype(IGenotype* g);
		
		/** Casts the genotype to an IOperation.
		 */
		void removeGenotype(IGenotype* g);
		
		void generationFinished(const std::vector<IGenotype*>&);
		
		void generationFinished(const std::set<IGenotype*>&);
		
		/** Sets this IOperation to be owned by the graph
		 */
		virtual void addOperation(IOperation* op);
		
		/** Removes this IOperation from control of the graph.
		 * This may result in the deletion of this operation and any defuct operations resulting from its deletion.
		 */
		virtual void removeOperation(IOperation* op);
		
		
		//void operationAttached(IOperation& parent, IOperation& child);
		//void operationRemoved(IOperation& parent, IOperation& child);
		
	private:
		OperationGraph(OperationGraph const&);
		OperationGraph& operator=(OperationGraph const&);
		
		ICompressionPolicy* _policy;
	};
}
#endif
