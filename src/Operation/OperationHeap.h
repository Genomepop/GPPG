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
	
	template <typename T> class OperationGraph : IGenotypeHeap {
	public:
		
		/**
		 Sets the policy for the @ref OperationHeap.  If there is already a policy in place, then this throws an error.
		 @param cp	The @ref CompressionPolicy.
		 */
		static void setCompressionPolicy(ICompressionPolicy& cp);
		
		/**
		 Retrieves the policy
		 */
		static ICompressionPolicy& compressionPolicy();
		
		
		//void operationAttached(IOperation& parent, IOperation& child);
		//void operationRemoved(IOperation& parent, IOperation& child);
				
	private:
		OperationGraph<T>();
		OperationGraph<T>(OperationGraph<T> const&);
		OperationGraph<T>& operator=(OperationGraph<T> const&);
		
		ICompressionPolicy* policy;
	};
}
#endif
