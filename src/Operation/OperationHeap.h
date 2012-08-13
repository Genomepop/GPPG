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

// This needs to be thread safe, since it is a singleton
// But I won't do that in v1.

#include "CompressionPolicy.h"
#include "Operation.h"

namespace GPPG {
	template <typename T> class OperationHeap {
	public:
		static OperationHeap<T>& Instance();
		static OperationHeap<T>* InstancePtr();
		
		/**
		 Sets the policy for the @ref OperationHeap.  If there is already a policy in place, then this throws an error.
		 @param cp	The @ref CompressionPolicy.
		 */
		static void setCompressionPolicy(ICompressionPolicy& cp);
		
		/**
		 Retrieves the policy
		 */
		static ICompressionPolicy& compressionPolicy();
		
		void setPoolSize(int k);
		int poolSize() const;
		
		void operationAttached(IOperation& parent, IOperation& child);
		void operationRemoved(IOperation& parent, IOperation& child);
		
		T* malloc();
		void free(T* chunk);
		
	private:
		OperationHeap<T>();
		OperationHeap<T>(OperationHeap<T> const&);
		OperationHeap<T>& operator=(OperationHeap<T> const&);
		static OperationHeap<T>* m_pInstance;
		
		ICompressionPolicy* policy;
	};
}
#endif
