/*
 *  BaseCompressionPolicy.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/23/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef BASE_COMPRESSION_POLICY_
#define BASE_COMPRESSION_POLICY_

#include "Operation/CompressionPolicy.h"

namespace GPPG {
	
	/** StoreFlag determines the simple storage policy.
	 */
	enum StoreFlag { 
		STORE_ROOT,		/* Store only the root as uncompressed */
		STORE_ALL,		/* Store all operations uncompressed */
		STORE_ACTIVE		/* Store active operations */
	};
	
	class BaseCompressionPolicy : public CompressionPolicy {
	public:
		BaseCompressionPolicy(StoreFlag flag);
		
		~BaseCompressionPolicy();
		
		void operationAdded(IOperation* op);
		
		void operationRemoved(IOperation* op);
		
		/** Gets the StoreFlag associated with this policy.
		 */
		StoreFlag storeFlag() const;
		
		void generationFinished( OperationGraph *heap, const std::set<IOperation*>& active );
		
	private:
		StoreFlag _flag;
		
	};
}
#endif
