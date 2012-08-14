/*
 *  CompressionPolicy.h
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef CORE_COMPRESSION_POLICY_
#define CORE_COMPRESSION_POLICY_

//#include "Operation.h"

namespace GPPG {

	class ICompressionPolicy {
	public:
		virtual ~ICompressionPolicy(){}
		
		
		//virtual void operationAdded( IOperation& op ) = 0;
		
		//virtual void operationRemoved( IOperation& op ) = 0;
		
		
		//virtual void compressionReleased( IOperation& op );
		
		//void decompressOperation(IOperation& op);
		//void compressOperation(IOperation& op);
		
	};
	
}

#endif