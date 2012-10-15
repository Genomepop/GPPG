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

#include <vector>
#include <set>

namespace GPPG {
	class IOperation;
	
	class ICompressionPolicy {
	public:
		virtual ~ICompressionPolicy(){}
		
		/** Called when an Operation is added to the OperationGraph
		 */
		virtual void operationAdded( IOperation* op ) = 0;
		
		/** Called when an Operation is removed from the OperationGraph
		 */
		virtual void operationRemoved( IOperation* op ) = 0;
		
		/** Called when an uncompressed Operation is deleted by the OperationGraph
		 */
		virtual void decompressionReleased( IOperation* op ) = 0;
		
		/** Called when a generation, or a time period equivalent to a generation, is complete.
		 * This provides the policy an opportunity to optimize the compression decisions.
		 */
		virtual void generationFinished( const std::set<IOperation*>& ) = 0;

		virtual void generationFinished( const std::vector<IOperation*>& ) = 0;
	};
	
	/** This class provides a bare-bones implementation of the CompressionPolicy
	 */
	class CompressionPolicy : public ICompressionPolicy {
	public:
		void operationAdded(IOperation* op);
		
		void operationRemoved( IOperation* op );
		
		void decompressionReleased( IOperation* op );
		
		void generationFinished( const std::vector<IOperation*>& );
		
		void generationFinished( const std::set<IOperation*>& );
	};
}

#endif