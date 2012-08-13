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

#include "CompressionPolicy.h"

namespace GPPG {
	class GreedyLoad : public ICompressionPolicy {
	public:
		GreedyLoad(int k);
		
		void apply( );
		
	private:
		GreedyLoad(GreedyLoad const&);
		GreedyLoad const& operator=(GreedyLoad const&);
		
		int k;
	};
}
#endif
