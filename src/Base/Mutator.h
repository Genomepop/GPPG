/*
 *  Mutator.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/15/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef MUTATOR_
#define MUTATOR_

namespace GPPG {
	class IGenotype;
	
	class IMutator {
	public:
		~IMutator() {}
		
		virtual IGenotype* mutate( IGenotype& geno) const = 0;
	};
	
	
}
#endif MUTATOR_