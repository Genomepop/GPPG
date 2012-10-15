/*
 *  Recombinator.h
 *  Demo
 *
 *  Created by Troy Ruths on 9/14/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef RECOMBINATOR_
#define RECOMBINATOR_

namespace GPPG {
	class IGenotype;
	
	class IRecombinator {
	public:
		IRecombinator() {}
		
		virtual int numMutants( IGenotype& geno1, IGenotype& geno2, long N) const = 0;
		
		virtual IGenotype* recombine( IGenotype& geno1, IGenotype& geno2) const = 0;
	};
	
}

#endif