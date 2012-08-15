/*
 *  GenotypeFactory.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/15/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */
#ifndef GENOTYPE_FACTORY_
#define GENOTYPE_FACTORY_

#include "Genotype.h"

namespace GPPG {

/** Interface for the GenotypeFactory.
 * A GenotypeFactory generates random Genotypes.
 */
class IGenotypeFactory {
public:
	virtual ~IGenotypeFactory(){}
	
	/** Create a random genotype.
	 * Implementing classes use this function to produce a random genotype.
	 */
	virtual IGenotype* random() const = 0;


};

template <typename T>
class GenotypeFactory : public IGenotypeFactory {
public:
	virtual Genotype<T>* random() const { return new Genotype<T>( randomData() ); }
	
	virtual T* randomData() const = 0;
};

}

#endif