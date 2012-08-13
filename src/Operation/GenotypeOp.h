/*
 *  GenotypeOp.h
 *  GPPG
 *
 *  Created by Troy Ruths on 8/10/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef GENOTYPEOP_
#define GENOTYPEOP_

#include "Genotype.h"
#include "Operation.h"

namespace GPPG {
	class IGenotypeOp : public IGenotype {
	public:
		virtual ~IGenotypeOp(){}
		
		virtual IOperation& operation() const;
		
		virtual bool isCompressed() const = 0;
		virtual bool setCompressed(bool doCompress) = 0;
		
	};
	
	/*
	 * The Genotype class ... should the cache be here?
	 * Subclass this genotype class or use the Template parameter to provide your own content.
	 * The template is the GenotypeData (e.g. DNA, pathway, etc.)
	 */
	template <typename T> class GenotypeOp : public IGenotypeOp, Genotype<T> {
		
	public:
		GenotypeOp<T>(Operation<T> &operation);
		
		Operation<T>& operation() const;
		
		bool isCompressed() const;
		bool setCompressed(bool doCompress);
		
		T* data() const;
		
	protected:
		Operation<T> _operation;
		
	private:
		void innerConstructor();
		
		// Disable copy-construction
		GenotypeOp<T>(GenotypeOp<T> const& g);
		GenotypeOp<T>& operator=(GenotypeOp<T> const& g);
		
	};
	
}

#endif