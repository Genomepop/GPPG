/*
 *  Genotype.h
 *  GPPG
 *
 *  Created by Troy Ruths on 8/10/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef GENOTYPE_
#define GENOTYPE_

namespace GPPG {
	class IGenotype {
	public:
		virtual ~IGenotype(){}
		
		virtual void configure() = 0;
		
		virtual double frequency() const = 0;
		virtual void setFrequency(double freq) = 0;
		
		virtual double total() const = 0;
		virtual void setTotal(double total) = 0;
		
		//virtual bool operator==(IGenotype const& other) const;
		
	};
	
	/*
	 * The Genotype class ... should the cache be here?
	 * Subclass this genotype class or use the Template parameter to provide your own content.
	 * The template is the GenotypeData (e.g. DNA, pathway, etc.)
	 */
	template <typename T> class Genotype : public IGenotype {
		
	public:
		Genotype<T>(T* genoData);
		
		void configure();
		
		double frequency() const;
		void setFrequency(double f);
		
		double total() const;
		void setTotal(double t);
		
		virtual T* data() const;

	protected:
		void setData(T* data);		
	
	private:
		// Disable copy-construction
		Genotype<T>(Genotype<T> const& g) {}
		Genotype<T>& operator=(Genotype<T> const& g) {}
		
		T* _data;				
		double _freq, _total;
		
	};

}
#endif
