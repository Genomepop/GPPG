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
		
		virtual bool isActive() const = 0;
		
		virtual int index() const = 0;
		virtual void setIndex(int i) = 0;
		
		virtual int order() const = 0;
		virtual void setOrder(int i) = 0;
	};
	
	class BaseGenotype : public IGenotype {
	public:
		BaseGenotype();
		
		void configure();
		double frequency() const;
		void setFrequency(double freq);
		
		double total() const;
		void setTotal(double total);
		
		bool isActive() const;
		
		int index() const;
		void setIndex(int i);
		
		int order() const;
		void setOrder(int i);
		
	private:
		double _freq, _total;
		int _index, _order;
	};
	
	/*
	 * The Genotype class ... should the cache be here?
	 * Subclass this genotype class or use the Template parameter to provide your own content.
	 * The template is the GenotypeData (e.g. DNA, pathway, etc.)
	 */
	template <typename T> class Genotype : public BaseGenotype {
		
	public:
		Genotype<T>(T* genoData) : BaseGenotype(), _data(genoData) {}
		
		
		virtual T* data() const { return _data; }

	protected:
		void setData(T* data) { _data = data; }
	
	private:
		// Disable copy-construction
		Genotype<T>(Genotype<T> const& g) {}
		Genotype<T>& operator=(Genotype<T> const& g) {}
		
		T* _data;				
	};

}
#endif
