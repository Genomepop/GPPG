/*
 *  GenotypeHeap.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/23/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef GENOTYPE_HEAP_
#define GENOTYPE_HEAP_

#include <vector>
#include <set>

namespace GPPG {

	class IGenotype;
	
enum HeapRemovalPolicy {
	DELETE, KEEP
};
	
class IGenotypeHeap {
public:
	virtual ~IGenotypeHeap() {};
	
	virtual void setHeapRemovalPolicy(HeapRemovalPolicy policy) = 0;
	virtual HeapRemovalPolicy heapRemovalPolicy() const = 0;
	
	/** Calling this function gives ownership of the \param g to the heap.  
	 * The lifespan of g is now tied to this heap.
	 */
	virtual void addGenotype(IGenotype* g) = 0;
	
	/** The simulator uses this function to notify the Heap that the genotype may be removed.
	 */
	virtual void removeGenotype(IGenotype* g) = 0;
	
	virtual void generationFinished(const std::vector<IGenotype*>&) = 0;
	virtual void generationFinished(const std::set<IGenotype*>&) = 0;
	
};

class BasicGenotypeHeap : public IGenotypeHeap {
public:
	BasicGenotypeHeap();
	
	~BasicGenotypeHeap();
	
	void setHeapRemovalPolicy( HeapRemovalPolicy policy);
	
	HeapRemovalPolicy heapRemovalPolicy() const;
	
	void addGenotype(IGenotype* g);
	
	void removeGenotype(IGenotype* g);
	
	void generationFinished(const std::vector<IGenotype*>&);
};
	
}
#endif