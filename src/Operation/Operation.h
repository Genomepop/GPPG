/*
 *  Operation.h
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef CORE_OPERATION_
#define CORE_OPERATION_

#include <set>

namespace GPPG {
	
	class IOperation {
	public:
		virtual ~IOperation(){}
		
		virtual void setCompressed(bool compress) = 0;
		virtual bool isCompressed() = 0;
		
		// Should this be a part of the interface?
		virtual void setLoad(double value) = 0;
		virtual double load() const = 0;
		
		/**
		 * Returns whether or not the genotype of this operation is active
		 */
		//virtual bool isActive() const = 0;
		
		virtual double cost() const = 0;
		
		virtual int numChildren() const = 0;
		
		virtual IOperation& parent(int i) = 0;
		virtual int numParents() const = 0;
		
		//virtual void addChild(IOperation& op) = 0;
	};
	
	template <typename T> class Operation : public IOperation {
	public:
		Operation<T>(double cost);
		Operation<T>(double cost, Operation<T> &parent);
		Operation<T>(double cost, Operation<T> &parent1, Operation<T> &parent2);
		
		~Operation<T>();
		
		Operation<T>& parent(int i) const;
		int numParents() const;
		
		/**
		 Returns the result of this operation. If the cache is full, then it will be returned quickly; 
		 otherwise, the operation will be evaluated, along with its parents.
		 */
		T* result() const;
		
		/**
		 Returns the number of children
		 */
		int numChildren() const;
		//set<Operation<T>& >& children() const;
		
		// Manage compression
		void setCompressed(bool compress);
		bool isCompressed() const;
		
		// Manage Load
		void setLoad(double value);
		double load() const;
		
		// Data Size
		int dataSize() const;
		
		/**
		 * The cost of applying the operation
		 */
		double cost() const;
		
	protected:
		virtual T* evaluate() const = 0;
		
	private:
		Operation<T>(Operation<T> const& op) {}
		Operation<T>& operator=(Operation<T> const& op) {}
		
		void inline innerConstructor();
		void addChild(Operation<T>& op);
		//void setParent(Operation<T>& op);
		//void setParent(int i, Operation<T>& op);
		
		T *cache;
		double m_load;
		double m_cost;
		bool compress;
		Operation<T> *parent1, *parent2;
		//set<Operation<T>& > operations;
	};
}

#endif
