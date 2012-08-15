/*
 *  Operation.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

// TODO: Move this template code into the header!

#include "Operation.h"

using namespace GPPG;

/*
template <typename T> Operation<T>::Operation(double cost): 
Genotype<T>(0) 
{ innerConstructor(cost, 0, 0); }

template <typename T> Operation<T>::Operation(double cost, Operation<T>& parent): 
Genotype<T>(0)
{ innerConstructor(cost, &parent, 0); }

template <typename T> Operation<T>::Operation(double cost, Operation<T>& parent, Operation<T>& parent2): 
Genotype<T>(0)
{ innerConstructor(cost, &parent, &parent2); }

template <typename T>
Operation<T>::~Operation() {
	// TODO: What to do about children and parents?

}

template <typename T> void Operation<T>::innerConstructor(double cost, Operation<T>* parent1, Operation<T>* parent2) {
	_cost = cost;
	_parent1 = parent1;
	_parent2 = parent2;
	
	if (_parent1 != 0) {
		_parent1.addChild( this );
	}
	if (_parent2 != 0) {
		_parent2.addChild( this );
	}
	
	_load = 0;
	//_genotype = new GenotypeOp<T>( this );
}


template <typename T> Operation<T>& Operation<T>::parent(int i) const {
	if( i<0 || i>=numParents()) { throw "Incorrect index"; }
	
	switch( i ){
	case 0: return _parent1; 
	case 1: return _parent2;
	}
	
	throw "Index is too high";
}

template <typename T> int Operation<T>::numParents() const {
	if (_parent2 != 0) {
		return 2;
	} else if (_parent1 != 0) {
		return 1;
	}
	return 0;
}

template <typename T> T* Operation<T>::data() const {
	if (isCompressed()) {
		return evaluate();
	}
	return Genotype<T>::data();
}

template <typename T>
void Operation<T>::addChild(Operation<T>* op) {
	_children.push_back( op );
}

template <typename T>
void Operation<T>::removeChild(Operation<T>* op) {
	_children.remove( op );
}

template <typename T>
int Operation<T>::numChildren() const {
	return _children.size();
}

template <typename T>
int Operation<T>::dataSize() const {
	return 1;
}

//template <typename T> GenotypeOp<T>& Operation<T>::genotype() const {
//	return _genotype;
//}

template <typename T> 
double Operation<T>::cost() const { return _cost; }

template <typename T> 
void Operation<T>::setCompressed(bool compress) {
	if (compress && !isCompressed()) {
		Genotype<T>::setData(0);
	} else if (!compress && isCompressed()) {
		// Fill the cache
		Genotype<T>::setData( evaluate() );
	}
}

template <typename T>
T* Operation<T>::evaluate() const { 
	if (isCompressed()) {
		return NULL;
	}
	return Genotype<T>::data()->copy();
}

template <typename T> bool Operation<T>::isCompressed() const {
	return data() == NULL;
}
*/