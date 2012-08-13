/*
 *  Operation.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */


#include "Operation.h"

using namespace GPPG;

template <typename T> Operation<T>::Operation(double cost): cache(0), m_cost(cost), parent1(0), parent2(0) 
{ innerConstructor(); }

template <typename T> Operation<T>::Operation(double cost, Operation<T>& parent): 
cache(0), m_cost(cost), parent1(&parent), parent2(0) { innerConstructor(); }

template <typename T> Operation<T>::Operation(double cost, Operation<T>& parent, Operation<T>& parent2): 
cache(0), m_cost(cost), parent1(&parent),	 parent2(&parent2) { innerConstructor(); }

template <typename T> void Operation<T>::innerConstructor() {
	if (parent1 != 0) {
		parent1.addChild( this );
	}
	if (parent2 != 0) {
		parent2.addChild( this );
	}
	
	m_load = 0;
	cache = NULL;
}


template <typename T> Operation<T>& Operation<T>::parent(int i) const {
	if( i<0 || i>=numParents()) { throw "Incorrect index"; }
	
	switch( i ){
	case 0: return parent1; 
	case 1: return parent2;
	}
	
	throw "Index is too high";
}

template <typename T> int Operation<T>::numParents() const {
	if (parent2 != 0) {
		return 2;
	} else if (parent1 != 0) {
		return 1;
	}
	return 0;
}

template <typename T> T* Operation<T>::result() const {
	if (cache != NULL) return cache;
	return evaluate();
}
		
template <typename T> double Operation<T>::cost() const { return m_cost; }

template <typename T> void Operation<T>::setCompressed(bool compress) {
	if (compress && !isCompressed()) {
		// TODO: deallocate the cache
	} else if (!compress && isCompressed()) {
		// Fill the cache
		cache = evaluate();
	}
}

template <typename T> bool Operation<T>::isCompressed() const {
	return cache == NULL;
}