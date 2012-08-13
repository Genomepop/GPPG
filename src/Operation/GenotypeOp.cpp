/*
 *  GenotypeOp.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 8/10/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "GenotypeOp.h"

using namespace GPPG;

template <typename T> GenotypeOp<T>::GenotypeOp(Operation<T> &operation): 
		Genotype<T>(0), _operation(operation)
{ innerConstructor(); }

template <typename T> Operation<T>& GenotypeOp<T>::operation() const {
	return _operation;
}

template <typename T> void GenotypeOp<T>::innerConstructor() {
	// Add the operation into the Genotype relationships
	// Here I assume _operation has been set by a constructor
	
}

template <typename T> bool GenotypeOp<T>::isCompressed() const {
	return this._operation.isCompressed();
}

template <typename T> bool GenotypeOp<T>::setCompressed(bool doCompress) {
	this._operation.setCompressed(doCompress);
	/*
	if (doCompress) {
		// Clear data
		if (this._data != 0) {
			// Free the data
			delete this._data;
			this._data = 0;
		}
	} else {
		// Uncompress data
		if (this._data == 0) {
			this._data = this.op.result();
		}
		
	}*/
	
	return true;
}

template <typename T> T* GenotypeOp<T>::data() const {
	return this._operation.result();
	/*
	if (this._data == 0) {
		return this._operation.result();
	} 
	return this._data;
	 */
}