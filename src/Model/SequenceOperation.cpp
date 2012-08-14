/*
 *  SequenceOperation.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 8/13/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "SequenceOperation.h"

using namespace GPPG::Model;

/*******************************************************************
 *				SEQUENCE OPERATION
 */
SequenceOperation::SequenceOperation( double cost, int length ) : 
	Operation<SequenceData>(cost), _length(length) {}

SequenceOperation::SequenceOperation( double cost, int length, SequenceOperation& parent1 ) : 
	Operation<SequenceData>(cost, parent1), _length(length) {}

SequenceOperation::SequenceOperation( double cost, int length, SequenceOperation& parent1, SequenceOperation& parent2 ) :
	Operation<SequenceData>(cost, parent1, parent2), _length(length) {}

STYPE SequenceOperation::get(int i) const { 
	if (isCompressed()) {
		if (numParents() > 0) {
			return ((SequenceOperation&)parent(0)).get(i);
		}
		throw "No Parent to propogate get() request";
	}
	return data()->sequence[i];
}

int SequenceOperation::length() const {
	return _length; 
}

/*******************************************************************
 *				SEQUENCE DATA ROOT
 */
SequenceOperationRoot::SequenceOperationRoot( SequenceData* data ) : SequenceOperation(0, data->length) {
	setData( data );
}

