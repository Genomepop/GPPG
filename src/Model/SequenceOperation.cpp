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
//using namespace GPPG; 

/*******************************************************************
 *				SEQUENCE OPERATION
 */
/*
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
*/

SequenceData* SequenceData::copy() const
{
	SequenceData* other = new SequenceData();
	other->length = length;
	other->sequence = (STYPE*)malloc(sizeof(STYPE)*length);
	memcpy( other->sequence, sequence, sizeof(STYPE)*length);
	return other;
}

SequenceOperationRoot::SequenceOperationRoot( SequenceData* data ) : Operation<SequenceData>(0) {
	setData( data );
}



SequencePointChange::SequencePointChange(Operation<SequenceData>& op, int* locs, int numLocs, STYPE* dest) : 
Operation<SequenceData>(numLocs, op), _loc(locs), _numlocs(numLocs), _c(dest) {}

SequencePointChange::~SequencePointChange() { delete _loc; delete _c; }

SequenceData* SequencePointChange::evaluate() const {
	
	SequenceData* sd = Operation<SequenceData>::evaluate();
	if (sd) return sd;
	
	// Get the sequence from the parent and add the point changes
	sd = parent(0)->evaluate();
	for (int i=0; i<_numlocs; i++) {
		sd->sequence[ _loc[i] ] = _c[i];
	}
	return sd;
}
