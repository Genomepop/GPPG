/*
 *  SequenceOperation.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 8/13/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "SequenceOperation.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/numeric/ublas/io.hpp>

// Seed for RNG
boost::mt19937 gen;

using namespace GPPG::Model;
using namespace GPPG; 

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
	SequenceData* other = new SequenceData(_length);
	memcpy( other->sequence(), _sequence, sizeof(STYPE)*_length);
	return other;
}

SequenceData::SequenceData(int length) : _length(length) {
	_sequence = (STYPE*)malloc(sizeof(STYPE)*_length);
}

STYPE* SequenceData::sequence() { return _sequence; }
int SequenceData::length() const { return _length; }

STYPE SequenceData::get(int i) const { return _sequence[i]; }

void SequenceData::set(int i, STYPE c) { _sequence[i] = c; }

ostream& operator<<(ostream& output, const SequenceData& s) {
	output << "[";
	for (int i=0; i<s.length(); i++) output << s.get(i);
	output << "]";
	return output;
}

ostream& operator<<(ostream& output, const SequencePointMutator& s) {
	output << "SequencePointMutator(rate=" << s.rate() << "): " << s.transition();
	return output;
}

SequenceFactory::SequenceFactory(int length, const ublas::vector<double>& distr ) : _length(length), _distr(distr) {}

SequenceData* SequenceFactory::randomData() const {
	SequenceData* sd = new SequenceData(_length);
	boost::random::discrete_distribution<> dist(_distr);
	gen.seed((unsigned int)time(0));
	for (int i=0; i<_length; i++) {
		STYPE c = (STYPE)dist(gen);
		sd->set(i, c);
	}
	return sd;
}

//SequenceOperationRoot::SequenceOperationRoot( SequenceData* data ) : Operation<SequenceData>(0) {
//	setData( data );
//}



SequencePointChange::SequencePointChange(Operation<SequenceData>& op, int* locs, int numLocs, STYPE* dest) : 
Operation<SequenceData>(numLocs, op), _loc(locs), _numlocs(numLocs), _c(dest) {}

SequencePointChange::~SequencePointChange() { delete _loc; delete _c; }

SequenceData* SequencePointChange::evaluate() const {
	
	SequenceData* sd = Operation<SequenceData>::evaluate();
	if (sd) return sd;
	
	// Get the sequence from the parent and add the point changes
	sd = parent(0)->evaluate();

	for (int i=0; i<_numlocs; i++) {
		sd->set(_loc[i], _c[i]);
	}
	return sd;
}

SequencePointMutator::SequencePointMutator(double rate, const ublas::matrix<double> &T) : 
	OperationMutator<SequenceData>(), _rate(rate), _M(T) {}

Operation<SequenceData>* SequencePointMutator::mutate( Operation<SequenceData>& g) const {
	
	int* locs;
	int numLocs = 3;
	STYPE* dest;
	SequencePointChange* spc = new SequencePointChange(g, locs, 5, dest);
	return spc;
}

double SequencePointMutator::rate() const { return _rate; }

const ublas::matrix<double>& SequencePointMutator::transition() const { return _M; }


