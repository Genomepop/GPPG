/*
 *  Operation.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 8/17/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "Operation.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/binomial_distribution.hpp>
//#include <boost/math/distributions/binomial.hpp>

// Seed for RNG
boost::mt19937 gen;

using namespace GPPG::Model;
using namespace GPPG; 

/*******************************************************************
 *				SEQUENCE OPERATION
 */
/*
 SequenceOperation::SequenceOperation( double cost, int length ) : 
 SequenceOperation(cost), _length(length) {}
 
 SequenceOperation::SequenceOperation( double cost, int length, SequenceOperation& parent1 ) : 
 SequenceOperation(cost, parent1), _length(length) {}
 
 SequenceOperation::SequenceOperation( double cost, int length, SequenceOperation& parent1, SequenceOperation& parent2 ) :
 SequenceOperation(cost, parent1, parent2), _length(length) {}
 
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

OpSequenceBase::OpSequenceBase( double cost, int length, OpSequence& parent1 ):
OpSequence(cost, parent1), _length(length) {}

OpSequenceBase::OpSequenceBase( double cost, int length, OpSequence& parent1, OpSequence& parent2 ):
OpSequence(cost, parent1, parent2), _length(length) {}

int OpSequenceBase::length() const { return _length; }

STYPE OpSequenceBase::get(int i) const {
	if (isCompressed()) {
		return proxyGet(i);
	}
	return data()->get(i);
}

SequenceRoot::SequenceRoot(SequenceData* d) : OperationRoot<SequenceData, ISequence>(d) {}

int SequenceRoot::length() const { return data()->length(); }
STYPE SequenceRoot::get(int i) const { return data()->get(i); }


SequenceRootFactory::SequenceRootFactory(int length, const ublas::vector<double>& distr ) : _length(length), _distr(distr) {}

SequenceRoot* SequenceRootFactory::random() const {
	SequenceData* sd = new SequenceData(_length);
	boost::random::discrete_distribution<> dist(_distr);
	gen.seed((unsigned int)time(0));
	for (int i=0; i<_length; i++) {
		STYPE c = (STYPE)dist(gen);
		sd->set(i, c);
	}
	return new SequenceRoot(sd);
}

//SequenceOperationRoot::SequenceOperationRoot( SequenceData* data ) : SequenceOperation(0) {
//	setData( data );
//}


SequencePointChange::SequencePointChange(OpSequence& op, int* locs, int numLocs, STYPE* dest) : 
OpSequenceBase(numLocs,op.length(), op), _loc(locs), _numlocs(numLocs), _c(dest) {}

SequencePointChange::~SequencePointChange() { delete _loc; delete _c; }

SequenceData* SequencePointChange::evaluate() const {
	SequenceData* sd = OpSequence::evaluate();
	if (sd != NULL) return sd;
	
	// Get the sequence from the parent and add the point changes
	sd = parent(0)->evaluate();
	
	for (int i=0; i<_numlocs; i++) {
		sd->set(_loc[i], _c[i]);
	}
	return sd;
}

std::string SequencePointChange::toString() const {
	std::ostringstream output;
	for (int i=0; i<numSites(); i++) {
		output << getSite(i) << "->" << getMutation(i);
		if (i < numSites()-1) output << ", ";
	}
	return output.str();
}

int SequencePointChange::numSites() const { return _numlocs; }

STYPE SequencePointChange::getMutation(int i) const { return _c[i]; }

int SequencePointChange::getSite(int i) const { return _loc[i]; }


STYPE SequencePointChange::proxyGet(int l) const {
	// See if the index is in the list
	for (int i=0; i<_numlocs; i++) {
		if (_loc[i] == l) return _c[i];
	}
	return parent(0)->get(l);
}

SequencePointMutator::SequencePointMutator(double rate, const ublas::matrix<double> &T) : 
OperationMutator<OpSequence>(), _rate(rate), _M(T) {
	// Create random discrete distributions for each character
	
	std::vector<double> weights( _M.size1() );
	
	for (int i=0; i<_M.size1(); i++) {
		for (int j=0; j<_M.size2(); j++) {
			weights[j] = _M(i,j);
		}
		_transition.push_back( boost::random::discrete_distribution<>( weights ) );
	}
}



OpSequence* SequencePointMutator::mutate( OpSequence& g) const {
	// Get the sequence data
	int isCopy;
	SequenceData* data = g.data(isCopy);
	
	// Calculate the number of sites to mutate (use binomial)
	boost::random::binomial_distribution<> dist( data->length(), _rate );
	
	int numLocs = dist(gen);
	if (numLocs == 0) numLocs = 1;
	boost::random::uniform_int_distribution<> idist( 0, data->length()-1 );	
	
	int* locs = (int*) malloc(sizeof(int)*numLocs);
	STYPE* dest = (STYPE*) malloc(sizeof(STYPE)*numLocs);
	for (int i=0; i<numLocs; i++) {
		int loc = idist(gen);
		STYPE c = data->get(loc);
		locs[i] = loc;
		dest[i] = (STYPE)( _transition[c](gen) );
		
	}
	SequencePointChange* spc = new SequencePointChange(g, locs, numLocs, dest);
	
	if (isCopy) delete data;
	
	return spc;
}

double SequencePointMutator::rate() const { return _rate; }

const ublas::matrix<double>& SequencePointMutator::transition() const { return _M; }

int SequencePointMutator::numMutants(OpSequence& g, long N, double f) const {
	boost::random::binomial_distribution<> dist( N*f, _rate*g.length() );
	return dist(gen);
}

ostream& operator<<(ostream& output, const SequencePointMutator& s) {
	output << "SequencePointMutator(rate=" << s.rate() << "): " << s.transition();
	return output;
}

/*
 ostream& operator<<(ostream& output, const SequencePointChange& s) {
 output << "Point Change(s):\n ";
 for (int i=0; i<s.numSites(); i++) {
 output << "[" << i << "] " << s.getSite(i) << " -> " << s.getMutation(i) << "\n";
 }
 return output;
 }
 
 
 ostream& operator<<(ostream& output, const GPPG::Model::SequenceOperation& s) {
 SequenceData* sd = s.evaluate();
 output << "SequenceOperation [" << typeid(s).name() << "]: " << *sd;
 delete sd;
 
 return output;
 }
 */
