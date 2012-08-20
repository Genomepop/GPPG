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

SequencePointMutator::SequencePointMutator(double rate, const ublas::matrix<double> &T) : 
OperationMutator<SequenceData>(), _rate(rate), _M(T) {
	// Create random discrete distributions for each character
	
	std::vector<double> weights( _M.size1() );
	
	for (int i=0; i<_M.size1(); i++) {
		for (int j=0; j<_M.size2(); j++) {
			weights[j] = _M(i,j);
		}
		_transition.push_back( boost::random::discrete_distribution<>( weights ) );
	}
}

int SequencePointMutator::numMutants(IGenotype& g, long N, double f) const {
	return 1;
}

Operation<SequenceData>* SequencePointMutator::mutate( Operation<SequenceData>& g) const {
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
