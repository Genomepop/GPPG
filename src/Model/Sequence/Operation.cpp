/*
 *  Operation.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 8/17/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "GPPG.h"
#include "Operation.h"
//#include "Util/Random.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/binomial_distribution.hpp>
//#include <boost/math/distributions/binomial.hpp>
#include <algorithm>

// Seed for RNG
boost::mt19937 gen;

using namespace GPPG::Model;
using namespace GPPG; 

/*******************************************************************
 *				OPSEQUENCEBASE OPERATION
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


SequenceData* randomSequenceData(int length, const ublas::vector<double>& distr) {
	SequenceData* sd = new SequenceData(length);
	boost::random::discrete_distribution<> dist(distr);
	gen.seed((unsigned int)time(0));
	for (int i=0; i<length; i++) {
		STYPE c = (STYPE)dist(gen);
		sd->set(i, c);
	}
	return sd;
}

SequenceRootFactory::SequenceRootFactory(int length, const ublas::vector<double>& distr ) : _length(length), _distr(distr) {}

SequenceRoot* SequenceRootFactory::random() const {
	
	return new SequenceRoot( randomSequenceData( _length, _distr) );
}



SequencePointChange::SequencePointChange(OpSequence& op, int* locs, int numLocs, STYPE* dest) : 
OpSequenceBase(numLocs,op.length(), op), _loc(locs), _numlocs(numLocs), _c(dest) {}

SequencePointChange::~SequencePointChange() { delete _loc; delete _c; }

SequenceData* SequencePointChange::evaluate() const {
	SequenceData* sd = OpSequenceBase::evaluate();
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
	//int isCopy;
	//SequenceData* data = g.data(isCopy);
	
	// Calculate the number of sites to mutate (use binomial)
	int length = g.length(); //data->length()
	boost::random::binomial_distribution<> dist( length, _rate );
	
	int numLocs = dist(gen);
	//if (numLocs == 0) return &g; 
	
	if (numLocs == 0) numLocs = 1;
	boost::random::uniform_int_distribution<> idist( 0, length-1 );	
	
#ifdef DEBUG_0
	std::cout << "SequencePointMutator: mutating..." << std::endl;
#endif
	
	int* locs = (int*) malloc(sizeof(int)*numLocs);
	STYPE* dest = (STYPE*) malloc(sizeof(STYPE)*numLocs);
	for (int i=0; i<numLocs; i++) {
		int loc = idist(gen);
		STYPE c = g.get(loc); //data->get(loc);
		locs[i] = loc;
		dest[i] = (STYPE)( _transition[c](gen) );
		
	}
	SequencePointChange* spc = new SequencePointChange(g, locs, numLocs, dest);
	
#ifdef DEBUG_0
	std::cout << "Created Operation (" << spc->numParents() << "): " << spc->toString() << std::endl;
#endif
	//if (isCopy) delete data;
	
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


SequenceDeletion::SequenceDeletion(OpSequence& op, int loc, int span): OpSequenceBase(1, op.length()-span, op), _loc(loc), _span(span) {}

SequenceData* SequenceDeletion::evaluate() const {
	SequenceData* sd = OpSequenceBase::evaluate();
	if (sd != NULL) return sd;
	
	// Get the sequence from the parent and add the point changes
	sd = parent(0)->evaluate();
	
	
	SequenceData* data = new SequenceData( length() );
	
	STYPE* sdata = data->sequence();
	STYPE* pdata = sd->sequence();
	if (_loc > 0) 
		memcpy(sdata, pdata, sizeof(STYPE)*_loc);
	memcpy(&sdata[_loc], &pdata[_loc+_span] , sizeof(STYPE)*(_span-_loc));

	delete sd;
	
	return data;
}

STYPE SequenceDeletion::proxyGet(int i) const {
	// See if the index is in the list
	if (i >= _loc) {
		return parent(0)->get(i+_span);
	}
	return parent(0)->get(i);
}

SequenceInsertion::SequenceInsertion(OpSequence& op, int loc, SequenceData* span): 
	OpSequenceBase(1, op.length()+span->length(), op), _loc(loc), _span(span) {}

SequenceInsertion::~SequenceInsertion() { delete _span; }

SequenceData* SequenceInsertion::evaluate() const {
	SequenceData* sd = OpSequenceBase::evaluate();
	if (sd != NULL) return sd;
	
	// Get the sequence from the parent and add the point changes
	sd = parent(0)->evaluate();
	
	
	SequenceData* data = new SequenceData( length() );
	int end = _loc+_span->length();
	
	STYPE* sdata = data->sequence();
	STYPE* pdata = sd->sequence();
	if (_loc > 0) 
		memcpy(sdata, pdata, sizeof(STYPE)*_loc);
	// copy span
	memcpy(&sdata[_loc], _span->sequence() , sizeof(STYPE)*(_span->length()));
	
	memcpy(&sdata[end], &pdata[_loc] , sizeof(STYPE)*(length()-end));
	delete sd;
	
	return data;
}

STYPE SequenceInsertion::proxyGet(int i) const {
	// See if the index is in the list
	int index = i;
	if (i >= _loc+_span->length()) {
		index = i-_span->length();
	} else if (i >= _loc && i < _loc+_span->length()) {
		return _span->get( i-_loc );
	}
	
	return parent(0)->get(index);
}

SequenceDeletionMutator::SequenceDeletionMutator(double rate, int minL, int maxL) :
OperationMutator<OpSequence>(), _rate(rate), _minL(minL), _maxL(maxL) {}

OpSequence* SequenceDeletionMutator::mutate( OpSequence& g) const {
	// Perform deletion	
	int length = g.length(); 
	
	boost::random::uniform_int_distribution<> dist( _minL, _maxL );
	int spanLength = dist(gen);
	
	boost::random::uniform_int_distribution<> idist( 0, length-1-spanLength );	
	int loc = idist(gen);
	
	SequenceDeletion *sd = new SequenceDeletion(g, loc, spanLength);
	return sd;
	
}

int SequenceDeletionMutator::numMutants(OpSequence& g, long N, double f) const {
	boost::random::binomial_distribution<> dist( N*f, _rate*g.length() );
	return dist(gen);
}

double SequenceDeletionMutator::rate() const { return _rate; }


SequenceInsertionMutator::SequenceInsertionMutator(double rate, int minL, int maxL, const ublas::vector<double>& distr) :
OperationMutator<OpSequence>(), _rate(rate), _minL(minL), _maxL(maxL), _distr(distr) {}

OpSequence* SequenceInsertionMutator::mutate( OpSequence& g) const {
	// Perform insertion
	int length = g.length(); 
	
	boost::random::uniform_int_distribution<> dist( _minL, _maxL );
	int spanLength = dist(gen);
	
	boost::random::uniform_int_distribution<> idist( 0, length-1 );	
	int loc = idist(gen);
	
	// Generate random sequence
	SequenceData* span = randomSequenceData(spanLength, _distr);
	
	SequenceInsertion *sd = new SequenceInsertion(g, loc, span);
	return sd;
}

int SequenceInsertionMutator::numMutants(OpSequence& g, long N, double f) const {
	boost::random::binomial_distribution<> dist( N*f, _rate*g.length() );
	return dist(gen);
}

double SequenceInsertionMutator::rate() const { return _rate; }


/*******************************************************************
 * SEQUENCE RECOMBINATION
 *******************************************************************/

SequenceCrossover::SequenceCrossover(OpSequence& op1, OpSequence& op2, const std::vector<int>& locs) :
OpSequenceBase(100, (locs.size()%2==0)? op1.length(): op2.length(), op1, op2), _locs(locs) {}


SequenceData* SequenceCrossover::evaluate() const {
	SequenceData* sd = OpSequenceBase::evaluate();
	if (sd != NULL) return sd;
	
	SequenceData* sd1 = parent(0)->evaluate();
	SequenceData* sd2 = parent(1)->evaluate();
	
	// choose the host
	SequenceData* result = (_locs.size() % 2 == 0) ? sd1 : sd2;
	SequenceData* to_delete = (_locs.size() % 2 == 0) ? sd2 : sd1;
	STYPE* data = result->sequence();
	STYPE* other = to_delete->sequence();
	
	int a,b,l;
	int start = (_locs.size() % 2 == 0) ? 1: 0;
	
	for (int i=start; i<_locs.size(); i+=2) {
		a = (i==0) ? 0 : _locs[i-1];
		b = _locs[i];	
		l = b-a;
		memcpy( &data[a], &other[a], l*sizeof(STYPE) );
	}
	
	delete to_delete;
	return result;
}

STYPE SequenceCrossover::proxyGet(int i) const {
	// See if the index is in the list
	for (int j=0; j<_locs.size(); j++) {
		if (i<j) {
			return (j%2==0) ? parent(0)->get(i) : parent(1)->get(i);
		}
	}
	return (_locs.size()%2==0) ? parent(0)->get(i) : parent(1)->get(i);
}

SequenceRecombinator::SequenceRecombinator(double rate) : _rate(rate) {}

int SequenceRecombinator::numMutants(OpSequence& g, OpSequence& g2, long N) const {
	double amt = N*g.frequency()*g2.frequency();
	if (amt <= 0) return 0;
	
	else if (amt < 1) {
		boost::random::uniform_01<> r;
		if (r(gen) < amt) {
			amt = 1;
		}
	}
	
	boost::random::binomial_distribution<> dist( amt, _rate*g.length() );
	return dist(gen);
}

OpSequence* SequenceRecombinator::recombine(OpSequence& g1, OpSequence& g2) const {
	if (g1.key() == g2.key()) return &g1;
	
	int length = (g1.length() < g2.length()) ? g1.length() : g2.length();
	boost::random::binomial_distribution<> dist( length , _rate );
	boost::random::uniform_int_distribution<> idist( 0, length-10 );
	int num_sites = dist(gen);
	if (num_sites == 0) return &g1;

	//if (num_sites == 0) num_sites = 1;
	std::vector<int> locs;
	for (int i=0; i<num_sites; i++) {
		locs.push_back( idist(gen) );
	}
	
	sort(locs.begin(), locs.end());
	
	return new SequenceCrossover( g1, g2, locs );
	
}

double SequenceRecombinator::rate() const { return _rate; }

