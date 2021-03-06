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
#include "Util/Random.h"
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/discrete_distribution.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/random/binomial_distribution.hpp>
//#include <boost/math/distributions/binomial.hpp>
#include <algorithm>
#include <sstream>

// Seed for RNG
//boost::mt19937 gen;

using namespace GPPG::Model;
using namespace GPPG; 

/*******************************************************************
 *				OPSEQUENCEBASE OPERATION
 */

OpSequenceBase::OpSequenceBase( int cost, int length, OpSequence& parent1 ):
OpSequence(cost, parent1), _length(length) {}

OpSequenceBase::OpSequenceBase( int cost, int length, OpSequence& parent1, OpSequence& parent2 ):
OpSequence(cost, parent1, parent2), _length(length) {}

int OpSequenceBase::length() const { return _length; }

STYPE OpSequenceBase::get(int i) {
	incrRequests(1);
	if (isCompressed()) {
		return proxyGet(i);
	}
	return data()->get(i);
}

const char* OpSequenceBase::exportFormat() {
	std::stringstream output;
	const char* alpha = "ACTG";
	SequenceData* sd = evaluate();
	for(int i=0; i<length(); i++) {
		output << alpha[sd->get(i)];
	}
	
	return output.str().c_str();
}

SequenceRoot::SequenceRoot(SequenceData* d) : OperationRoot<SequenceData, ISequence>(d) { setCost(1); }

int SequenceRoot::length() const { return data()->length(); }
STYPE SequenceRoot::get(int i) { incrRequests(1); return data()->get(i); }


int discreteDistributionRandom(const std::vector<double>& distr) {
	double v = random01();
	for (int i=0;i<distr.size(); i++) 
		if (v <= distr[i] ) return i;
	return distr.size()-1;
}

SequenceData* randomSequenceData(int length, const std::vector<double>& distr) {
	SequenceData* sd = new SequenceData(length);
	
	for (int i=0; i<length; i++) {
		STYPE c = (STYPE)discreteDistributionRandom(distr);
		sd->set(i, c);
	}
	return sd;
}

void cumSum(std::vector<double>& d) {
	double csum=0;
	for (int i=0; i<d.size(); i++) {
		csum += d[i];
		d[i] = csum;
	}
}

SequenceRootFactory::SequenceRootFactory(int length, const std::vector<double>& distr ) : _length(length), _distr(distr) {
	cumSum(_distr);
}

SequenceRoot* SequenceRootFactory::random() const {
	
	return new SequenceRoot( randomSequenceData( _length, _distr) );
}



SequencePointChange::SequencePointChange(OpSequence& op, int* locs, int numLocs, STYPE* dest) : 
OpSequenceBase(numLocs,op.length(), op), _loc(locs), _numlocs(numLocs), _c(dest) {}

SequencePointChange::~SequencePointChange() { delete _loc; delete _c; }

SequenceData* SequencePointChange::evaluate()  {
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


STYPE SequencePointChange::proxyGet(int l)  {
	// See if the index is in the list
	for (int i=0; i<_numlocs; i++) {
		if (_loc[i] == l) return _c[i];
	}
	return parent(0)->get(l);
}

SequencePointMutator::SequencePointMutator(int cost, double rate, const std::vector<double> &T) : 
OperationMutator<OpSequence>(cost), _rate(rate), _M(T) {
	// Create random discrete distributions for each character
	int size = _M.size() >> 2;
	
	std::vector<double> weights( size );
	
	for (int i=0; i< size; i++) {
		for (int j=0; j< size; j++) {
			weights[j] = _M[i*size+j];
		}
		cumSum( weights );
		_transition.push_back( weights );
	}
}



OpSequence* SequencePointMutator::mutate( OpSequence& g) const {
	// Get the sequence data
	//int isCopy;
	//SequenceData* data = g.data(isCopy);
	
	// Calculate the number of sites to mutate (use binomial)
	int length = g.length(); //data->length()
	
	int numLocs = binomial( length, _rate ); 
	//if (numLocs == 0) return &g; 
	
	if (numLocs == 0) return &g;
	
#ifdef DEBUG_0
	std::cout << "SequencePointMutator: mutating..." << std::endl;
#endif
	
	int* locs = (int*) malloc(sizeof(int)*numLocs);
	STYPE* dest = (STYPE*) malloc(sizeof(STYPE)*numLocs);
	for (int i=0; i<numLocs; i++) {
		int loc = (int)(random01()*length);
		STYPE c = g.get(loc); //data->get(loc);
		locs[i] = loc;
		dest[i] = (STYPE)( discreteDistributionRandom(_transition[c]) );
		
	}
	SequencePointChange* spc = new SequencePointChange(g, locs, numLocs, dest);
	
#ifdef DEBUG_0
	std::cout << "Created Operation (" << spc->numParents() << "): " << spc->toString() << std::endl;
#endif
	//if (isCopy) delete data;
	spc->setCost( cost() );
	return spc;
}

double SequencePointMutator::rate() const { return _rate; }

const std::vector<double>& SequencePointMutator::transition() const { return _M; }

int SequencePointMutator::numMutants(OpSequence& g, long N, double f) const {
	return binomial(N*f, _rate*g.length() );
}

ostream& operator<<(ostream& output, const SequencePointMutator& s) {
	output << "SequencePointMutator(rate=" << s.rate() << "): ?"; // << s.transition();
	return output;
}


SequenceDeletion::SequenceDeletion(OpSequence& op, int loc, int span): OpSequenceBase(10, op.length()-span, op), _loc(loc), _span(span) {}

SequenceData* SequenceDeletion::evaluate() {
	SequenceData* sd = OpSequenceBase::evaluate();
	if (sd != NULL) return sd;
	
	// Get the sequence from the parent and add the point changes
	sd = parent(0)->evaluate();
	
	
	SequenceData* data = new SequenceData( length() );
	
	STYPE* sdata = data->sequence();
	STYPE* pdata = sd->sequence();
	if (_loc > 0) 
		memcpy(sdata, pdata, sizeof(STYPE)*_loc);
	memcpy(&sdata[_loc], &pdata[_loc+_span] , sizeof(STYPE)*(length()-_loc));

	delete sd;
	
	return data;
}

STYPE SequenceDeletion::proxyGet(int i)  {
	// See if the index is in the list
	if (i >= _loc) {
		return parent(0)->get(i+_span);
	}
	return parent(0)->get(i);
}

std::string SequenceDeletion::toString() const {
	std::ostringstream output;
	output << _span << " nt @ " << _loc;
	
	return output.str();
}

SequenceInsertion::SequenceInsertion(OpSequence& op, int loc, SequenceData* span): 
	OpSequenceBase(10, op.length()+span->length(), op), _loc(loc), _span(span) {}

SequenceInsertion::~SequenceInsertion() { delete _span; }

SequenceData* SequenceInsertion::evaluate() {
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

STYPE SequenceInsertion::proxyGet(int i)  {
	// See if the index is in the list
	int index = i;
	if (i >= _loc+_span->length()) {
		index = i-_span->length();
	} else if (i >= _loc && i < _loc+_span->length()) {
		return _span->get( i-_loc );
	}
	
	return parent(0)->get(index);
}

std::string SequenceInsertion::toString() const {
	std::ostringstream output;
	output << _span->length() << " nt @ " << _loc;
	
	return output.str();
}

SequenceDeletionMutator::SequenceDeletionMutator(int cost, double rate, int minL, int maxL) :
OperationMutator<OpSequence>(cost), _rate(rate), _minL(minL), _maxL(maxL) {}

OpSequence* SequenceDeletionMutator::mutate( OpSequence& g) const {
	if (binomial(g.length(), _rate) == 0) return &g;
	
	// Perform deletion	
	int length = g.length(); 
	
	int spanLength = (int)(random01()*(_maxL-_minL))+_minL;
	
	int loc = (int)(random01()*(length-spanLength));
	
	SequenceDeletion *sd = new SequenceDeletion(g, loc, spanLength);
	sd->setCost(cost());
	return sd;
	
}

int SequenceDeletionMutator::numMutants(OpSequence& g, long N, double f) const {
	return binomial(N*f, _rate*g.length());
}

double SequenceDeletionMutator::rate() const { return _rate; }


SequenceInsertionMutator::SequenceInsertionMutator(int cost, double rate, int minL, int maxL, const std::vector<double>& distr) :
OperationMutator<OpSequence>(cost), _rate(rate), _minL(minL), _maxL(maxL), _distr(distr) {
	cumSum(_distr);
}

OpSequence* SequenceInsertionMutator::mutate( OpSequence& g) const {
	// See if insertion occurs
	if (binomial(g.length(), _rate) == 0) return &g;
	
	// Perform insertion
	int length = g.length(); 
	
	int spanLength = (int)(random01()*(_maxL-_minL))+_minL;
	
	int loc = (int)(random01()*length);
	
	// Generate random sequence
	SequenceData* span = randomSequenceData(spanLength, _distr);
	
	SequenceInsertion *sd = new SequenceInsertion(g, loc, span);
	sd->setCost( cost() );
	return sd;
}

int SequenceInsertionMutator::numMutants(OpSequence& g, long N, double f) const {
	return binomial(N*f, _rate*g.length());
}

double SequenceInsertionMutator::rate() const { return _rate; }


/*******************************************************************
 * SEQUENCE RECOMBINATION
 *******************************************************************/

SequenceCrossover::SequenceCrossover(OpSequence& op1, OpSequence& op2, const std::vector<int>& locs) :
OpSequenceBase(100, (locs.size()%2==0)? op1.length(): op2.length(), op1, op2), _locs(locs) {}


SequenceData* SequenceCrossover::evaluate() {
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

STYPE SequenceCrossover::proxyGet(int i)  {
	// See if the index is in the list
	parent(0)->touch();
	parent(1)->touch();	
	for (int j=0; j<_locs.size(); j++) {
		if (i<j) {
			return (j%2==0) ? parent(0)->get(i) : parent(1)->get(i);
		}
	}
	return (_locs.size()%2==0) ? parent(0)->get(i) : parent(1)->get(i);
}

std::string SequenceCrossover::toString() const {
	std::ostringstream output;
	output << "points=";
	for (int j=0; j<_locs.size(); j++) {
		if(j>0) output << ",";
		output << _locs[j];
	}
	
	return output.str();
}

SequenceRecombinator::SequenceRecombinator(int cost, double rate) : OperationRecombinator<OpSequence>(cost), _rate(rate) {}

int SequenceRecombinator::numMutants(OpSequence& g, OpSequence& g2, long N) const {
	double amt = N*g.frequency()*g2.frequency();
	if (amt <= 0) return 0;
	
	else if (amt < 1) {
		if (random01() < amt) {
			amt = 1;
		}
	}
	
	return binomial(amt, _rate*g.length());
}

OpSequence* SequenceRecombinator::recombine(OpSequence& g1, OpSequence& g2) const {
	if (g1.key() == g2.key()) return &g1;
	
	int length = (g1.length() < g2.length()) ? g1.length() : g2.length();

	int num_sites = binomial(length, _rate);
	if (num_sites == 0) return &g1;

	//if (num_sites == 0) num_sites = 1;
	std::vector<int> locs;
	for (int i=0; i<num_sites; i++) {
		locs.push_back( (int)(random01()*(length-10)) );
	}
	
	sort(locs.begin(), locs.end());
	
	SequenceCrossover* sc = new SequenceCrossover( g1, g2, locs );
	sc->setCost( cost() );
	return sc;
	
}

double SequenceRecombinator::rate() const { return _rate; }

