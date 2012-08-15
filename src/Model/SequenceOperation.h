/*
 *  SequenceOperation.h
 *  GPPG
 *
 *  Created by Troy Ruths on 8/13/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef SEQUENCE_OPERATION_
#define SEQUENCE_OPERATION_

#include "Operation.h"
#include "Mutator.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/array.hpp>
#include <iostream>

// Use only the end of the boost namespace to avoid collisions
using namespace boost::numeric;
using std::ostream;

namespace GPPG {

namespace Model {
	typedef short STYPE;	
	
	/**
	 * A simple data structure for holding sequence information.
	 */
	class SequenceData {
	public:
		/** Allocates \param length sequence.
		 *
		 */
		SequenceData(int length);
		
		/** Deep-copy
		 */
		SequenceData* copy() const;
		
		STYPE* sequence();
		
		int length() const;
		
		/** Gets the item at location i
		 */
		inline STYPE get(int i) const;
		
		/** Sets the item at location i
		 */
		inline void set(int i, STYPE c);
		
	private:
		STYPE* _sequence;
		int _length;
	};
	
	
	
	class SequenceFactory : public OperationFactory<SequenceData> {
	public:
		/** Creates a SequenceFactory that returns sequences of length \param length and character distribution \param distr
		 * This factory returns SequenceRoot objects.
		 */
		SequenceFactory(int length, const ublas::vector<double>& distr );
		
		SequenceData* randomData() const;
		
	private:
		int _length;
		ublas::vector<double> _distr;
	};
	
	/** Provides a root (non-op) storage for SequenceData.
	 * We just need to typedef it because the template does all the work.
	 */
	typedef OperationRoot<SequenceData> SequenceRoot;
	
	class SequencePointChange: public Operation<SequenceData> {
	public:
		SequencePointChange(Operation<SequenceData>& op, int* locs, int numLocs, STYPE* dest);
		
		~SequencePointChange(); 
		
		SequenceData* evaluate() const;
		
	private:
		int* _loc;		/* Locations array */
		int _numlocs;	/* Number of locations to change */
		STYPE* _c;		/* Characters to be changed to */
	};
	
	class SequencePointMutator : public OperationMutator<SequenceData, SequencePointChange> {
	public:
		/** Generates point mutations with \param rate and transition matrix \T.
		 */
		SequencePointMutator(double rate, const ublas::matrix<double> &T);
		
		Operation<SequenceData>* mutate( const Operation<SequenceData>& g ) const; 
		
		double rate() const;
		const ublas::matrix<double>& transition() const;
		
	private:
		double _rate; 
		ublas::matrix<double> _M; /* Transition matrix */
	};
	
	class SequenceDeletion: public Operation<SequenceData> {
		
		SequenceDeletion(Operation<SequenceData>& op, int loc, int span);
		
		SequenceData* evaluate() const;
		
	private:
		int _loc, _span;
	};
	
	class SequenceInsertion: public Operation<SequenceData> {
	public:
		SequenceInsertion(Operation<SequenceData>& op, int loc, SequenceData* span);
		
		SequenceData* evaluate() const;
		
	private:
		int _loc;
		SequenceData* _span;
	};
	
}
}

ostream& operator<<(ostream& output, const GPPG::Model::SequenceData& s);	
ostream& operator<<(ostream& output, const GPPG::Model::SequencePointChange& s);	
ostream& operator<<(ostream& output, const GPPG::Model::SequencePointMutator& s);	
#endif

//class SequenceOperation : public Operation<SequenceData> {
//public:
//	SequenceOperation( double cost, int length );
//	SequenceOperation( double cost, int length, SequenceOperation& parent1 );
//	SequenceOperation( double cost, int length, SequenceOperation& parent1, SequenceOperation& parent2 );
//	
//	
//	/** Gets the character at position i.
//	 * @param i - location to retrieve
//	 */
//	virtual STYPE get(int i) const;
//	
//	/** Returns the length of the sequence.
//	 *
//	 */
//	int length() const;
//	
//	
//private:
//	int _length;
//};


