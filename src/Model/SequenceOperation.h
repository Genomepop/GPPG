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

namespace GPPG {

namespace Model {
	typedef short STYPE;	
	/**
	 * A simple data structure for holding sequence information.
	 */
	class SequenceData {
	public:
		STYPE* sequence;
		int length;
		
		SequenceData* copy() const;
		
	};
	
	
	class SequenceOperationRoot: public Operation<SequenceData> {
	public:
		SequenceOperationRoot( SequenceData* data );
	};
	
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


