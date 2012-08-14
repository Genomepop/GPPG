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
	};
	
	class SequenceOperation : public Operation<SequenceData> {
	public:
		SequenceOperation( double cost, int length );
		SequenceOperation( double cost, int length, SequenceOperation& parent1 );
		SequenceOperation( double cost, int length, SequenceOperation& parent1, SequenceOperation& parent2 );
		
		
		/** Gets the character at position i.
		 * @param i - location to retrieve
		 */
		virtual STYPE get(int i) const;
		
		/** Returns the length of the sequence.
		 *
		 */
		int length() const;

		
	private:
		int _length;
	};
	
	class SequenceOperationRoot: public SequenceOperation {
	public:
		SequenceOperationRoot( SequenceData* data );
	};
	
	class SequencePointChange: public SequenceOperation {
	public:
		SequencePointChange(SequenceOperation& op, int* locs, int numLocs, STYPE* dest);
		
		~SequencePointChange();
		
		STYPE get(int i) const;
		
	protected:
		SequenceData* evaluate() const;
		
	private:
		int* _loc;		/* Locations array */
		int _numlocs;	/* Number of locations to change */
		STYPE* _c;		/* Characters to be changed to */
	};
	
	class SequenceDeletion: public SequenceOperation {
		
		SequenceDeletion(SequenceOperation& op, int loc, int span);
		
		STYPE get(int i) const;

	protected:
		SequenceData* evaluate() const;
		
	private:
		int _loc, _span;
	};
	
	class SequenceInsertion: public SequenceOperation {
	public:
		SequenceInsertion(SequenceOperation& op, int loc, SequenceData* span);
		
	protected:
		SequenceData* evaluate() const;
		
	private:
		int _loc;
		SequenceData* _span;
	};
	
}
}

#endif

