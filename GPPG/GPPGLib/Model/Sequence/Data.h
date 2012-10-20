/*
 *  SequenceData.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/17/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef SEQUENCE_DATA_
#define SEQUENCE_DATA_

namespace GPPG {
	
	namespace Model {
		typedef short STYPE;	
		
		class ISequence {
		public:
			/** Retrieve the length of the sequence.
			 */
			virtual int length() const = 0;
			
			/** Gets the item at location i
			 */
			virtual STYPE get(int i)  = 0;
		};
		
		/**
		 * A simple data structure for holding sequence information.
		 */
		class SequenceData : ISequence {
		public:
			/** Allocates \param length sequence.
			 *
			 */
			SequenceData(int length);
			
			~SequenceData();
			
			/** Deep-copy
			 */
			SequenceData* copy() const;
			
			/** Get a pointer to the raw sequence data.
			 */
			STYPE* sequence();
			
			int length() const;
			
			/** Gets the item at location i
			 */
			STYPE get(int i) ;
			
			STYPE get(int i) const;
			
			/** Sets the item at location i
			 */
			void set(int i, STYPE c);
			
		private:
			STYPE* _sequence;
			int _length;
		};
	}
}

#endif
