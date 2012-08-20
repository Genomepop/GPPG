/*
 *  Operation.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/17/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef SEQUENCE_OPERATION_
#define SEQUENCE_OPERATION_

#include <Operation/Operation.h>
#include <Base/Mutator.h>
#include <Model/Sequence/Data.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/array.hpp>

#include <iostream>

// Use only the end of the boost namespace to avoid collisions
using namespace boost::numeric;
using std::ostream;

namespace GPPG {
	
	namespace Model {
		
		typedef Operation<SequenceData> SequenceOperation;	
		
		/** Provides a root (non-op) storage for SequenceData.
		 * We just need to typedef it because the template does all the work.
		 */
		typedef OperationRoot<SequenceData> SequenceRoot;
		
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
		

		
		
		class SequencePointChange: public Operation<SequenceData> {
		public:
			SequencePointChange(Operation<SequenceData>& op, int* locs, int numLocs, STYPE* dest);
			
			~SequencePointChange(); 
			
			std::string toString() const;
			
			SequenceData* evaluate() const;
			
			/** Get the number of mutated sites
			 */
			int numSites() const;
			
			/** Retrieve the \param i'th mutation.
			 */
			STYPE getMutation(int i) const;
			
			/** Retrieve the \param i'th site.
			 */
			int getSite(int i) const;
			
		private:
			int* _loc;		/* Locations array */
			int _numlocs;	/* Number of locations to change */
			STYPE* _c;		/* Characters to be changed to */
		};
		
		class SequencePointMutator : OperationMutator<SequenceData> {
		public:
			/** Generates point mutations with \param rate and transition matrix \T.
			 */
			SequencePointMutator(double rate, const ublas::matrix<double> &T);
			
			Operation<SequenceData>* mutate( Operation<SequenceData>& g ) const; 
			
			int numMutants(IGenotype& g, long N, double f) const;
			
			double rate() const;
			const ublas::matrix<double>& transition() const;
			
		private:
			double _rate; 
			ublas::matrix<double> _M; /* Transition matrix */
			std::vector<boost::random::discrete_distribution<> > _transition;
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