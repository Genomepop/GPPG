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
				
		typedef Operation<SequenceData, ISequence> OpSequence;		
		
		
		class OpSequenceBase : public OpSequence {
		public:
			OpSequenceBase(double cost, int length, OpSequence& parent1);
			OpSequenceBase(double cost, int length, OpSequence& parent1, OpSequence& parent2);
			
			int length() const;
			
			STYPE get(int i) const;
			
		protected:
			virtual STYPE proxyGet(int i) const = 0;
			
		private:
			int _length;
		};
		

		/** Provides a root (non-op) storage for SequenceData.
		 * We just need to typedef it because the template does all the work.
		 */
		class SequenceRoot : public OperationRoot<SequenceData,ISequence> {
		public:
			SequenceRoot( SequenceData* d);
			int length() const;
			STYPE get(int i) const;
			
		}; 
		//OperationRoot<SequenceData, ISequenceData> SequenceRoot;
		
		class SequenceRootFactory : public GenotypeFactory<SequenceRoot> {
		public:
			/** Creates a SequenceFactory that returns sequences of length \param length and character distribution \param distr
			 * This factory returns SequenceRoot objects.
			 */
			SequenceRootFactory(int length, const ublas::vector<double>& distr );
			
			SequenceRoot* random() const;
			
		private:
			int _length;
			ublas::vector<double> _distr;
		};
		

		
		
		class SequencePointChange: public OpSequenceBase {
		public:
			SequencePointChange(OpSequence& op, int* locs, int numLocs, STYPE* dest);
			
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
			
		protected:
			STYPE proxyGet(int i) const;
			
		private:
			int* _loc;		/* Locations array */
			int _numlocs;	/* Number of locations to change */
			STYPE* _c;		/* Characters to be changed to */
			int _length;
		};
		
		class SequencePointMutator : public OperationMutator< OpSequence > {
		public:
			/** Generates point mutations with \param rate and transition matrix \T.
			 */
			SequencePointMutator(double rate, const ublas::matrix<double> &T);
			
			OpSequence* mutate( OpSequence& g ) const; 
			
			int numMutants(OpSequence& g, long N, double f) const;
			
			double rate() const;
			const ublas::matrix<double>& transition() const;
			
		private:
			double _rate; 
			ublas::matrix<double> _M; /* Transition matrix */
			std::vector<boost::random::discrete_distribution<> > _transition;
		};
		
		
		
		class SequenceDeletion: public OpSequenceBase {
		public:
			SequenceDeletion(OpSequence& op, int loc, int span);
			
			SequenceData* evaluate() const;
			
		protected:
			STYPE proxyGet(int i) const;
			
		private:
			int _loc, _span;
		};
		
		class SequenceDeletionMutator : public OperationMutator< OpSequence > {
		public:
			/** Generates sequence deletions mutations with \param rate and length between \param minL and \param maxL.
			 */
			SequenceDeletionMutator(double rate, int minL, int maxL);
			
			OpSequence* mutate( OpSequence& g ) const; 
			
			int numMutants(OpSequence& g, long N, double f) const;
			
			double rate() const;
			
			
		private:
			int _minL, _maxL;
			double _rate;
		};
		
		
		class SequenceInsertion: public OpSequenceBase {
		public:
			SequenceInsertion(OpSequence& op, int loc, SequenceData* span);
			~SequenceInsertion();
			
			SequenceData* evaluate() const;
			
		protected:
			STYPE proxyGet(int i) const;
			
		private:
			int _loc;
			SequenceData* _span;
		};
		
		class SequenceInsertionMutator : public OperationMutator< OpSequence > {
		public:
			/** Generates sequence deletions mutations with \param rate and length between \param minL and \param maxL.
			 */
			SequenceInsertionMutator(double rate, int minL, int maxL, const ublas::vector<double>& distr);
			
			OpSequence* mutate( OpSequence& g ) const; 
			
			int numMutants(OpSequence& g, long N, double f) const;
			
			double rate() const;
			
			
		private:
			int _minL, _maxL;
			double _rate;
			ublas::vector<double> _distr;
		};
		
		class SequenceCrossover: public OpSequenceBase {
		public:
			SequenceCrossover(OpSequence& op1, OpSequence& op2, const std::vector<int>& locs);
			SequenceCrossover();
			
			SequenceData* evaluate() const;
			
		protected:
			STYPE proxyGet(int i) const;
			
		private:
			std::vector<int> _locs;
		};
		
		class SequenceRecombinator : public OperationRecombinator< OpSequence > {
		public:
			SequenceRecombinator(double rate);
			
			int numMutants(OpSequence& g, OpSequence& g2, long N) const;
			
			OpSequence* recombine(OpSequence& g1, OpSequence& g2) const;
			
			double rate() const;
			
		private:
			double _rate;
		};
	}
}



#endif