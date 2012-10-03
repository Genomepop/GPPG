/*
 *  Operation.h
 *  Demo
 *
 *  Created by Troy Ruths on 9/24/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef PATHWAY_OPERATION_
#define PATHWAY_OPERATION_

#include <Operation/Operation.h>
#include <Model/Pathway/Data.h>

namespace GPPG {
	namespace Model {
		namespace TransReg {
			typedef Operation<PromoterData, ITransRegPathway> OpPathway;
			
			class OpPathwayBase : public OpPathway {
			public:
				OpPathwayBase(double cost, const GlobalInfo& info, OpPathway& parent1);
				OpPathwayBase(double cost, const GlobalInfo& info, OpPathway& parent1, OpPathway& parent2);
				
				int numGenes() const;
				
				int numTFs() const;
				
				int numMotifs() const;
				
				int totalRegions() const;
				
				PTYPE get(int i) const;
				
				PTYPE getBinding(int i, int j) const;
				
				const GlobalInfo& info() const;
				
			protected:
				virtual PTYPE proxyGet(int i) const = 0;
				
				const GlobalInfo& _info;
			};
			
			class PathwayRoot : public OperationRoot< PromoterData, ITransRegPathway > {
			public:
				PathwayRoot( PromoterData* p );
				int numGenes() const;
				
				int numTFs() const;
				
				int numMotifs() const;
				
				int totalRegions() const;
				
				PTYPE get(int i) const;
				
				PTYPE getBinding(int i, int j) const;
				
				const GlobalInfo& info() const;
			};
			
			class PathwayRootFactory : public GenotypeFactory<PathwayRoot> {
			public:
				PathwayRootFactory( const GlobalInfo& info );
				
				PathwayRoot* random() const;
				
				static GlobalInfo* randomInfo(int numGenes, int numTFs, int minRegions, int maxRegions);
				
			private:
				const GlobalInfo& _info;
			};
			
			class BindingSiteChange : public OpPathwayBase {
			public:
				BindingSiteChange(OpPathway& op, std::vector<int>* locs, std::vector<PTYPE>* dest);
				
				~BindingSiteChange(); 
				
				std::string toString() const;
				
				PromoterData* evaluate() const;
				
				/** Get the number of mutated sites
				 */
				int numSites() const;
				
				/** Retrieve the \param i'th mutation.
				 */
				PTYPE getMutation(int i) const;
				
				/** Retrieve the \param i'th site.
				 */
				int getSite(int i) const;
				
				
			protected:
				PTYPE proxyGet(int i) const;
				
			private:
				std::vector<int>* _locs;		/* Locations array */
				std::vector<PTYPE>* _c;		/* Characters to be changed to */
				
			};
			
			class PromoterCrossover : public OpPathwayBase {
				
			};
			
			/** The BindingSiteMutator mutates binding sites in a promoter.
			 * Binding sites correspond to motifs.  Motifs have different spontaneous loss and gain rates.
			 * This mutator simulates both motif loss and gain for a set of motifs.
			 */
			class BindingSiteMutator : public OperationMutator< OpPathway > {
			public:
				BindingSiteMutator( double u, int motifOverlap, const std::vector<double>& motifGainRates, const std::vector<double>& motifLossProb);
				~BindingSiteMutator();
				
				OpPathway* mutate( OpPathway& g ) const;
				
				int numMutants(OpPathway& g, long N, double f) const;
				
				double rate() const;
				
			private:
				
				double _u;
				int _overlap;
				std::vector<double> _gainRates, _lossProb;
			};
			
			class PromoterRecombinator : public OperationRecombinator< OpPathway > {
				
			};
			
		}
	}
}
#endif