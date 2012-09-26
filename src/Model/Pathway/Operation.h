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
				OpPathwayBase(double cost, int length, const GlobalInfo& info, OpPathway& parent1);
				OpPathwayBase(double cost, int length, const GlobalInfo& info, OpPathway& parent1, OpPathway& parent2);
				
				int numGenes() const;
				
				int numTFs() const;
				
				int numMotifs() const;
				
				int totalRegions() const;
				
				PTYPE get(int i) const;
				
				PTYPE getBinding(int i, int j) const;
				
			protected:
				virtual PTYPE proxyGet(int i) const = 0;
				
			private:
				int _length;
				const GlobalInfo& _info;
			};
			
			class PathwayRoot : public OpPathway {
			public:
				PathwayRoot( PromoterData* p );
				int length() const;
				PTYPE get(int i) const;
			};
			
			class PathwayRootFactory : public GenotypeFactory<PathwayRoot> {
			public:
				PathwayRootFactory( const GlobalInfo& info );
				
				PathwayRoot* random() const;
				
			private:
				const GlobalInfo& _info;
			}
			
			
		}
	}
}
#endif