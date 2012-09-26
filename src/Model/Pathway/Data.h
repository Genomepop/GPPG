/*
 *  Data.h
 *  Demo
 *
 *  Created by Troy Ruths on 9/19/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef PATHWAY_DATA_
#define PATHWAY_DATA_

#include <vector>
#include <string>
#include <map>

namespace GPPG {
	
	namespace Model {
		typedef unsigned short PTYPE;	
		
		
		
		class ITransRegPathway {
		public:
			/** Retrieve the number of genes
			 */
			virtual int numGenes() const = 0;
			
			/** Retrieve the number of TFs
			 */
			virtual int numTFs() const = 0;
			
			/** Retrieve the number of motifs
			 */
			virtual int numMotifs() const = 0;
			
			virtual int totalRegions() const = 0;
			
			/** Gets the item at location i
			 */
			virtual PTYPE get(int i) const = 0;
			
			virtual PTYPE getBinding(int i, int j) const = 0;
			
		};
		
		namespace TransReg {

			class GlobalInfo  {
			public:
				GlobalInfo(const std::vector<std::string>& genes,
						   const std::vector<int>& regions,
						   const std::vector<std::string>& motifs,
						   const std::vector<int>& tfs,
						   const std::map< int, std::vector<int> >& binding);
				
				int numTFs() const;
				
				int numGenes() const;
				
				int numMotifs() const;
				
				const std::string& getGeneName(int i) const;
				
				const std::string& getMotifPWM(int i) const;				
				
				int getTF(int i) const;
				
				int getGeneForRegion(int i) const;
				
				/** Retrieve network structure from promoter data
				 */
				
			private:
				std::vector<std::string> _genes, _motifs;
				std::vector<int> _regions, _tfs, _offset;
				std::map< int, std::vector<int> > _binding;
			};
			
			/**
			 * A simple data structure for holding Pathway (promoter) information.
			 */
			class PromoterData : ITransRegPathway {
			public:
				
				PromoterData(const GlobalInfo& info, int totalRegions);
				
				~PromoterData();
				
				/** Deep-copy
				 */
				PromoterData* copy() const;
				
				/** Get a pointer to the raw sequence data.
				 */
				PTYPE* sequence();
				
				/** Gets the number of regions
				 */
				int totalRegions() const;
				
				/** Gets the BS located at region i
				 */
				PTYPE get(int i) const;
				
				/** Sets the item at location i
				 */
				void set(int i, PTYPE c);
				
			private:
				PTYPE *_pool;
				int _totalRegions;
				const GlobalInfo& _info;
			};
		}
	}
}

#endif