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
		namespace TransReg {
			
		typedef unsigned short PTYPE;	
		
			class GlobalInfo  {
			public:
				GlobalInfo(const std::vector<std::string>& genes,
						   const std::vector<int>& regions,
						   const std::vector<std::string>& motifs,
						   const std::vector<int>& tfs,
						const std::map< int, std::vector<int> >& binding); 
						
				GlobalInfo(const std::vector<std::string>& genes,
						   const std::vector<int>& regions,
						   const std::vector<std::string>& motifs,
						   const std::vector<int>& tfs,
						   const std::map< int, std::vector<int> >& binding,
							const std::map<std::string, std::string>& motifSeq
							);
				
				int numTFs() const;
				
				int numGenes() const;
				
				int numMotifs() const;
				
				const std::string& getGeneName(int i) const;
				
				const std::string& getMotifName(int i) const;				
				
				const std::string& getMotifSequence(const std::string& motifName) const;
				
				int getTF(int i) const;
				
				int getGeneForRegion(int i) const;
				
				int totalRegions() const;
				
				int numRegions(int i) const;
				
				int offset(int i) const;
				
				std::vector<int> binding(int i) const;
				
				/** Retrieve network structure from promoter data
				 */
				
			private:
				std::vector<std::string> _genes, _motifs;
				std::vector<int> _regions, _tfs, _offset;
				// Maps MotifID -> TFIDs, used in creating interaction network
				std::map< int, std::vector<int> > _binding;
				std::map< std::string, std::string > _motifSeq;
				int _totalRegions;
			};
			
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
			
			/** Retrieve number of regions which may contain binding sites
			 */
			virtual int totalRegions() const = 0;
			
			/** Gets the item at location i
			 */
			virtual PTYPE get(int i)  = 0;
			
			/** Get the binding motif for promoter \param i , location \param j.
			 */
			virtual PTYPE getBinding(int i, int j)  = 0;
			
			/** Gets the number of binding sites in a region
			 */
			virtual short numSitesForGene(int i) = 0;
			
			virtual const GlobalInfo& info() const = 0;
			
		};
		

			
			/**
			 * A simple data structure for holding Pathway (promoter) information.
			 */
			class PromoterData : ITransRegPathway {
			public:
				
				PromoterData(const GlobalInfo& info);
				
				~PromoterData();
				
				/** Deep-copy
				 */
				PromoterData* copy() const;
				
				/** Clears all binding site information
				*/
				void clearData();
				
				/** Get a pointer to the raw sequence data.
				 */
				PTYPE* data();
				
				/** Gets the number of regions
				 */
				int totalRegions() const;
				
				/** Gets the BS located at region i
				 */
				PTYPE get(int i) ;
				
				/** Sets the item at location i
				 */
				void set(int i, PTYPE c);
				
				void set(int gene, int region, PTYPE c);
				
				int numGenes() const;
				
				int numTFs() const;
				
				int numMotifs() const;
			
				PTYPE getBinding(int i, int j) ;
				
				short numSitesForGene(int i);
				
				const GlobalInfo& info() const;
				
			private:
				PTYPE *_pool;
				short *_numSites;
				const GlobalInfo& _info;
			};
		}
	}
}

#endif