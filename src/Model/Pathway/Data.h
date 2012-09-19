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
			
			/** Gets the item at location i
			 */
			virtual PTYPE get(int i) const = 0;
		};
		
		
		/**
		 * A simple data structure for holding sequence information.
		 */
		class TransRegPathwayData : ITransRegPathway {
		public:
			
			TransRegPathwayData();
			
			/** Deep-copy
			 */
			TransRegPathwayData* copy() const;
			
			/** Get a pointer to the raw sequence data.
			 */
			PTYPE* sequence();
			
			int length() const;
			
			/** Gets the item at location i
			 */
			PTYPE get(int i) const;
			
			/** Sets the item at location i
			 */
			void set(int i, PTYPE c);
			
		private:
			int *_offset, *_tfs, *_regions, *_numSites;
			PTYPE *_pool;
			int _numGenes, _numTFs, _numMotifs, _totalRegions;
		};
	}
}

#endif