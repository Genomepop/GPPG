/*
 *  Data.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 9/19/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "Data.h"

using namespace GPPG::Model::TransReg;

GlobalInfo::GlobalInfo(const std::vector<std::string>& genes,
		   const std::vector<int>& regions,
		   const std::vector<std::string>& motifs,
		   const std::vector<int>& tfs,
			const std::map< int, std::vector<int> >& binding) :
_genes(genes), _regions(regions), _motifs(motifs), _tfs(tfs), _binding(binding) {
	
	// Calculate offset
	int offset = 0;
	for (int i=0; i<_regions.size(); i++) {
		_offset.push_back( offset );
		offset += _regions[i];
	}
}

int GlobalInfo::numTFs() const { return _tfs.size(); }

int GlobalInfo::numGenes() const { return _genes.size(); }

int GlobalInfo::numMotifs() const { return _motifs.size(); }

const std::string& GlobalInfo::getGeneName(int i) const { return _genes[i]; }


const std::string& GlobalInfo::getMotifPWM(int i) const { return _motifs[i]; }

int GlobalInfo::getTF(int i) const { return _tfs[i]; }

int GlobalInfo::getGeneForRegion(int i) const {
	return 0;
}

PromoterData::PromoterData(const GlobalInfo& info, int totalRegions): _info(info), _totalRegions(totalRegions) {
	_pool = (PTYPE*)malloc(sizeof(PTYPE)*_totalRegions);
}

PromoterData::~PromoterData() {
	delete _pool;
}

PromoterData* PromoterData::copy() const {
	PromoterData* pd = new PromoterData(_info, _totalRegions);
	memcpy( pd->_pool, _pool, sizeof(PTYPE)*_totalRegions);
	return pd;
}

PTYPE PromoterData::get(int i) const { return _pool[i]; }

int PromoterData::totalRegions() const { return _totalRegions; }



