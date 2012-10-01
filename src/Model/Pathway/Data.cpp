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
using namespace GPPG::Model;

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
	_totalRegions = offset;
}

int GlobalInfo::numTFs() const { return _tfs.size(); }

int GlobalInfo::numGenes() const { return _genes.size(); }

int GlobalInfo::numMotifs() const { return _motifs.size(); }

int GlobalInfo::totalRegions() const { return _totalRegions; }

const std::string& GlobalInfo::getGeneName(int i) const { return _genes[i]; }


const std::string& GlobalInfo::getMotifPWM(int i) const { return _motifs[i]; }

int GlobalInfo::getTF(int i) const { return _tfs[i]; }

int GlobalInfo::getGeneForRegion(int i) const {
	throw "Not Implemented";
	return 0;
}

PromoterData::PromoterData(const GlobalInfo& info): _info(info) {
	_pool = (PTYPE*)malloc(sizeof(PTYPE)*_info.totalRegions());
}

PromoterData::~PromoterData() {
	delete _pool;
}

PromoterData* PromoterData::copy() const {
	PromoterData* pd = new PromoterData(_info);
	memcpy( pd->_pool, _pool, sizeof(PTYPE)*_info.totalRegions());
	return pd;
}

void PromoterData::set(int i, PTYPE c) {
	_pool[i] = c;
}

PTYPE PromoterData::get(int i) const { return _pool[i]; }

int PromoterData::totalRegions() const { return _info.totalRegions(); }

int PromoterData::numGenes() const { return _info.numGenes(); }

int PromoterData::numTFs() const { return _info.numTFs(); }

int PromoterData::numMotifs() const { return _info.numMotifs(); }

PTYPE PromoterData::getBinding(int i, int j) const { 
	throw "Not Implemented";
	return 0; 
}

const GlobalInfo& PromoterData::info() const { return _info; }


