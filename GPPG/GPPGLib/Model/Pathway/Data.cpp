/*
 *  Data.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 9/19/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "Data.h"
#include <algorithm>
#include <iostream>

using namespace GPPG::Model::TransReg;
using namespace GPPG::Model;
using namespace std;

GlobalInfo::GlobalInfo(const std::vector<std::string>& genes,
		   const std::vector<int>& regions,
		   const std::vector<std::string>& motifs,
		   const std::vector<int>& tfs,
		const std::map< int, std::vector<int> >& binding): 
		_genes(genes), _regions(regions), _motifs(motifs), _tfs(tfs), _binding(binding)
{
	initialize();
}

GlobalInfo::GlobalInfo(const std::vector<std::string>& genes,
		   const std::vector<int>& regions,
		   const std::vector<std::string>& motifs,
		   const std::vector<int>& tfs,
			const std::map< int, std::vector<int> >& binding,
			const std::map<std::string, std::string>& motifSeq) :
_genes(genes), _regions(regions), _motifs(motifs), _tfs(tfs), _binding(binding), _motifSeq(motifSeq) {
	
	initialize();
}

void GlobalInfo::initialize() {
	// Calculate offset
	int offset = 0;
	for (int i=0; i<_regions.size(); i++) {
		_offset.push_back( offset );
		offset += _regions[i];
	}
	_totalRegions = offset;
	
	// Create TF->Motif
	for( map<int, vector<int> >::iterator it=_binding.begin(); it!=_binding.end(); it++ ) {
		vector<int>& v = it->second;
		int motifID = it->first;
		for( vector<int>::iterator vit=v.begin(); vit!=v.end();vit++) {
			int tfID = *vit;
			_bindingTF[tfID].push_back(motifID);
		}
	}
}

int GlobalInfo::numTFs() const { return _tfs.size(); }

int GlobalInfo::numGenes() const { return _genes.size(); }

int GlobalInfo::numMotifs() const { return _motifs.size(); }

int GlobalInfo::totalRegions() const { return _totalRegions; }

int GlobalInfo::offset(int i) const { return _offset[i]; }

int GlobalInfo::numRegions(int i) const { return _regions[i]; }

const std::string& GlobalInfo::getGeneName(int i) const { return _genes[i]; }


const std::string& GlobalInfo::getMotifName(int i) const { return _motifs[i]; }

const std::string& GlobalInfo::getMotifSequence(const std::string& motifName) const { return _motifSeq.at(motifName); }

int GlobalInfo::getTF(int i) const { return _tfs[i]; }

int GlobalInfo::getGeneForRegion(int i) const {
	return lower_bound( _offset.begin(), _offset.end(), i )-_offset.begin();
}

const std::vector<int>& GlobalInfo::binding(int i) const {
	return _binding.find(i)->second;
}


const std::vector<int>& GlobalInfo::bindingTFsForTF(int motifID) const {
	return binding(motifID);
}

const std::vector<int>& GlobalInfo::bindingMotifsForTF(int tfID) const {
	return _bindingTF.find(tfID)->second;
}

PromoterData::PromoterData(const GlobalInfo& info): _info(info), _pool(0), _numSites(0) {
	_pool = (PTYPE*)malloc(sizeof(PTYPE)*_info.totalRegions());
	_numSites = (STYPE*)malloc(sizeof(STYPE)*_info.numGenes());
}

PromoterData::~PromoterData() {
	if(_pool) delete _pool;
	if(_numSites) delete _numSites;
}

PromoterData* PromoterData::copy() const {
	PromoterData* pd = new PromoterData(_info);
	memcpy( pd->_pool, _pool, sizeof(PTYPE)*_info.totalRegions());
	memcpy( pd->_numSites, _numSites, sizeof(short)*_info.numGenes());
	return pd;
}

void PromoterData::clearData() {
	for(int i=0; i<_info.numGenes(); i++)
		_numSites[i] = 0;
	for(int i=0; i<_info.totalRegions(); i++) {
		_pool[i] = 0;
	}
}

void PromoterData::set(int i, PTYPE c) {
	//cout << "PromoterData::set " << i << ", " << c << " -> "<< _info.getGeneForRegion(i) << endl;
	if( c == 0 && _pool[i] > 0) { // Losing a binding site
		_numSites[_info.getGeneForRegion(i)] -= 1;		
	} else if( c > 0 && _pool[i] == 0) { // Gaining a binding site
		_numSites[_info.getGeneForRegion(i)] += 1;
	}
	
	_pool[i] = c;
	
}

void PromoterData::set(int i, int j, PTYPE c) {
	set( _info.offset(i) + j, c);
}

PTYPE PromoterData::get(int i)  { return _pool[i]; }

int PromoterData::totalRegions() const { return _info.totalRegions(); }

int PromoterData::numGenes() const { return _info.numGenes(); }

int PromoterData::numTFs() const { return _info.numTFs(); }

int PromoterData::numMotifs() const { return _info.numMotifs(); }

PTYPE PromoterData::getBinding(int i, int j)  { 
	return _pool[ _info.offset(i)+j]; 
}

STYPE PromoterData::numSitesForGene(int i) {
	return _numSites[i];
}

const GlobalInfo& PromoterData::info() const { return _info; }


