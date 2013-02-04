/*
 *  Operation.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 9/24/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "Operation.h"



#include <iostream>
#include <algorithm>
#include <sstream>
#include "Util/Random.h"

using namespace GPPG;
using namespace GPPG::Model;
using namespace GPPG::Model::TransReg;
using std::vector;
using std::string;
using std::map;

template <class T> std::string TToStr( const T &t )
{
    std::ostringstream oss;
    oss << t;
    return std::string (oss.str());
}

OpPathwayBase::OpPathwayBase(double cost, const GlobalInfo& info, OpPathway& parent1) : 
	OpPathway(cost, parent1), _info(info) {
		
}

OpPathwayBase::OpPathwayBase(double cost, const GlobalInfo& info, OpPathway& parent1, OpPathway& parent2) :
OpPathway(cost, parent1, parent2), _info(info)  {
		
}

int OpPathwayBase::numGenes() const { return _info.numGenes(); }

int OpPathwayBase::numTFs() const { return _info.numTFs(); }

int OpPathwayBase::numMotifs() const { return _info.numMotifs(); }

int OpPathwayBase::totalRegions() const { return _info.totalRegions(); }

PTYPE OpPathwayBase::get(int i) {
	incrRequests(1);
	if (isCompressed()) {
		return proxyGet(i);
	}
	return data()->get(i);
}

short OpPathwayBase::numSitesForGene(int i) {
	incrRequests(1);
	if (isCompressed()) {
		return proxyNumSitesForGene(i);
	}
	return data()->numSitesForGene(i);
}

PTYPE OpPathwayBase::getBinding(int i, int j)  {
	return get( _info.offset(i) + j );
}

const GlobalInfo& OpPathwayBase::info() const { return _info; }

const char* OpPathwayBase::exportFormat() {
	std::ostringstream output;
	PromoterData* pd = evaluate();
	for(int i=0; i<numGenes(); i++) {
		// print gene
		output << _info.getGeneName(i) << "\t[" << i+1 << "]\t";
		/*
		std::vector<int> motifs = _info.binding(i);
		for (int j=0; j<motifs.size(); j++) {
			if(j>0) output << ",";
			output << motifs[i];
		}
		output << "]\t";
		*/
		for(int j=0; j<_info.numRegions(i); j++) {
			output << "|" << pd->getBinding(i,j);
		}
		output << std::endl;
	}
	
	return output.str().c_str();
}

/**
 ********************************** PATHWAY ROOT ****************************************
 */

PathwayRoot::PathwayRoot( PromoterData* p) : OperationRoot<PromoterData,ITransRegPathway>(p) {}

int PathwayRoot::numGenes() const { return data()->numGenes(); }

int PathwayRoot::numTFs() const { return data()->numTFs(); }

int PathwayRoot::numMotifs() const { return data()->numMotifs(); }

int PathwayRoot::totalRegions() const { return data()->totalRegions(); }

PTYPE PathwayRoot::get(int i)  { incrRequests(1); return data()->get(i); }

PTYPE PathwayRoot::getBinding(int i, int j)  { incrRequests(1); return data()->getBinding(i,j); }

short PathwayRoot::numSitesForGene(int i) { incrRequests(1); return data()->numSitesForGene(i); }

const GlobalInfo& PathwayRoot::info() const { return data()->info(); }

/**
 ********************************** PATHWAY FACTORY ****************************************
 */

PromoterData* randomPromoter(const GlobalInfo& info) {
	PromoterData* data = new PromoterData( info );
	data->clearData();
	int numMotifs = data->numMotifs();
	for (int i=0; i<info.numGenes(); i++) {
		data->set(i, 0, (PTYPE)(random01()*numMotifs)+1);
	}
	
	return data;
}

PathwayRootFactory::PathwayRootFactory( const GlobalInfo& info ) : _info(info) {}

PathwayRoot* PathwayRootFactory::random() const {
	return new PathwayRoot( randomPromoter( _info ) );
}

GlobalInfo* PathwayRootFactory::randomInfo(int numGenes, int numTFs, int minRegions, int maxRegions) {
	vector<int> regions, tfs;
	vector<string> genes, motifs;
	map< int, vector<int> > binding;
	
	for (int i=0; i<numGenes; i++) {
		regions.push_back( random01()*(maxRegions-minRegions) + minRegions );
		genes.push_back( "G" + TToStr<int>(i) );
	}
	
	// Populate a 1 to 1 mapping of motif to TF
	for (int i=0; i<numTFs; i++) {
		motifs.push_back( "M" + TToStr<int>(i) );
		tfs.push_back(i);
		binding[i].push_back(i);
	}
	
	GlobalInfo* info = new GlobalInfo( genes, regions, motifs, tfs, binding);
	return info;
}

/**
 ********************************** OPERATIONS *********************************************
 */
BindingSiteChange::BindingSiteChange(OpPathway& op, std::vector<int>* locs, std::vector<PTYPE>* dest) :
OpPathwayBase(1, op.info(), op), _locs(locs), _c(dest), _deltaSites(), _noSiteExists(false) {
	vector<int>::iterator it_loc = _locs->begin();
	vector<PTYPE>::iterator it_c = _c->begin();
	int gene;
	while (it_loc != _locs->end() ) {
		gene = info().getGeneForRegion( *it_loc );
		if( _deltaSites.count(gene) == 0) 
			_deltaSites[gene] = op.numSitesForGene(gene);
		
		if( *it_c == 0)
			_deltaSites[gene] -= 1;
		else
			_deltaSites[gene] += 1;
		it_loc++; it_c++;
	}
	for(map<int,short>::iterator it_map = _deltaSites.begin(); it_map != _deltaSites.end(); it_map++) {
		if(it_map->second == 0) {
			_noSiteExists = true;
		}
	}
	
}


BindingSiteChange::~BindingSiteChange() {
	delete _locs;
	delete _c;
}

PromoterData* BindingSiteChange::evaluate() {
	PromoterData* sd = OpPathwayBase::evaluate();
	if (sd != NULL) return sd;
	
	// Get the sequence from the parent and add the point changes
	sd = parent(0)->evaluate();
	
	vector<int>::iterator it_loc = _locs->begin();
	vector<PTYPE>::iterator it_c = _c->begin();
	while (it_loc != _locs->end() ) {
		sd->set( *it_loc, *it_c);
		it_loc++; it_c++;
	}
	return sd;
}

std::string BindingSiteChange::toString() const {
	std::ostringstream output;
	for (int i=0; i<numSites(); i++) {
		output << getSite(i) << "->" << getMutation(i);
		if (i < numSites()-1) output << ", ";
	}
	return output.str();
}

int BindingSiteChange::numSites() const { return _locs->size(); }

PTYPE BindingSiteChange::getMutation(int i) const { return (*_c)[i]; }

int BindingSiteChange::getSite(int i) const { return (*_locs)[i]; }

const std::map<int, short>& BindingSiteChange::deltaSites() const { return _deltaSites; }

PTYPE BindingSiteChange::proxyGet(int l)  {
	// See if the index is in the list
	vector<int>::iterator it_loc = _locs->begin();
	vector<PTYPE>::iterator it_c = _c->begin();
	while (it_loc != _locs->end() ) {
		if (*it_loc == l) return *it_c;
		it_loc++; it_c++;
	}
	return parent(0)->get(l);
}

short BindingSiteChange::proxyNumSitesForGene(int i) {
	if( _deltaSites.count(i) )
		return _deltaSites[i];

	return parent(0)->numSitesForGene(i);
}
bool BindingSiteChange::areAllGenesRegulated() const {
	return _noSiteExists;
}

BindingSiteMutator::BindingSiteMutator( double cost, double u, int motifOverlap, const vector<double>& motifGainRates, const vector<double>& motifProbLoss) :
OperationMutator< OpPathway >(cost), _u(u), _overlap(motifOverlap), _gainRates(motifGainRates), _lossProb(motifProbLoss) {
	
	
}

BindingSiteMutator::~BindingSiteMutator() {

}


OpPathway* BindingSiteMutator::mutate( OpPathway& g ) const {
	bool isCompressed = g.isCompressed();
	if( isCompressed ) g.setCompressed(false);
	
	int totalRegions = g.totalRegions();
	int numMotifs = g.numMotifs();
	const GlobalInfo& info = g.info();
	// Calculate the losses
	int numLosses = binomial( totalRegions, _u );
	int loc, minSite, maxSite;
	vector<int>* locs = new vector<int>();
	vector<PTYPE>* sites = new vector<PTYPE>();
	
	PTYPE c;
	int g_i, g_offset, g_numRegions;
	for (int i=0; i<numLosses; i++) {
		loc = (int)(random01()*totalRegions);
		minSite = loc-_overlap;
		maxSite = loc+_overlap;
		if (_overlap > 0) {
			g_i = info.getGeneForRegion(loc);
			g_offset = info.offset(	g_i );
			g_numRegions = info.numRegions(g_i)-1;
			if (minSite < g_offset) minSite = g_offset;
			if (maxSite > g_offset+g_numRegions) maxSite = g_offset+g_numRegions;
		}
		for (int site_i=minSite; site_i<maxSite+1; site_i++) {
			c = g.get( site_i );
			if (c>0 && random01() <= _lossProb[c]) {
				// Save site_i, ->0
				sites->push_back((PTYPE)0);
				locs->push_back(site_i);
			}
		}
	}
	// Loop through the motifs and calculate the number of gains
	int numGains;
	for (int i=0; i<numMotifs; i++) {
		numGains = binomial( totalRegions, _gainRates[i] );
		for (int j=0; j<numGains; j++) {
			loc = (int)(random01()*totalRegions);
			c = g.get( loc );
			if (c != (PTYPE)i) {
				// Save loc, ->i
				sites->push_back((PTYPE)i);
				locs->push_back(loc);
			}
		}
	}
	
	if( isCompressed )  g.setCompressed(true); 
	
	if( sites->size() == 0) {
		delete locs;
		delete sites;
		return &g;
	}
	
	// Create mutation
	BindingSiteChange* bsc = new BindingSiteChange(g, locs, sites);
	bsc->setCost( cost() );
	return bsc;
}

int BindingSiteMutator::numMutants(OpPathway& g, long N, double f) const {
	throw "Not Yet Implemented";
	return 1;
}

double BindingSiteMutator::rate() const {
	return _u;
}
