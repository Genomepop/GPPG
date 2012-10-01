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

PTYPE OpPathwayBase::get(int i) const {
	if (isCompressed()) {
		return proxyGet(i);
	}
	return data()->get(i);
}

PTYPE OpPathwayBase::getBinding(int i, int j) const {
	throw "Not Implemented";
	return 0;
}

const GlobalInfo& OpPathwayBase::info() const { return _info; }

/**
 ********************************** PATHWAY ROOT ****************************************
 */

PathwayRoot::PathwayRoot( PromoterData* p) : OperationRoot<PromoterData,ITransRegPathway>(p) {}

int PathwayRoot::numGenes() const { return data()->numGenes(); }

int PathwayRoot::numTFs() const { return data()->numTFs(); }

int PathwayRoot::numMotifs() const { return data()->numMotifs(); }

int PathwayRoot::totalRegions() const { return data()->totalRegions(); }

PTYPE PathwayRoot::get(int i) const { return data()->get(i); }

PTYPE PathwayRoot::getBinding(int i, int j) const { return data()->getBinding(i,j); }

const GlobalInfo& PathwayRoot::info() const { return data()->info(); }

/**
 ********************************** PATHWAY FACTORY ****************************************
 */

PromoterData* randomPromoter(const GlobalInfo& info) {
	PromoterData* data = new PromoterData( info );
	int numMotifs = data->numMotifs();
	for (int i=0; i<info.numGenes(); i++) {
		for (int j=0; j<info.numRegions(i); j++) {
			if (j == 0) 
				data->set(i, j, (PTYPE)(random01()*numMotifs)+1);
			else
				data->set(i, j, 0);
		}
	}
	return data;
}

PathwayRootFactory::PathwayRootFactory( const GlobalInfo& info ) : _info(info) {}

PathwayRoot* PathwayRootFactory::random() const {
	return new PathwayRoot( randomPromoter( _info ) );
}

/**
 ********************************** OPERATIONS *********************************************
 */
BindingSiteChange::BindingSiteChange(OpPathway& op, int* locs, int numLocs, PTYPE* dest) :
OpPathwayBase(1, op.info(), op), _loc(locs), _numlocs(numLocs), _c(dest) {}


BindingSiteChange::~BindingSiteChange() {
	delete _loc;
	delete _c;
}

PromoterData* BindingSiteChange::evaluate() const {
	PromoterData* sd = OpPathwayBase::evaluate();
	if (sd != NULL) return sd;
	
	// Get the sequence from the parent and add the point changes
	sd = parent(0)->evaluate();
	
	for (int i=0; i<_numlocs; i++) {
		sd->set(_loc[i], _c[i]);
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

int BindingSiteChange::numSites() const { return _numlocs; }

PTYPE BindingSiteChange::getMutation(int i) const { return _c[i]; }

int BindingSiteChange::getSite(int i) const { return _loc[i]; }


PTYPE BindingSiteChange::proxyGet(int l) const {
	// See if the index is in the list
	for (int i=0; i<_numlocs; i++) {
		if (_loc[i] == l) return _c[i];
	}
	return parent(0)->get(l);
}


BindingSiteMutator::BindingSiteMutator( double u, int motifOverlap, const vector<double>& motifGainRates, const vector<double>& motifProbLoss) :
_u(u), _overlap(motifOverlap), _gainRates(motifGainRates), _lossProb(motifLossProb) {
	
	_bufSize = 10000;
	_bufLoc = (int*)malloc(sizeof(int)*_bufSize);
	_bufC = (PTYPE*)malloc(sizeof(PTYPE)*_bufSize);
}

BindingSiteMutator::~BindingSiteMutator() {
	delete _bufLoc;
	delete _bufC;
}

OpPathway* BindingSiteMutator::mutate( OpPathway& g ) const {
	int totalRegions = g.totalRegions();
	int numMotifs = g.numMotifs();
	// Calculate the losses
	int numLosses = binomial( totalRegions, _u );
	int loc, minSite, maxSite;
	PTYPE c;
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
		for (int site_i=minSite; site_i<max_site+1; site_i++) {
			c = g->get( site_i );
			if (c>0 && random01() <= _lossProb[c]) {
				// TODO: Save site_i, ->0

			}
		}
	}
	// Loop through the motifs and calculate the number of gains
	int numGains;
	for (int i=0; i<numMotifs; i++) {
		numGains = binomial( totalRegions, _gainRates[i] );
		for (int j=0; j<numGains; j++) {
			loc = (int)(random01()*totalRegions);
			c = g->get( loc );
			if (c != (PTYPE)i) {
				// TODO: Save loc, ->i
			}
		}
	}
}

int BindingSiteMutator::numMutants(OpPathway& g, long N, double f) const {
	throw "Not Yet Implemented";
	return 1;
}

double BindingSiteMutator::rate() const {
	return _u;
}
