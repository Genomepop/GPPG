/*
 *  Operation.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 9/24/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "Operation.h"

using namespace GPPG::Model;
using namespace GPPG::Model::TransReg;

#include <iostream>
#include <algorithm>
#include <sstream>

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

