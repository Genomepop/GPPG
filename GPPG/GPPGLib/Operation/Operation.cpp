/*
 *  Operation.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "GPPG.h"
#include "Operation.h"
#include <sstream>
#include <string>
#include <iomanip>
#include <iostream>


#ifdef UBIGRAPH
extern "C" {
#include <Util/Ubigraph/ubiclient.h>
}
#endif

using namespace GPPG;
using std::string;

string ConvertRGBtoHex(int num) {
	static string hexDigits = "0123456789ABCDEF";
	string rgb;
	for (int i=(3*2) - 1; i>=0; i--) {
		rgb += hexDigits[((num >> i*4) & 0xF)];
	}
	return rgb;
}

string ConvertRGBtoHex(int r, int g, int b) {
	if (r > 255) r = 255;
	if (g > 255) g = 255;
	if (b > 255) b = 255;
	int rgbNum = ((r & 0xff) << 16)
	| ((g & 0xff) << 8)
	| (b & 0xff);
	
	return ConvertRGBtoHex(rgbNum);
}

template <class T> std::string TToStr( const T &t )
{
    std::ostringstream oss;
	oss << t;
    return std::string (oss.str());
}

std::ostream& operator<<(std::ostream& output, const GPPG::IOperation& op) {
	output << "Op (" << op.key() << ", " << op.frequency() << ") P:[";
	for (int i =0; i<op.numParents(); i++) {
		output << op.parent(i)->key() << ", ";
	}
	output << "]   C:[";
	const std::set<IOperation*>& children = op.children();
	for (std::set<IOperation*>::iterator it = children.begin(); it != children.end(); it++) {
		output << (*it)->key() << ", ";
	}
	output << "]";
	return output;
}

BaseOperation::BaseOperation(int cost) : 
	_key(-1), _index(-1), _state(-1), _freq(0.0), _total(0), _order(-1), _fitness(1.0), _cost(cost), _load(0), _loadFreq(0), _loadCost(0), _requests(0), _touch(0) {
#ifdef UBIGRAPH
	//ubigraph_new_vertex_w_id( (long)this );
#endif
}

BaseOperation::~BaseOperation() {
#ifdef UBIGRAPH
	ubigraph_remove_vertex( key() );
	//ubigraph_set_vertex_attribute( key(), "shape", "cone" );
#endif
	//std::cout << "DELETING " << key() << std::endl;
}

void BaseOperation::configure() {
#ifdef UBIGRAPH
	ubigraph_new_vertex_w_id( key() );
	for (int i=0; i<numParents(); i++) {
		int eid = (key() << 16) | parent(i)->key();
		ubigraph_new_edge_w_id(eid, parent(i)->key(), key() );
		//ubigraph_set_edge_attribute( eid, "width", TToStr<double>(this->cost()).c_str() );
	}
	
	if (isCompressed()) {
		
		if (_freq == 0) {
			ubigraph_change_vertex_style( key(), 1);	
		} else if( _freq == 0) {
			ubigraph_change_vertex_style( key(), 0);
		}
		
	}
	
#endif
}

int BaseOperation::key() const { return _key; }

void BaseOperation::setKey(int k) { _key = k; }

double BaseOperation::frequency() const { return _freq; }


void BaseOperation::setFrequency(double f) { 
	 
#ifdef UBIGRAPH
	if (isCompressed() && f != _freq) {
		
		if (f == 0) {
			ubigraph_change_vertex_style( key(), 1);	
		} else if( _freq == 0) {
			ubigraph_change_vertex_style( key(), 0);
		}

	}
	//ubigraph_set_vertex_attribute( key(), "size", TToStr<double>(10*_freq+1.0).c_str() );
	//int i = 255*_freq;
	//std::string rc = "#" + ConvertRGBtoHex( i, 255-i, 0);
	//ubigraph_set_vertex_attribute( key(), "color", rc.c_str() );
#endif
	_freq = f;
}


double BaseOperation::total() const { return _total; }


void BaseOperation::setTotal(double t) { _total = t; }

bool BaseOperation::isActive() const { return _freq > 0; }

int BaseOperation::index() const { return _index; }

void BaseOperation::setIndex(int i) { _index = i; }

int BaseOperation::state() const { return _state; }

void BaseOperation::setState(int i) { _state = i; }

int BaseOperation::order() const {
	return _order;
}

void BaseOperation::setOrder(int i) {
	_order = i;
}

void BaseOperation::setCompressed( bool c ) {
#ifdef UBIGRAPH
	if(key() < 0) return;
	if(c) {
		if (isActive()) ubigraph_change_vertex_style( key(), 0);
		else ubigraph_change_vertex_style( key(), 1);
		ubigraph_set_vertex_attribute( key(), "label", "" );
	} else {
		ubigraph_change_vertex_style( key(), 2);
		ubigraph_set_vertex_attribute( key(), "label", "X" );
	}
#endif
}

double BaseOperation::fitness() const { return _fitness; }
void BaseOperation::setFitness(double f) { _fitness = f; }

int BaseOperation::cost() const { return _cost; }

void BaseOperation::setCost(int v) { _cost = v; if(_cost<=0) _cost=1; }

int BaseOperation::requests() const { return _requests + _touch; }

void BaseOperation::clearRequests() {
	setRequests(0);
	_touch = 0;
}

void BaseOperation::clearDescendentRequests() {
	if(_requests == 0 && _touch == 0) return;
	
	clearRequests();
	const std::set<IOperation*>& childs = children();
	for(std::set<IOperation*>::iterator it=childs.begin(); it!=childs.end(); it++) {
		if( (*it)->isCompressed() && (*it)->requests()>0)
			(*it)->clearDescendentRequests();
	}
}

const char* BaseOperation::exportFormat() {
	return "Operation has no export formats";
}

void BaseOperation::setRequests(int i) {
	_requests = i;
	if(_requests < 0) _requests = 0;
	#ifdef UBIGRAPH
	double v = _requests/1.0;
	v = (v > 5) ? 5 : v+1;
	if(key() >= 0) ubigraph_set_vertex_attribute(key(), "size", TToStr<double>(v).c_str());
	
	#endif
}
void BaseOperation::incrRequests(int i) {
	setRequests(_requests + i*_cost);
}
void BaseOperation::decrRequests(int i) { 
	setRequests(_requests -i); 
}

void BaseOperation::touch() {
	if( _touch > 0) return;
	_touch = 1;
	if(isCompressed()) {
		int p = numParents();
		if(p>0) parent(0)->touch();
		if(p>1) parent(1)->touch();
	}
}
std::string BaseOperation::toString() const {
	return "<No Content>";
}