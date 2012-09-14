/*
 *  GreedyLoad.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "GreedyLoad.h"
#include "Operation/Operation.h"
#include "GPPG.h"
#include <sstream>

#ifdef UBIGRAPH
extern "C" {
#include <Util/Ubigraph/ubiclient.h>
}
#endif

using namespace GPPG;

using std::set;

typedef set<IOperation*>::iterator OpIter;

void annotate( const std::set<IOperation*>& active );
void innerAnnotate( IOperation* op, double freq, double cost );
void resetAnnotation( const std::set<IOperation*>& active);
void reset(IOperation* op);
IOperation* findMaxAdvance( IOperation* op, bool doReset );
IOperation* uncoveredChild( IOperation* op );

template <class T> std::string TToStr( const T &t )
{
    std::ostringstream oss;
    oss << t;
    return std::string (oss.str());
}

GreedyLoad::GreedyLoad(int maxExplicit, int numGens) : _root(0),_maxExplicit(maxExplicit), _waitGens(numGens), _elapsedGens(0), _numExplicit(0) {}

void GreedyLoad::decompressionReleased( IOperation* op ) {
#ifdef DEBUG_0
	std::cout << "Releasing decompression... " << op->key() << " count=" << _U.count( op ) << std::endl;
	for (OpIter it=_U.begin(); it!=_U.end(); it++) {
		std::cout << (*it)->key() << ",";
	}
	std::cout << std::endl;
#endif
	_U.erase( op );
}

void GreedyLoad::operationAdded( IOperation* op) {

}



void GreedyLoad::generationFinished( const std::set<IOperation*>& active ) {
	_elapsedGens++;
	if (_elapsedGens > _waitGens) {
		apply( active );
	}
}

void GreedyLoad::apply( const std::set<IOperation*>& active ) {
	_elapsedGens = 0;
	
#ifdef UBIGRAPH_GL
	ubigraph_set_vertex_attribute( 0, "label", "Applying Greedy Load" );
#endif
	// Step 0: Initialize with root
	if (_root == 0) {
		IOperation* op = *active.begin();
		while (op->numParents() > 0) op = op->parent(0);
		_root = op;
	}
	if (_U.size() == 0) {
		add( _root );
	}
	
	// Step 1: Annotate the Tree
#ifdef UBIGRAPH_GL
	ubigraph_set_vertex_attribute( 0, "label", "Annotating Load" );
#endif
	annotate( active );
	
	// Step 2: Remove defunct nodes
#ifdef UBIGRAPH_GL
	ubigraph_set_vertex_attribute( 0, "label", "Remove Defunct" );
#endif
	OpIter it=_U.begin();
	while (it!=_U.end()) {
		IOperation* op = *it;
		it++;
		if (op->load() == 0 && op != _root) 
			remove(op, false, false);
	}
	
	update();
	
	// Step 3: Advance down
#ifdef UBIGRAPH_GL
	ubigraph_set_vertex_attribute( 0, "label", "Advance" );
#endif
	it=_U.begin();
	while (it!=_U.end()) {
		IOperation* op = *it;
		it++;
		if (op != _root) 
			advance( op );
	}
	
	update();
	
	// Step 4: Split
	// Make a copy of the set
	set<IOperation*> C = _U;
	int s1, s2;
	IOperation* g1, *g2;
	
#ifdef UBIGRAPH_GL
	ubigraph_set_vertex_attribute( 0, "label", "Split" );
#endif
	while (_U.size() < _maxExplicit && C.size() > 0) {
		IOperation* op = getMaxItem( C, false );
		if (op == 0) break;
		
#ifdef UBIGRAPH_GL
		ubigraph_set_vertex_attribute( op->key(), "label", "Spliting..." );
#endif
		split( op, s1, s2, g1, g2 );
		
		if (s1 == 0) {
			C.erase( op );
		} else {
			C.insert( g1 );
			if (s2 == 1) {
				C.erase( op );
				C.insert( g2 );
			}
		}
#ifdef UBIGRAPH_GL
		//usleep(1000000);		
		ubigraph_set_vertex_attribute( op->key(), "label", "" );
#endif
	
	}
	
	
	// Step 5: Apply compression
	
	// Step 6: Reset
	resetAnnotation( active );
	for (OpIter it=_U.begin(); it!=_U.end(); it++) {
		reset( *it );
	
		
#ifdef UBIGRAPH_GL
		//ubigraph_set_vertex_attribute( (*it)->key(), "label", TToStr<int>((*it)->key()).c_str() );
		//usleep(1000000);
#endif
	}
	
#ifdef UBIGRAPH_GL
	ubigraph_set_vertex_attribute( 0, "label", "" );
#endif
}

void GreedyLoad::split( IOperation* op, int& s1, int& s2, IOperation*& g1, IOperation*& g2 ) {
	g2 = 0;
	s2 = 0;
	
	IOperation* c = getMaxItem( op->children(), true );
	if (c == 0) {
		s1 = 0;
		g1 = 0;
		return;
	}
	double r = c->load();
	c = findMaxAdvance( c, true );
	add( c ); //false
	op->setLoad( op->load() - r );
	g1 = c;
	s1 = 1;
	
	if (!op->isActive()) {
		c = uncoveredChild( op );
		if (c != 0) {
			reset( op );
			c = findMaxAdvance( c, true );
			move( op, c );
			g2 = c;
			s2 = 1;
		}
	}
}

void GreedyLoad::add( IOperation* op ) {
	_U.insert( op );
	op->setCompressed(false);
#ifdef UBIGRAPH
	ubigraph_set_vertex_attribute( op->key(), "label", "U" );
#endif
}

void GreedyLoad::remove( IOperation* op, bool cache, bool doRecurse ) {
	_U.erase( op );
	op->setCompressed(true);
#ifdef UBIGRAPH
	ubigraph_set_vertex_attribute( op->key(), "label", "" );
#endif
}

void GreedyLoad::update( ) {
	
}

void GreedyLoad::clearCache() {
	
}

void GreedyLoad::move( IOperation* a, IOperation* b ) {
	add( b );
	if (a != _root) {
		remove( a, true, false );
	}
}

IOperation* uncoveredChild( IOperation* op ) {
	IOperation* comp = 0;
	int numCompressed = 0;
	IOperation* child;
	const set<IOperation*>& children = op->children();
	for (OpIter it = children.begin(); it!=children.end(); it++) {
		child = *it;
		if (child->isCompressed() && (child->load() > 0 || child->isActive())) {
			numCompressed++;
			comp = child;
			if (numCompressed > 1) break;
		}
	}
	
	if (numCompressed == 1) return comp;

	return 0;
}

IOperation* GreedyLoad::getMaxItem( const std::set<IOperation*>& items, bool compare ) {
	double amt;
	double maxval = 0;
	IOperation* maxitem = 0;
	IOperation* op;
	for (OpIter it=items.begin(); it!=items.end(); it++) {
		op = *it;
		if (compare && _U.count( op ) > 0) continue;
		amt = op->load();
		if(amt > maxval) {
			maxval = amt;
			maxitem = op;
		}
	}
	return maxitem;
}

IOperation* findMaxAdvance( IOperation* op, bool doReset ) {
	IOperation* child = uncoveredChild( op );
	while (!op->isActive() && child != 0) {
		// Move it down
		if (doReset) reset(op);
		op = child;
		child = uncoveredChild( op );
	}
	return op;
}

void GreedyLoad::advance( IOperation* op ) {
	IOperation* t = findMaxAdvance( op, true );
	if (t != op) {
		reset( op );
		move( op, t );
	}
}

void annotate( const std::set<IOperation*>& active) {
	for (OpIter it=active.begin(); it!=active.end(); it++) {
		IOperation* op = *it;
		innerAnnotate( op, op->frequency(), 1 );
	}
}

void innerAnnotate( IOperation* op, double freq, double cost ) {
	op->setLoad( op->load()+freq*cost );

	if (op->isCompressed() && op->numParents() > 0) {
		innerAnnotate( op->parent(0), freq, cost+op->cost() );
	}
}

void reset(IOperation* op) {
	if (op->load() > 0 || op->isActive()) {
		op->setLoad(0);
		for (int i=0; i<op->numParents(); i++) {
			reset( op->parent(i) );
		}
	}
}

void resetAnnotation( const std::set<IOperation*>& active) {

	for (OpIter it=active.begin(); it!=active.end(); it++) {
		reset( *it );
	}
	
	
}

int GreedyLoad::maxUncompressed() const { return _maxExplicit; }

int GreedyLoad::numGenerations() const { return _waitGens; }