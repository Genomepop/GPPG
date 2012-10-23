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
#include "Operation/OperationHeap.h"
#include "GPPG.h"
#include <sstream>
#include <map>
#include <Util/Random.h>

#ifdef UBIGRAPH
extern "C" {
#include <Util/Ubigraph/ubiclient.h>
}
#endif

using namespace GPPG;

using std::set;
using std::map;

typedef set<IOperation*>::iterator OpIter;

template <class T> std::string TToStr( const T &t )
{
    std::ostringstream oss;
    oss << t;
    return std::string (oss.str());
}


GreedyLoad::GreedyLoad(int maxExplicit, int numGens) : 
	_root(0),_maxExplicit(maxExplicit), _waitGens(numGens), _elapsedGens(0), _numExplicit(0), _runs(0) {}

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



void GreedyLoad::generationFinished( OperationGraph* heap, const std::set<IOperation*>& active ) {
	_elapsedGens++;

	if (_elapsedGens == _waitGens) {
		//heap->clearRequests();
	}
	else if (_elapsedGens > _waitGens) {
		apply( active );
		// Clear data requests
		//heap->clearRequests();
	}
	/*
	if( true ) {
		apply( active );
		heap->clearRequests();
	}*/
}

void GreedyLoad::apply( const std::set<IOperation*>& active ) {
	_elapsedGens = 0;
	
#ifdef UBIGRAPH_GL
	ubigraph_set_vertex_attribute( 0, "label", "Applying Greedy Load" );
#endif
	// Copy the initial U set to do the compression more efficiently
	_V = _U;
	
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
		if (load(op) == 0 && op != _root) 
			remove(op);
	}
	
	// Step 3: Advance down
#ifdef UBIGRAPH_GL
	//ubigraph_set_vertex_attribute( 0, "label", "Advance" );
#endif
	it=_U.begin();
	while (it!=_U.end()) {
		IOperation* op = *it;
		it++;
		if (op != _root) 
			advance( op );
	}
	
	// Step 4: Split
	// Make a copy of the set
	set<IOperation*> C = _U;
	int s1, s2;
	IOperation* g1, *g2;
	
#ifdef UBIGRAPH_GL
	//ubigraph_set_vertex_attribute( 0, "label", "Split" );
#endif
	while (_U.size() < _maxExplicit && C.size() > 0) {
		IOperation* op = getMaxItem( C, false );
		if (op == 0) {
			std::cout << "No max item found\n";
			break;
		}
		
#ifdef UBIGRAPH_GL
		//ubigraph_set_vertex_attribute( op->key(), "label", "Spliting..." );
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
		//ubigraph_set_vertex_attribute( op->key(), "label", "" );
#endif
	
	}
	
	
	// Step 5: Apply compression
	/*
	for (OpIter it=_V.begin(); it!=_V.end(); it++) 
		ubigraph_set_vertex_attribute( (*it)->key(), "label", "V" );
			
	for (OpIter it=_U.begin(); it!=_U.end(); it++) {
		if( _runs > 5 ) usleep(100000);
		ubigraph_set_vertex_attribute( (*it)->key(), "label", "Z" );		
		(*it)->setCompressed(false);
		if (_runs > 5 && _root->requests() > 0)
			ubigraph_set_vertex_attribute( (*it)->key(), "label", "P" );
		else 
			ubigraph_set_vertex_attribute( (*it)->key(), "label", "W" );
	}
	
	for (OpIter it=_V.begin(); it!=_V.end(); it++) 
		if (_U.count( *it ) == 0) {
			(*it)->setCompressed(true);
			ubigraph_set_vertex_attribute( (*it)->key(), "label", "" );
		}
	*/

	for (OpIter it=_U.begin(); it!=_U.end(); it++)
		(*it)->setCompressed(false);
	for (OpIter it=_V.begin(); it!=_V.end(); it++) 
		if (_U.count( *it ) == 0) (*it)->setCompressed(true);

	// Step 6: Reset
	clearLoadMap();
	
	_runs++;
	
	//resetAnnotation( active );
	//resetAnnotation(_U);
	//for (OpIter it=_U.begin(); it!=_U.end(); it++) {
	//	reset( *it );
	
		
#ifdef UBIGRAPH_GL
		//ubigraph_set_vertex_attribute( (*it)->key(), "label", TToStr<int>((*it)->key()).c_str() );
		//usleep(1000000);
#endif
	//}
	
#ifdef UBIGRAPH_GL
	ubigraph_set_vertex_attribute( 0, "label", "Simulating" );
#endif
}

void GreedyLoad::clearLoadMap() {
	for (OpIter it=_U.begin(); it!=_U.end(); it++)
		(*it)->clearDescendentRequests();
}

void GreedyLoad::setLoad(IOperation* op, int v) {
	if( v == 0) op->clearRequests();
	else op->setRequests(v);
}

void GreedyLoad::incrLoad(IOperation* op, int v) {
	op->incrRequests( v );
}

void GreedyLoad::decrLoad(IOperation* op, int v) {
	op->decrRequests( v );
}


int GreedyLoad::load(IOperation* op) {
	return op->requests();
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
	
	//op->setLoad( op->load() - c->load() - op->cost()*c->loadFreq(), op->loadFreq()-c->loadFreq() );
	int l = load(op)-load(c);
	op->setRequests( (l<=0)?1:l);
	
	c = findMaxAdvance( c, true );
	add( c ); //false
	g1 = c;
	s1 = 1;

	if (!op->isActive()) {
		c = uncoveredChild( op );
		if (c != 0) {
			reset(op);
			//resetAnnotation( op, true );
			c = findMaxAdvance( c, true );
			move( op, c );
			g2 = c;
			s2 = 1;
		}
	}
}

void GreedyLoad::add( IOperation* op ) {
	//resetAnnotation(op, false);
	_U.insert( op );
	//op->setCompressed(false);
#ifdef UBIGRAPH
	ubigraph_set_vertex_attribute( op->key(), "label", "U" );
#endif
}

void GreedyLoad::remove( IOperation* op ) {
	_U.erase( op );
	//op->setCompressed(true);
#ifdef UBIGRAPH
	ubigraph_set_vertex_attribute( op->key(), "label", "" );
#endif
}

void GreedyLoad::move( IOperation* a, IOperation* b ) {
	if( a == b) return;
	
	add( b );
	if (a != _root) {
		remove( a );
	}
}

bool GreedyLoad::isCompressed(IOperation* op) {
	return _U.count(op) == 0;
}

IOperation* GreedyLoad::uncoveredChild( IOperation* op ) {
	IOperation* comp = 0;
	int numCompressed = 0;
	IOperation* child;
	const set<IOperation*>& children = op->children();
	for (OpIter it = children.begin(); it!=children.end(); it++) {
		child = *it;
		if (isCompressed(child) && (load(child) > 0 || child->isActive())) {
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
		amt = load(op);
		if(amt > maxval) {
			maxval = amt;
			maxitem = op;
		}
	}
	return maxitem;
}

IOperation* GreedyLoad::findMaxAdvance( IOperation* op, bool doReset ) {
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
		//resetAnnotation( op, true );
		move( op, t );
	}
}


void GreedyLoad::annotate(const std::set<IOperation*>& active) {
	for (OpIter it=active.begin(); it!=active.end(); it++) {
		IOperation* op = *it;
		op->touch();
	}
}


/*
void GreedyLoad::resetOp(IOperation* op) {
	setLoad(op,0,0);
}
*/
void GreedyLoad::reset(IOperation* op) {
	op->clearRequests();
}


int GreedyLoad::maxUncompressed() const { return _maxExplicit; }

int GreedyLoad::numGenerations() const { return _waitGens; }