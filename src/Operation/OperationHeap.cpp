/*
 *  OperationHeap.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "GPPG.h"
#include "OperationHeap.h"
#include "Operation/Operation.h"
#include "Base/Genotype.h"
#include "Operation/CompressionPolicy.h"

#include <iostream>
#include <list>

#ifdef UBIGRAPH
extern "C" {
#include <Util/Ubigraph/ubiclient.h>
}
#endif

using namespace GPPG;
using std::cout;
using std::endl;
using std::list;

OperationGraph::OperationGraph(ICompressionPolicy* p) : _policy(p) {
	
#ifdef UBIGRAPH
	ubigraph_clear();
	// Set Properties
	ubigraph_set_edge_style_attribute(0, "oriented", "true");
	//ubigraph_set_edge_style_attribute(0, "color", "#C5892F");
	//ubigraph_set_edge_style_attribute(0, "color", "#CCCC2F");
	ubigraph_set_edge_style_attribute(0, "color", "#FFFFFF");
	ubigraph_set_edge_style_attribute(0, "width", "2.0");
	ubigraph_set_vertex_style_attribute( 0, "shape", "none" );

	// for compressed, active Operations
	ubigraph_set_vertex_style_attribute( 0, "shape", "sphere" );
	//ubigraph_set_vertex_style_attribute( 0, "color", "#80219C" );
	ubigraph_set_vertex_style_attribute( 0, "color", "#00FF00" );
	
	// For compressed,inactive Operations
	ubigraph_new_vertex_style_w_id(1, 0);
	ubigraph_set_vertex_style_attribute( 1, "shape", "none" );
	ubigraph_set_vertex_style_attribute( 1, "color", "#00FF00" );
	
	// For uncompressed Operations
	ubigraph_new_vertex_style_w_id(2, 0);
	ubigraph_set_vertex_style_attribute( 2, "shape", "cube" );
	ubigraph_set_vertex_style_attribute( 2, "color", "#FF0000" );
#endif
}

OperationGraph::~OperationGraph() {
	delete _policy;
}

ICompressionPolicy& OperationGraph::compressionPolicy() {
	return *_policy;
}

void OperationGraph::setHeapRemovalPolicy(HeapRemovalPolicy policy) {}

HeapRemovalPolicy OperationGraph::heapRemovalPolicy() const {
	return DELETE;
}


void OperationGraph::addGenotype(IGenotype* g) {
	addOperation( (IOperation*)g );
}


void OperationGraph::removeGenotype(IGenotype* g) {
	removeOperation( (IOperation*)g );
}


void OperationGraph::addOperation(IOperation* op) {
	//cout << "Add Operation: " << op << endl;
#ifdef UBIGRAPH
	//ubigraph_change_vertex_style( (long)op, 0);
#endif
	_policy->operationAdded( op );
}

/*
void OperationGraph::removeOperation(IOperation* op) {
	_policy->operationRemoved( op );
	
	IOperation* start = op;

#ifdef UBIGRAPH
	//ubigraph_change_vertex_style( (long)op, 1);
#endif
	
	// If this operation has children, abort
	if (start->numChildren() > 0 || start->numParents() == 0) return;
		
	if (1) { //! (start->isCompressed())) {
		// Release compression of start
		_policy->decompressionReleased( start );
	}
	IOperation* p_end = start;
	IOperation* end = start->parent(0);
	
	while (end != 0 && end->numChildren() == 1 && end->index() < 0 && end->numParents() > 0) {

		if (1) { //! (end->isCompressed())) {
			// Release compression
			_policy->decompressionReleased( end );
		}
		
		p_end = end;
		end = end->parent(0);
		
	}
	
	delete p_end;
}
 */

#define QUEUED -2
#define INACTIVE -1

void OperationGraph::removeOperation(IOperation* op) {
	_policy->operationRemoved( op );
	
	list<IOperation*> ops;
	ops.push_front(op);
	op->setState(QUEUED);
	
	IOperation* wop, *pop;
	
	while (ops.size() > 0) {
		wop = ops.front();
		ops.pop_front();
		wop->setState(INACTIVE);
#ifdef UBIGRAPH
		//ubigraph_set_vertex_attribute( wop->key(), "shape", "cone" );
		//sleep(0.1);
#endif
		// If this operation has children, is active, or is the root, bail!
		//if (wop != 0 && wop->numChildren() == 0 && !wop->isActive() && wop->numParents() > 0) {		
		if (wop != 0 && wop->numChildren() == 0 && wop->index() < 0 && wop->numParents() > 0) {
			for (int i=0; i<wop->numParents(); i++) {
				pop = wop->parent(i);
				if (pop->state() != QUEUED) {
					pop->setState(QUEUED);
					ops.push_front(pop);
				}
			}
			
			_policy->decompressionReleased(wop);
			delete wop;
		}
	}	
}

void OperationGraph::generationFinished(const std::vector<IGenotype*>& genos) {
	//_policy->generationFinished( (const std::vector<IOperation*>&) genos );
	// Convert vector to set...
}

void OperationGraph::generationFinished(const std::set<IGenotype*>& genos) {
	_policy->generationFinished( (const std::set<IOperation*>&) genos );
}
