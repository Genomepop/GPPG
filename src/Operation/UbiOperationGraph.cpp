/*
 *  UbiOperationGraph.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 8/23/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "UbiOperationGraph.h"
#include <Operation/Operation.h>
#include <Operation/CompressionPolicy.h>
#include <iostream>

extern "C" {
#include <Util/Ubigraph/ubiclient.h>
}

using namespace GPPG;
using std::cout;
using std::endl;

UbiOperationGraph::UbiOperationGraph(ICompressionPolicy* p) : OperationGraph(p) {
	ubigraph_clear();
}

UbiOperationGraph::~UbiOperationGraph() {}

void UbiOperationGraph::addOperation(IOperation* op) {
	OperationGraph::addOperation( op );
	//cout << op->toString() << endl;
	
	// Add Vertex
	ubigraph_new_vertex_w_id( (long)op );
	//cout << "Operation is " << typeid(op).name() << endl;
	//cout << "Done creating vertex, about to add " << op->numParents() << " edges.\n";

	// Connect parents
	for (int i=0; i<op->numParents(); i++) {
		ubigraph_new_edge( (long)op->parent(i), (long)op );
		
	}
}

void UbiOperationGraph::removeOperation(IOperation* op) {
	OperationGraph::removeOperation( op );
	cout << "Ubi RemoveOperation: " << op << endl;
}