/*
 *  UbiOperationGraph.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/23/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef UBIOPERATION_GRAPH_
#define UBIOPERATION_GRAPH_

#include <Operation/OperationHeap.h>

namespace GPPG {
	
	class IOperation;
	class ICompressionPolicy;
	
class UbiOperationGraph : public OperationGraph {
public:
	UbiOperationGraph(ICompressionPolicy* p);
	~UbiOperationGraph();
	
	void addOperation(IOperation* op);
	
	void removeOperation(IOperation* op);
};
	
}
#endif