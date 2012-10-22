/*
 *  CompressionPolicy.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "CompressionPolicy.h"
#include "Operation/Operation.h"
#include "Operation/OperationHeap.h"

using namespace GPPG;

void CompressionPolicy::operationAdded(IOperation* op) {}

void CompressionPolicy::operationRemoved( IOperation* op ) {}

void CompressionPolicy::decompressionReleased( IOperation* op ) {}

void CompressionPolicy::generationFinished( OperationGraph* heap, const std::vector<IOperation*>& ) {}

void CompressionPolicy::generationFinished( OperationGraph* heap, const std::set<IOperation*>& ) {}