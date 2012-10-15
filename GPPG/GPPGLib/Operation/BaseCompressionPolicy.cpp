/*
 *  BaseCompressionPolicy.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 8/23/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "BaseCompressionPolicy.h"
#include "Operation/Operation.h"

using namespace GPPG;

BaseCompressionPolicy::BaseCompressionPolicy(StoreFlag flag) : _flag(flag) {}

BaseCompressionPolicy::~BaseCompressionPolicy() {
	
}

StoreFlag BaseCompressionPolicy::storeFlag() const {
	return _flag;
}

void BaseCompressionPolicy::operationAdded(IOperation* op) {

	switch (_flag) {
		case STORE_ROOT:
			op->setCompressed(true);
			break;
		default:
			op->setCompressed(false);
			break;
	}

}

void BaseCompressionPolicy::operationRemoved(IOperation* op) {
	if (_flag == STORE_ACTIVE) {
		op->setCompressed(true);
	}
}