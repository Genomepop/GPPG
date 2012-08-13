/*
 *  OperationHeap.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */


#include "OperationHeap.h"

using namespace GPPG;

template <typename T> OperationHeap<T>* OperationHeap<T>::m_pInstance = NULL;

template <typename T> OperationHeap<T>* OperationHeap<T>::InstancePtr() {
	if(!m_pInstance)
		m_pInstance = new OperationHeap<T>();
	return m_pInstance;
}

template <typename T> OperationHeap<T>& OperationHeap<T>::Instance() {
	return OperationHeap<T>::InstancePtr();
}

template <typename T> OperationHeap<T>::OperationHeap() : policy(0) {}

template <typename T> OperationHeap<T>::OperationHeap(OperationHeap<T> const&) {}

template <typename T> OperationHeap<T>& OperationHeap<T>::operator=(OperationHeap<T> const&) {}

