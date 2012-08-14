/*
 *  Genotype.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 8/10/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "Genotype.h"

using namespace GPPG;

template <typename T>
Genotype<T>::Genotype(T* genoData) : _data(genoData), _freq(0), _total(0) {}

template <typename T>
double Genotype<T>::frequency() const { return _freq; }

template <typename T>
void Genotype<T>::setFrequency(double f) { _freq = f; }

template <typename T>
double Genotype<T>::total() const { return _total; }

template <typename T>
void Genotype<T>::setTotal(double t) { _total = t; }

template <typename T>
T* Genotype<T>::data() const { return _data; }

// TODO: Add some checks here.
template <typename T>
void Genotype<T>::setData(T* data) { _data = data; }