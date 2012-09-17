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

BaseGenotype::BaseGenotype() : _index(-1), _freq(0), _total(0), _order(-1), _state(-1), _fitness(1.0) {}

int BaseGenotype::key() const { return (long)this; }

void BaseGenotype::configure() {}

double BaseGenotype::frequency() const { return _freq; }


void BaseGenotype::setFrequency(double f) { _freq = f; }


double BaseGenotype::total() const { return _total; }


void BaseGenotype::setTotal(double t) { _total = t; }

bool BaseGenotype::isActive() const { return _freq > 0; }

int BaseGenotype::index() const { return _index; }

void BaseGenotype::setIndex(int i) { _index = i; }

int BaseGenotype::state() const { return _state; }
void BaseGenotype::setState(int i) { _state = i; }

int BaseGenotype::order() const {
	return _order;
}

void BaseGenotype::setOrder(int i) {
	_order = i;
}

double BaseGenotype::fitness() const { return _fitness; }
void BaseGenotype::setFitness(double f) { _fitness = f; }