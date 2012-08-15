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

BaseGenotype::BaseGenotype() {}

void BaseGenotype::configure() {}

double BaseGenotype::frequency() const { return _freq; }


void BaseGenotype::setFrequency(double f) { _freq = f; }


double BaseGenotype::total() const { return _total; }


void BaseGenotype::setTotal(double t) { _total = t; }

bool BaseGenotype::isActive() const { return _freq > 0; }

