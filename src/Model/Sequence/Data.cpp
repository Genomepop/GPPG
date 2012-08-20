/*
 *  SequenceData.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 8/17/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "Model/Sequence/Data.h"
#include <iostream>

using namespace GPPG::Model;


SequenceData* SequenceData::copy() const
{
	SequenceData* other = new SequenceData(_length);
	memcpy( other->sequence(), _sequence, sizeof(STYPE)*_length);
	return other;
}

SequenceData::SequenceData(int length) : _length(length) {
	_sequence = (STYPE*)malloc(sizeof(STYPE)*_length);
}

STYPE* SequenceData::sequence() { return _sequence; }

int SequenceData::length() const { return _length; }

STYPE SequenceData::get(int i) const { return _sequence[i]; }

void SequenceData::set(int i, STYPE c) { _sequence[i] = c; }

