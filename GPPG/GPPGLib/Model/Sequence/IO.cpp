/*
 *  IO.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 8/17/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include "IO.h"

using namespace std;

ostream& operator<<(ostream& output, const GPPG::Model::SequenceData& s) {
	for (int i=0; i<s.length(); i++) output << s.get(i);
	return output;
}