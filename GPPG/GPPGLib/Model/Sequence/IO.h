/*
 *  IO.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/17/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef SEQUENCE_IO_
#define SEQUENCE_IO_

#include <iostream>
#include "Model/Sequence/Data.h"
#include "Model/Sequence/Operation.h"

std::ostream& operator<<(std::ostream& output, const GPPG::Model::SequenceData& s);	
std::ostream& operator<<(std::ostream& output, const GPPG::Model::SequencePointMutator& s);	

#endif