/*
 *  IO.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/17/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef PATHWAY_IO_
#define PATHWAY_IO_

#include <iostream>
#include "Model/Pathway/Data.h"
#include "Model/Pathway/Operation.h"

std::ostream& operator<<(std::ostream& output, const GPPG::Model::TransReg::GlobalInfo& info);	


#endif