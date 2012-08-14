/*
 *  DNAOperation.h
 *  GPPG
 *
 *  Created by Troy Ruths on 5/14/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#ifndef DNA_OPERATION_
#define DNA_OPERATION_

#include "Operation.h"

using namespace GPPG;

namespace Model {

/**
 * A simple data structure for holding sequence information.
 */
struct DNAData {
	short* seq;
	int length;
};

class DNAOperation : public Operation<DNAData> {
public:
	DNAOperation();
	
protected:
	DNAData* evaluate() const;
};

}
#endif