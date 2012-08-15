/*
 *  main.cpp
 *  GPPG
 *
 *  Created by Troy Ruths on 5/14/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */


#include <iostream>
using namespace std;

//#include "GreedyLoad.h"
#include "SequenceOperation.h"

using namespace GPPG;
using namespace GPPG::Model;

int main ()
{
	//GreedyLoad gl(3);
	SequenceData* sd = new SequenceData();

	SequenceOperationRoot dop(sd);
	
	SequenceData* sd2 = dop.evaluate();
	
	cout << "Hello World!";
	return 0;
}