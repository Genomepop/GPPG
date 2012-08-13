/*
 *  GPPG.cp
 *  GPPG
 *
 *  Created by Troy Ruths on 5/11/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */

#include <iostream>
#include "GPPG.h"
#include "GPPGPriv.h"

void GPPG::HelloWorld(const char * s)
{
	 GPPGPriv *theObj = new GPPGPriv;
	 theObj->HelloWorldPriv(s);
	 delete theObj;
};

void GPPGPriv::HelloWorldPriv(const char * s) 
{
	std::cout << s << std::endl;
};

