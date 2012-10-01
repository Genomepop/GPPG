/*
 *  Random.cpp
 *  Demo
 *
 *  Created by Troy Ruths on 10/1/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */


#include "Random.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <Util/binomial.h>

// Seed for RNG
		
boost::mt19937 gen2;
	
int GPPG::binomialb(int n, double r) {
	boost::random::binomial_distribution<> dist( n, r );
	return dist(gen2);
}
	
double GPPG::random01() { return ranmar(); }

long GPPG::binomial(long n, double pp) { return ignbin(n, pp); }

void GPPG::initRandom2(int a, int b) { 
	rmarin(a,b); 
}

void GPPG::initRandom() {
	rmarin((int)time(0), (int)time(0));
	gen2.seed((unsigned int)time(0));
}
