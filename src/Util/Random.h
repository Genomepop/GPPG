/*
 *  Random.h
 *  Demo
 *
 *  Created by Troy Ruths on 8/21/12.
 *  Copyright 2012 Rice University. All rights reserved.
 *
 */
#ifndef RANDOM_
#define RANDOM_

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <Util/binomial.h>

// Seed for RNG
namespace GPPG {
	

boost::mt19937 gen;

inline int binomialb(int n, double r) {
	boost::random::binomial_distribution<> dist( n, r );
	return dist(gen);
};
	
	inline double random01() { return ranmar(); }
	inline long binomial(long n, double pp) { return ignbin(n, pp); }
	void initRandom(int a, int b) { rmarin(a,b); }
	void initRandom() {
		initRandom((int)time(0), (int)time(0));
	}
}
#endif
