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

// Seed for RNG
namespace GPPG {
	

boost::mt19937 gen;

inline int binomial(int n, double r) {
	boost::random::binomial_distribution<> dist( n, r );
	return dist(gen);
};

}
#endif
