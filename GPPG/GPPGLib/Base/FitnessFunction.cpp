#include "Base/FitnessFunction.h"
#include "Util/Random.h"

using namespace GPPG;

double UniformRandomFitnessFunction::calculate(IGenotype* g) {
	return random01();
}