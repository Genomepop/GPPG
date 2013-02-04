#ifndef PATHWAY_FITNESS_
#define PATHWAY_FITNESS_

#include <Base/FitnessFunction.h>
#include <Model/Pathway/Operation.h>

namespace GPPG {
	
	namespace Model {
		namespace TransReg {
			
			class ConnectedFitness : public FitnessFunction<OpPathway> {
			public:
				virtual double calculate( OpPathway* g);
			};
		}
	}
}
 
#endif