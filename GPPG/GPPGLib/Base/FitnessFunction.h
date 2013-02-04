#ifndef FITNESS_FUNCTION_
#define FITNESS_FUNCTION_

namespace GPPG {
	class IGenotype;
	
	class IFitnessFunction {
	public:
		virtual ~IFitnessFunction(){}
		
		virtual double calculate( IGenotype* g) = 0;
		
	};
	
	class UniformRandomFitnessFunction : public IFitnessFunction {
	public:
		virtual double calculate( IGenotype* g);
	};
	
	template <typename T> class FitnessFunction : public IFitnessFunction {
		
	public:
		virtual double calculate( IGenotype* g) { calculate( (T*) g); }
		virtual double calculate(T* g) = 0;
	};
	
}
#endif