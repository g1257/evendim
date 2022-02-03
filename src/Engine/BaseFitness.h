#ifndef BASEFITNESS_H
#define BASEFITNESS_H
#include "PsimagLite.h"

namespace Gep {

class NullClass {};

template<typename EvolutionType>
class BaseFitness {
public:

	typedef NullClass FitnessParamsType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	const SizeType status() const { return 0; }

	template<typename SomeChromosomeType>
	PsimagLite::String info(const SomeChromosomeType&) const
	{
		return "";
	}

};
}
#endif // BASEFITNESS_H
