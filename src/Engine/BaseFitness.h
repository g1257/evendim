#ifndef BASEFITNESS_H
#define BASEFITNESS_H
#include "PsimagLite.h"
#include "MersenneTwister.h"

namespace Gep {

class NullClass {};

template<typename EvolutionType>
class BaseFitness {
public:

	typedef NullClass FitnessParamsType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::Vector<long unsigned int>::Type VectorLongUnsignedIntType;

	BaseFitness(long unsigned int seed) : rng_(seed) {}

	BaseFitness() : rng_(12344) {}

	const SizeType status() const { return 0; }

	template<typename SomeChromosomeType>
	PsimagLite::String info(const SomeChromosomeType&) const
	{
		return "";
	}

	VectorLongUnsignedIntType createSeeds(SizeType total)
	{
		VectorLongUnsignedIntType seeds(total);
		for (SizeType i = 0; i < total; ++i)
			seeds[i] = rng_.random();
		return seeds;
	}

	double rng() const { return rng_(); }

private:

	mutable PsimagLite::MersenneTwister rng_;
};
}
#endif // BASEFITNESS_H
