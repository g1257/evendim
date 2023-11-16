#ifndef BASEFITNESS_H
#define BASEFITNESS_H
#include "MersenneTwister.h"
#include "PsimagLite.h"

namespace Gep {

class NullClass {
};

/* PSIDOC BaseFitness
The BaseFitness class provides an interface that fitness classes must follow.
It does also provide some basic non-virtual functionality.
PSIDOCCOPY getFitness
PSIDOCCOPY maxFitness
*/
template <typename ChromosomeType>
class BaseFitness {
public:

	typedef NullClass FitnessParamsType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::Vector<long unsigned int>::Type VectorLongUnsignedIntType;
	typedef double RealType;

	BaseFitness(long unsigned int seed)
	    : rng_(seed)
	{
	}

	BaseFitness()
	    : rng_(12344)
	{
	}

	/* PSIDOC getFitness
PSIDOCCOPY $FirstProtoBelow
Returns the fitness of the passed chromosome, where seed is a seed
for a random number generator, and threadNum is the number of the
thread in case fitness is computed in parallel.
*/
	virtual RealType getFitness(const ChromosomeType& chromosome,
	                            long unsigned int seed,
	                            SizeType threadNum)
	    = 0;

	/* PSIDOC maxFitness
PSIDOCCOPY $FirstProtoBelow
Returns the maximum fitness for this problem.
*/
	virtual RealType maxFitness() const = 0;

	const SizeType status() const { return 0; }

	virtual PsimagLite::String info(const ChromosomeType&) const
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
