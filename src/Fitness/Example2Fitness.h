/*
Copyright (c) 2017, UT-Battelle, LLC

evendim, Version 0.

This file is part of evendim.
evendim is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
evendim is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with evendim. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef EXAMPLE_2_FITNESS_H
#define EXAMPLE_2_FITNESS_H
#include "BaseFitness.h"

namespace Gep {

template <typename ChromosomeType>
class Example2Fitness : public BaseFitness<ChromosomeType> {

public:

	typedef BaseFitness<ChromosomeType> BaseType;
	typedef typename BaseType::FitnessParamsType FitnessParamsType;
	typedef typename ChromosomeType::EvolutionType EvolutionType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	Example2Fitness(SizeType samples, EvolutionType& evolution, FitnessParamsType*)
	    : samples_(samples)
	    , evolution_(evolution)
	{
		if (evolution.nodeHelper().numberOfInputs() != 1) {
			throw PsimagLite::RuntimeError("Example2Fitness::ctor(): 1 input expected\n");
		}
	}

	RealType getFitness(const ChromosomeType& chromosome,
	                    unsigned long int seed,
	                    SizeType threadNum)
	{
		bool verbose = evolution_.verbose();
		RealType sum = 0;
		PsimagLite::MersenneTwister rng(seed);
		for (SizeType i = 0; i < maxFitness(); i++) {
			SizeType x = static_cast<SizeType>(rng() * 1000);
			RealType fOfX = f(x);
			evolution_.nodeHelper().setInput(0, x, threadNum);
			if (verbose)
				evolution_.nodeHelper().printInputs(std::cout);

			RealType tmp = fabs((chromosome.exec(0) - fOfX) / fOfX);

			sum += (1.0 - fabs(tmp));
		}
		return sum;
	}

	RealType maxFitness() const { return samples_; }

private:

	RealType f(const SizeType& x)
	{
		if (x < 4)
			return 1;

		SizeType sqrtX = 1 + static_cast<SizeType>(sqrt(x));

		for (SizeType i = 2; i < sqrtX; i++) {
			if (x % i == 0)
				return -1;
		}

		return 1;
	}

	SizeType samples_;
	EvolutionType& evolution_;
}; // class Example2Fitness

} // namespace Gep

#endif // EXAMPLE_2_FITNESS_H
