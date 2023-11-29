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
#ifndef EXAMPLE_3_FITNESS_H
#define EXAMPLE_3_FITNESS_H
#include "BaseFitness.h"
#include "Vector.h"

namespace Gep {

/* PSIDOC Example3FitnessClass
 Example3Fitness illustrates the case of a training function with many variables
 and consists of a function f(x0, x1, ..., x5) of six variables.
 The variables are in the space of valid alphanumeric characters.
 If x0 is a digit then the function returns that digit plus one.
 If not, but if x1 is a digit then the function returns that digit plus one.
 And so on until all arguments to f are evaluated. If none of them are digits,
 then the function returns -1.
 */
template <typename ChromosomeType>
class Example3Fitness : public BaseFitness<ChromosomeType> {

	static const SizeType stringLength_ = 6;

public:

	typedef BaseFitness<ChromosomeType> BaseType;
	typedef typename BaseType::FitnessParamsType FitnessParamsType;
	typedef typename ChromosomeType::EvolutionType EvolutionType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

	Example3Fitness(SizeType samples, const EvolutionType& evolution, FitnessParamsType*)
	    : samples_(samples)
	    , evolution_(evolution)
	{
		if (evolution.nodeHelper().numberOfInputs() != stringLength_) {
			throw PsimagLite::RuntimeError("Example3Fitness::ctor(): " + ttos(stringLength_) + " inputs expected\n");
		}
	}

	RealType getFitness(const ChromosomeType& chromosome,
	                    long unsigned int,
	                    SizeType threadNum)
	{
		if (threadNum > 0)
			err("Threading not supported yet (sorry)\n");

		bool verbose = evolution_.verbose();
		RealType sum = 0;

		VectorRealType r(stringLength_);
		for (SizeType i = 0; i < samples_; i++) {
			for (SizeType j = 0; j < stringLength_; ++j)
				r[j] = static_cast<SizeType>(128 * evolution_.rng());

			evolution_.nodeHelper().setInput(r);

			if (verbose)
				evolution_.nodeHelper().printInputs(std::cout);
			RealType fOfX = f(r);
			RealType tmp = fabs((chromosome.exec(0) - fOfX) / fOfX);

			sum += (1.0 - fabs(tmp));
		}

		return sum;
	}

	RealType maxFitness() const { return samples_; }

private:

	RealType f(const VectorRealType& r) const
	{
		RealType sum = 0;
		for (SizeType i = 0; i < r.size(); ++i)
			sum += r[i];

		return sum;
	}

	SizeType samples_;
	const EvolutionType& evolution_;
}; // class Example3Fitness

} // namespace Gep

#endif // EXAMPLE_3_FITNESS_H
