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
#ifndef EXAMPLE_3_H
#define EXAMPLE_3_H
#include "BaseFitness.h"
#include "Vector.h"

namespace Gep {

/* PSIDOC Example3class
 Example3 illustrates the case of a training function with many variables
 and consists of a function f(x0, x1, ..., x5) of six variables.
 The variables are in the space of valid alphanumeric characters.
 If x0 is a digit then the function returns that digit plus one.
 If not, but if x1 is a digit then the function returns that digit plus one.
 And so on until all arguments to f are evaluated. If none of them are digits,
 then the function returns -1.
 */
template<typename EvolutionType_>
class Example3 : public BaseFitness<EvolutionType_> {

	static const SizeType stringLength_ = 6;

public:

	typedef BaseFitness<EvolutionType_> BaseType;
	typedef typename BaseType::FitnessParamsType FitnessParamsType;
	typedef EvolutionType_ EvolutionType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

	Example3(SizeType samples, const EvolutionType& evolution, FitnessParamsType*)
	    : samples_(samples),evolution_(evolution)
	{
		if (evolution.numberOfInputs() != stringLength_) {
			throw PsimagLite::RuntimeError("Example3::ctor(): 1 input expected\n");
		}
	}

	template<typename SomeChromosomeType>
	RealType getFitness(const SomeChromosomeType& chromosome,
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
				r[j] = validLetter();

			evolution_.setInput(r);

			if (verbose) evolution_.printInputs(std::cout);
			RealType fOfX = f(r);
			RealType tmp = fabs((chromosome.exec(0)-fOfX)/fOfX);

			sum += (1.0 - fabs(tmp));
		}

		return sum;
	}

	RealType maxFitness() const { return samples_; }

private:

	RealType f(const VectorRealType& r) const
	{
		for (SizeType i = 0; i < r.size(); ++i)
			if (r[i] >= 65 && r[i] <= 90) return i+1;

		return -1.0;
	}

	SizeType validLetter() const
	{
		SizeType l = 0;
		while ((l < 32) || (l > 126)) {
			l = static_cast<SizeType>(128*evolution_.rng());
		}

		return l;
	}

	SizeType samples_;
	const EvolutionType& evolution_;
}; // class Example3

} // namespace Gep

#endif // EXAMPLE_3_H
