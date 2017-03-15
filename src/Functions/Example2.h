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
#ifndef EXAMPLE_2_H
#define EXAMPLE_2_H

namespace Gep {

template<typename EvolutionType_>
class Example2 {

public:

	typedef EvolutionType_ EvolutionType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	Example2(SizeType samples, const EvolutionType& evolution)
	    : samples_(samples),evolution_(evolution)
	{
		if (evolution.inputs() != 1) {
			throw PsimagLite::RuntimeError("Example2::ctor(): 1 input expected\n");
		}
	}

	template<typename SomeChromosomeType>
	RealType getFitness(const SomeChromosomeType& chromosome)
	{
		bool verbose = evolution_.verbose();
		RealType sum = 0;
		const PrimitivesType& primitives = evolution_.primitives();
		for (SizeType i = 0; i < maxFitness(); i++) {
			SizeType x = static_cast<SizeType>(primitives.rng()*1000);
			RealType fOfX = f(x);
			evolution_.setInput(0,x);
			if (verbose) evolution_.printInputs(std::cout);

			RealType tmp = fabs((chromosome.exec(0)-fOfX)/fOfX);

			sum += (1.0 - fabs(tmp));
		}
		return sum;
	}

	RealType maxFitness() const { return samples_; }

private:

	RealType f(const SizeType& x)
	{
		if (x<4) return 1;

		SizeType sqrtX = 1 + static_cast<SizeType>(sqrt(x));

		for (SizeType i = 2; i < sqrtX; i++) {
			if (x % i == 0) return -1;
		}

		return 1;
	}

	SizeType samples_;
	const EvolutionType& evolution_;
}; // class Example2

} // namespace Gep

#endif // EXAMPLE_2_H
