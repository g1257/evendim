/*
Copyright (c) 2017-2021, UT-Battelle, LLC

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
#ifndef QUANTUM_ORACLE_H
#define QUANTUM_ORACLE_H
#include "PsimagLite.h"

namespace Gep {

template<typename EvolutionType_>
class QuantumOracle {

public:

	typedef EvolutionType_ EvolutionType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType VectorRealType; // or complex FIXME TODO
	typedef typename VectorRealType::value_type RealType;

	QuantumOracle(SizeType samples, const EvolutionType& evolution)
	    : samples_(samples),
	      evolution_(evolution),
	      inVector_((1 << evolution.primitives().numberOfBits())),
	      outVector_(inVector_.size())
	{
		if (evolution.inputs() != 1)
			err("QuantumOracle::ctor(): 1 input expected\n");
	}

	template<typename SomeChromosomeType>
	RealType getFitness(const SomeChromosomeType& chromosome)
	{
		bool verbose = evolution_.verbose();
		RealType sum = 0;
		for (SizeType i = 0; i < maxFitness(); i++) {
			fillRandomVector();
			evolution_.setInput(0, inVector_);
			if (verbose) evolution_.printInputs(std::cout);

			// RealType tmp = fabs((chromosome.exec(0)-fOfX)/fOfX);

			RealType tmp = 0;
			err("getFitness: unimplemented\n");
			sum += (1.0 - fabs(tmp));
		}
		return sum;
	}

	RealType maxFitness() const { return samples_; }

private:

	void fillRandomVector()
	{
		const SizeType n = inVector_.size();
		for (SizeType i = 0; i < n; ++i)
			inVector_[i] = 2.0*evolution_.primitives().rng() - 1.0;
	}

	SizeType samples_;
	const EvolutionType& evolution_;
	VectorRealType inVector_;
	VectorRealType outVector_;
}; // class QuantumOracle

} // namespace Gep

#endif // QUANTUM_ORACLE_H
