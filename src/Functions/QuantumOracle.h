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
#include "Minimizer.h"
#include "MinimizerParams.h"
#include "BaseFitness.h"

namespace Gep {

template<typename ComplexOrRealType>
class FunctionToMinimize {
public:
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef RealType FieldType;
};

template<typename EvolutionType_>
class QuantumOracle {

public:

	typedef EvolutionType_ EvolutionType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType VectorType;
	typedef typename VectorType::value_type ComplexType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef FunctionToMinimize<ComplexType> FunctionToMinimizeType;
	typedef typename PsimagLite::Minimizer<RealType, FunctionToMinimizeType> MinimizerType;
	typedef MinimizerParams<RealType> MinimizerParamsType;
	typedef MinimizerParamsType ExtraParamsType;

	QuantumOracle(SizeType samples,
	              const EvolutionType& evolution,
	              const MinimizerParamsType& minParams)
	    : samples_(samples),
	      evolution_(evolution),
	      minParams_(minParams),
	      inVector_((1 << evolution.primitives().numberOfBits())),
	      outVector_(inVector_.size())
	{
		if (evolution.inputs() != 1)
			err("QuantumOracle::ctor(): 1 input expected\n");
	}

	template<typename SomeChromosomeType>
	RealType getFitness(const SomeChromosomeType& chromosome)
	{
		const bool verbose = evolution_.verbose();

		RealType sum = 0;
		VectorType outputFromIndividual(inVector_.size());
		for (SizeType i = 0; i < maxFitness(); i++) {
			fillRandomVector();
			evolution_.setInput(0, inVector_);
			if (verbose) evolution_.printInputs(std::cout);

			functionF(outVector_, inVector_);

			findAnglesMax(outputFromIndividual, chromosome);
			RealType tmp = vectorDiff2(outputFromIndividual, outVector_);

			sum += (1.0 - fabs(tmp));
		}

		return sum/samples_;
	}

	RealType maxFitness() const { return samples_; }

private:

	template<typename SomeChromosomeType>
	int findAnglesMax(VectorType& output, const SomeChromosomeType& chromosome)
	{
		FunctionToMinimizeType f;
		MinimizerType min(f, minParams_.maxIter, minParams_.verbose);

		int used = 0;
		VectorRealType angles;
		if (minParams_.algo == MinimizerParamsType::SIMPLEX) {
			used = min.simplex(angles,
			                   minParams_.delta,
			                   minParams_.tol);
		} else {
			used = min.conjugateGradient(angles,
			                             minParams_.delta,
			                             minParams_.delta2,
			                             minParams_.tol,
			                             minParams_.saveEvery);
		}

		const bool printFooter = minParams_.verbose;
		const int returnStatus = (used > 0) ? 0 : 1;

		if (!printFooter) return returnStatus; // <--- EARLY EXIT HERE

		std::cerr<<"QuantumOracle::minimize(): ";
		if (min.status() == MinimizerType::GSL_SUCCESS) {
			std::cerr<<" converged after ";
		} else {
			std::cerr<<"NOT CONVERGED after ";
		}

		++used;
		std::cerr<<used<<" iterations.\n";

		return returnStatus;
	}

	// Flip the first bit
	// 0.1*|0000> -0.2|1110>
	// src[0] = 0.1;   src[14] = -0.2 src[..] = 0
	// 0.1*|0001> - 0.2|1111>
	// dest[1] = 0.1;  dest[15] = -0.2 dest[...] = 0
	// Flip the first bit
	// 1111 <--- 15 --> i
	// 0001 <--- 1
	// 1110 <--- 14 --> j
	static void functionF(VectorType& dest, const VectorType& src)
	{
		const SizeType n = dest.size();
		assert(n == src.size());
		for (SizeType i = 0; i < n; ++i) {
			SizeType j = i ^ 1;
			dest[j] = src[i];
		}
	}

	void fillRandomVector()
	{
		const SizeType n = inVector_.size();
		ComplexType sum = 0;
		for (SizeType i = 0; i < n; ++i) {
			inVector_[i] = 2.0*evolution_.primitives().rng() - 1.0;
			sum += inVector_[i]*PsimagLite::conj(inVector_[i]);
		}

		RealType factor = 1.0/sqrt(PsimagLite::real(sum));
		for (SizeType i = 0; i < n; ++i)
			inVector_[i] *= factor;
	}

	static RealType vectorDiff2(const VectorType& v1, const VectorType& v2)
	{
		const SizeType n = v1.size();
		assert(n == v2.size());
		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i)
			sum += fabs(v2[i] - v1[i]);

		return sum/n;
	}

	SizeType samples_;
	const EvolutionType& evolution_;
	const MinimizerParamsType minParams_;
	VectorType inVector_;
	VectorType outVector_;
}; // class QuantumOracle

} // namespace Gep

#endif // QUANTUM_ORACLE_H
