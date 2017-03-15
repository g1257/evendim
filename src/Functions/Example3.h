#ifndef EXAMPLE_3_H
#define EXAMPLE_3_H

#include "Vector.h"

namespace Gep {

template<typename EvolutionType_>
class Example3 {

	static const SizeType stringLength_ = 6;

public:

	typedef EvolutionType_ EvolutionType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

	Example3(SizeType samples, const EvolutionType& evolution)
	    : samples_(samples),evolution_(evolution)
	{
		if (evolution.inputs() != stringLength_) {
			throw PsimagLite::RuntimeError("Example3::ctor(): 1 input expected\n");
		}
	}

	template<typename SomeChromosomeType>
	RealType getFitness(const SomeChromosomeType& chromosome)
	{
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
		while (l < 32 || l > 126) {
			l = static_cast<SizeType>(128*evolution_.primitives().rng());
		}

		return l;
	}

	SizeType samples_;
	const EvolutionType& evolution_;
}; // class Example3

} // namespace Gep

#endif // EXAMPLE_3_H
