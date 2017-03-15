#ifndef EXAMPLE_1_H
#define EXAMPLE_1_H

namespace Gep {

template<typename EvolutionType_>
class Example1 {

public:

	typedef EvolutionType_ EvolutionType;
	typedef typename EvolutionType::PrimitivesType PrimitivesType;
	typedef typename PrimitivesType::ValueType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	Example1(SizeType samples, const EvolutionType& evolution)
	    : samples_(samples),evolution_(evolution)
	{
		if (evolution.inputs() != 1) {
			throw PsimagLite::RuntimeError("Example1::ctor(): 1 input expected\n");
		}
		for (SizeType i = 0; i < samples; i++)
			samples_[i] = evolution_.primitives().rng() * 2.0 - 1.0;
	}

	template<typename SomeChromosomeType>
	RealType getFitness(const SomeChromosomeType& chromosome)
	{
		bool verbose = evolution_.verbose();
		RealType sum = 0;

		for (SizeType i = 0; i < maxFitness(); i++) {
			bool b = true; //(evolution_.rng() < 0.5) ? true : false;
			RealType r = evolution_.primitives().rng() * 2.0 - 1.0;
			RealType x = (b) ? r : samples_[i];
			samples_[i] = r;
			RealType fOfX = f(x);
			evolution_.setInput(0,x);

			if (verbose) evolution_.printInputs(std::cout);

			RealType tmp = fabs((chromosome.exec(0)-fOfX)/fOfX);

			sum += (1.0 - fabs(tmp));
		}
		return sum;
	}

	RealType maxFitness() const { return samples_.size(); }

private:

	RealType f(const RealType& x)
	{
		RealType tmp = x * (x - 1) * (x + 1);
		return tmp*tmp;
	}

	VectorRealType samples_;
	const EvolutionType& evolution_;
}; // class Example1

} // namespace Gep

#endif // EXAMPLE_1_H
