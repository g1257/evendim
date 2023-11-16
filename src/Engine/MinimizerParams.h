#ifndef MINIMIZERPARAMS_H
#define MINIMIZERPARAMS_H
#include "InputCheck.h"
#include "InputNg.h"
#include "Vector.h"

namespace Gep {

template <typename RealType_>
struct MinimizerParams {

	typedef RealType_ RealType;

	typedef PsimagLite::InputNg<InputCheck>::Readable InputNgReadableType;

	enum EnumAlgo { NONE,
		        SIMPLEX,
		        CONJUGATE_GRADIENT };

	MinimizerParams(EnumAlgo algo_,
	                SizeType maxIter_,
	                RealType delta_,
	                RealType delta2_,
	                RealType tol_,
	                SizeType saveEvery_,
	                bool verbose_)
	    : verbose(verbose_)
	    , algo(algo_)
	    , maxIter(maxIter_)
	    , saveEvery(saveEvery_)
	    , delta(delta_)
	    , delta2(delta2_)
	    , tol(tol_)
	{
	}

	MinimizerParams(InputNgReadableType& io, SizeType /*numberOfThreads*/)
	    : verbose(false)
	    , algo(EnumAlgo::CONJUGATE_GRADIENT)
	    , maxIter(100)
	    , saveEvery(0)
	    , delta(0.1)
	    , delta2(0.1)
	    , tol(1e-3)
	{
		try {
			PsimagLite::String algoString;
			io.readline(algoString, "MinimizerAlgorithm=");
			setAlgo(algoString);
		}
		catch (std::exception&) {
		}

		try {
			io.readline(maxIter, "MinimizerMaxIterations=");
		}
		catch (std::exception&) {
		}

		try {
			io.readline(saveEvery, "MinimizerSaveEvery=");
		}
		catch (std::exception&) {
		}

		try {
			io.readline(delta, "MinimizerDelta=");
		}
		catch (std::exception&) {
		}

		try {
			io.readline(delta2, "MinimizerDelta2=");
		}
		catch (std::exception&) {
		}

		try {
			io.readline(tol, "MinimizerTolerance=");
		}
		catch (std::exception&) {
		}

		try {
			int tmp = 0;
			io.readline(tmp, "MinimizerVerbose=");
			verbose = (tmp > 0);
		}
		catch (std::exception&) {
		}
	}

	void setAlgo(PsimagLite::String algoString)
	{
		if (algoString == "ConjugateGradient")
			algo = EnumAlgo::CONJUGATE_GRADIENT;
		else if (algoString == "Simplex")
			algo = EnumAlgo::SIMPLEX;
		else if (algoString == "None")
			algo = EnumAlgo::NONE;
		else
			err("setAlgo(): Unknown minimizer algorithm " + algoString + "\n");
	}

	bool verbose;
	EnumAlgo algo;
	SizeType maxIter;
	SizeType saveEvery;
	RealType delta;
	RealType delta2;
	RealType tol;
};

template <typename RealType>
std::ostream& operator<<(std::ostream& os, const MinimizerParams<RealType>& m)
{
	os << "algo= " << m.algo << "\n";
	os << "maxIter= " << m.maxIter << "\n";
	os << "delta= " << m.delta << "\n";
	os << "delta2= " << m.delta2 << "\n";
	os << "tolerance= " << m.tol << "\n";
	os << "verbose= " << m.verbose << "\n";
	return os;
}
} // namespace Gep

#endif // MINIMIZERPARAMS_H
