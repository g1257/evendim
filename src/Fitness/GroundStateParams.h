#ifndef GROUNDSTATEPARAMS_H
#define GROUNDSTATEPARAMS_H
#include "../Primitives/QuasiVector.hh"
#include "InputCheck.h"
#include "InputNg.h"
#include "MinimizerParams.h"
#include "ProgramGlobals.h"

namespace Gep {

template <typename HamiltonianType, typename ComplexType_, bool HAS_XACC>
struct GroundStateParams {

	typedef ComplexType_ ComplexType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef PsimagLite::InputNg<InputCheck> InputNgType;
	typedef MinimizerParams<RealType> MinimizerParamsType;
	typedef QuasiVector<ComplexType> QuasiVectorType;

	GroundStateParams(typename InputNgType::Readable& io, SizeType numberOfThreads)
	    : useXaccOptimizer(false)
	    , minParams(io, numberOfThreads)
	    , hamiltonian(io, numberOfThreads)
	{
		PsimagLite::String vectorFilename;
		io.readline(vectorFilename, "InVectorFile=");
		inVector.fromFile(vectorFilename);

		SizeType bits = 0;
		io.readline(bits, "NumberOfBits=");
		const SizeType hilbert = (1 << bits);
		if (hilbert != inVector.size())
			err("Initial vector has " + ttos(inVector.size()) + " entries, but I was expecting " + ttos(hilbert) + "\n");

		try {
			std::string tmp;
			io.readline(tmp, "UseXaccOptimizer=");
			if (tmp == "yes" || tmp == "true") {
				useXaccOptimizer = true;
				if (!HAS_XACC) {
					err("UseXaccOptimizer cannot be true without XACC support\n");
				}
			}
		}
		catch (std::exception&) {
		}
	}

	bool useXaccOptimizer;
	MinimizerParamsType minParams;
	HamiltonianType hamiltonian;
	QuasiVectorType inVector;
};

}
#endif // GROUNDSTATEPARAMS_H
