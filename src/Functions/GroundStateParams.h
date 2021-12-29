#ifndef GROUNDSTATEPARAMS_H
#define GROUNDSTATEPARAMS_H
#include "InputCheck.h"
#include "InputNg.h"
#include "MinimizerParams.h"
#include "ProgramGlobals.h"

namespace Gep {

template<typename HamiltonianType, typename ComplexType_>
struct GroundStateParams {

	typedef ComplexType_ ComplexType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef PsimagLite::InputNg<InputCheck> InputNgType;
	typedef MinimizerParams<RealType> MinimizerParamsType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;

	GroundStateParams(typename InputNgType::Readable& io)
	    : minParams(io), hamiltonian(io)
	{
		PsimagLite::String vectorFilename;
		io.readline(vectorFilename, "InVectorFile=");
		Gep::readVector(inVector, vectorFilename);
		SizeType bits = 0;
		io.readline(bits, "NumberOfBits=");
		const SizeType hilbert = (1<<bits);
		if (hilbert != inVector.size())
			err("Initial vector has " + ttos(inVector.size()) +
			    " entries, but I was expecting " + ttos(hilbert) + "\n");

	}

	MinimizerParamsType minParams;
	HamiltonianType hamiltonian;
	VectorType inVector;
};

}
#endif // GROUNDSTATEPARAMS_H
