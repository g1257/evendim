#ifndef GROUNDSTATEPARAMS_H
#define GROUNDSTATEPARAMS_H
#include "InputCheck.h"
#include "InputNg.h"
#include "MinimizerParams.h"

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
		PsimagLite::String infile;
		io.readline(infile, "InVectorFile=");
		fillInVector(infile);
	}

	void fillInVector(PsimagLite::String)
	{
		err("fillInVector unimplemented, sorry!\n");
	}

	MinimizerParamsType minParams;
	HamiltonianType hamiltonian;
	VectorType inVector;
};

}
#endif // GROUNDSTATEPARAMS_H
