#ifndef HAMILTONIANEXAMPLE_H
#define HAMILTONIANEXAMPLE_H
#include "InputNg.h"
#include "InputCheck.h"
#include "PsimagLite.h"

namespace Gep {

template<typename ComplexType>
class HamiltonianExample {

public:

	typedef PsimagLite::InputNg<InputCheck> InputNgType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;

	HamiltonianExample(typename InputNgType::Readable& io)
	{}

	RealType energy(const VectorType& v) const
	{
		throw PsimagLite::RuntimeError("energy(): unimplemented yet, sorry\n");
	}
};
}
#endif // HAMILTONIANEXAMPLE_H
