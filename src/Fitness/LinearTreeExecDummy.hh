#ifndef LINEARTREEEXEC_DUMMY_HH
#define LINEARTREEEXEC_DUMMY_HH
#include "Vector.h"

namespace Gep {

template <typename ComplexOrRealType>
class LinearTreeExec {

public:

	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VecStringType = std::vector<std::string>;

	LinearTreeExec(const VecStringType& vecStr, SizeType /* threadNum */)
	{
		throw std::runtime_error("Please enable XACC to use LinearTreeExec\n");
	}

	RealType energy() const
	{
		throw std::runtime_error("Please enable XACC to use LinearTreeExec\n");
	}
};
}
#endif // LINEARTREEEXEC_DUMMY_HH
