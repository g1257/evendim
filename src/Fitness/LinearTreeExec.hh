#ifndef LINEARTREEEXEC_HH
#define LINEARTREEEXEC_HH
#include "Vector.h"

namespace Gep {

template<typename ComplexOrRealType>
class LinearTreeExec {

public:

	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VecStringType = std::vector<std::string>;

	LinearTreeExec(const VecStringType& vecStr, SizeType /* threadNum */)
	{

	}

	RealType energy() const
	{
		return 0;
	}


};
}
#endif // LINEARTREEEXEC_HH
