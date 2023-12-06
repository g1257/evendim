#ifndef HAMILTONIANDUMMY_HH
#define HAMILTONIANDUMMY_HH

#include "../Fitness/HamiltonianBase.hh"

namespace Gep {

template <typename ComplexType>
class HamiltonianDummy : public HamiltonianBase<std::vector<ComplexType>> {

	using VectorType = std::vector<ComplexType>;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;

	RealType energy(const VectorType& y, SizeType threadNum) const
	{
		return 0;
	}
};

}
#endif // HAMILTONIANDUMMY_HH
