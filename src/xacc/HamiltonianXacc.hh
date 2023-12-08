#ifndef HAMILTONIAN_XACC_H
#define HAMILTONIAN_XACC_H

#include "../Fitness/HamiltonianBase.hh"

namespace Gep {

template <typename ComplexType>
class HamiltonianXacc : public HamiltonianBase<std::vector<ComplexType>> {

public:

	using VectorType = std::vector<ComplexType>;
	using RealType = typename PsimagLite::Real<ComplexType>::Type;

	RealType energy(const VectorType& y, SizeType threadNum) const
	{
		return 0;
	}

	SizeType numberOfSites() const
	{
		return 0;
	}
};

}
#endif // HAMILTONIAN_XACC_H
