#ifndef EVENDIM_HAMILTONIAN_BASE_H
#define EVENDIM_HAMILTONIAN_BASE_H
#include "AllocatorCpu.h"

namespace Gep {

template <typename QuasiVectorType>
class HamiltonianBase {

public:

	using ComplexOrRealType = typename QuasiVectorType::value_type;
	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;

	virtual RealType energy(const QuasiVectorType& y, SizeType threadNum) const = 0;

	virtual SizeType numberOfSites() const = 0;
};
}
#endif // EVENDIM_HAMILTONIAN_BASE_H
