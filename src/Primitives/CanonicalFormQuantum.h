#ifndef CANONICALFORMQUANTUM_H
#define CANONICALFORMQUANTUM_H
#include "Vector.h"

namespace Gep {

class CanonicalFormQuantum {

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	CanonicalFormQuantum(const VectorStringType&) {}

	void changeIfNeeded(VectorStringType&, VectorStringType&) const {}
};
}
#endif // CANONICALFORMQUANTUM_H
