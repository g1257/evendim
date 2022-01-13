#ifndef CANONICAL_FORM_EMPTY_H
#define CANONICAL_FORM_EMPTY_H
#include "Vector.h"

namespace Gep {

class CanonicalFormEmpty {

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	CanonicalFormEmpty(const VectorStringType&) {}

	void changeIfNeeded(VectorStringType&, VectorStringType&) const {}

};
}
#endif // CANONICAL_FORM_EMPTY_H
