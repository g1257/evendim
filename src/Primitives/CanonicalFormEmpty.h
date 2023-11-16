#ifndef CANONICAL_FORM_EMPTY_H
#define CANONICAL_FORM_EMPTY_H
#include "Vector.h"

namespace Gep {

class CanonicalFormEmpty {

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	template <typename T>
	CanonicalFormEmpty(const VectorStringType&, const T&) { }

	void changeIfNeeded(VectorStringType&) const { }
};
}
#endif // CANONICAL_FORM_EMPTY_H
