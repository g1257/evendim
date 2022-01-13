#ifndef CANONICALFORMQUANTUM_H
#define CANONICALFORMQUANTUM_H
#include "Vector.h"
#include "Node.h"

namespace Gep {

template<typename ValueType_, typename AnglesType>
class CanonicalFormQuantum {

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Vector<ValueType_>::Type VectorValueType;
	typedef Node<VectorValueType, AnglesType> NodeType;
	typedef typename PsimagLite::Vector<NodeType*>::Type VectorNodeType;

	CanonicalFormQuantum(const VectorStringType& data,
	                     const VectorNodeType& nodes)
	    : data_(data) , needsChange_(false)
	{
//		needsChange_ |= orderGatesByBit();
	}

	void changeIfNeeded(VectorStringType&, VectorStringType&) const {}

private:

//	bool orderGatesByBit()
//	{
//		for (auto it = data_.begin(); it != data_.end(); ++it) {

//			SizeType bit = getBit(*it);
//			v[bit].push_back(it - data_.begin());
//		}
//	}

	VectorStringType data_;
	bool needsChange_;
};
}
#endif // CANONICALFORMQUANTUM_H
