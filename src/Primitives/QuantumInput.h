#ifndef QUANTUMINPUT_H
#define QUANTUMINPUT_H
#include "Vector.h"

namespace Gep {

template<typename VectorValueType>
class QuantumInput : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

public:

	QuantumInput(SizeType numberOfBits)
	{
		 numberOfBits_ = numberOfBits;
	}

	virtual PsimagLite::String code() const { return "0"; }

	virtual SizeType arity() const { return 0; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 0);
		return input_;
	}

	virtual void set(const ValueType& x) const { input_ = x; }

	virtual bool isInput() const  { return true; }

private:

	static SizeType numberOfBits_;
	mutable ValueType input_;
};

template<typename T>
SizeType QuantumInput<T>::numberOfBits_ = 0;

}
#endif // QUANTUMINPUT_H
