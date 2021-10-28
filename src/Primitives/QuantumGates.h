#ifndef QUANTUMGATES_H
#define QUANTUMGATES_H
#include "Node.h"
#include "Matrix.h"

namespace Gep {

template<typename VectorValueType>
class Hadamard : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;
	typedef typename ValueType::value_type ComplexOrRealType;

public:

	Hadamard(SizeType bitNumber, SizeType numberOfBits)
	    : bitNumber_(bitNumber), gateMatrix_(2, 2), w_(1 << numberOfBits)
	{
		static const ComplexOrRealType oneOverSqrt2 = 1/sqrt(2.);

		numberOfBits_ = numberOfBits;
		gateMatrix_(0, 1) = oneOverSqrt2;
		gateMatrix_(0, 2) = oneOverSqrt2;
		gateMatrix_(1, 0) = oneOverSqrt2;
		gateMatrix_(1, 1) = -oneOverSqrt2;
	}

	virtual PsimagLite::String code() const { return "H"; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		const int n = v.size();
		assert(n == (1 << numberOfBits_));  // 2^N

		std::fill(w_.begin(), w_.end(), 0);
		for (int i = 0; i < n; ++i) {
			SizeType j = findBasisState(i);
			SizeType bitI = getBitForIndex(i);
			SizeType bitJ = getBitForIndex(j);
			w_[i] += gateMatrix_(bitI, bitI);
			w_[j] += gateMatrix_(bitI, bitJ);
		}

		return w_;
	}

private:

	SizeType findBasisState(SizeType ind) const
	{
		const SizeType mask = (1 << bitNumber_);
		return ind ^ mask;
	}

	SizeType getBitForIndex(SizeType ind) const
	{
		const SizeType mask = (1 << bitNumber_);
		const SizeType result = ind & mask;
		return (result > 0) ? 1 : 0;
	}

	static SizeType numberOfBits_;
	SizeType bitNumber_;
	PsimagLite::Matrix<ComplexOrRealType> gateMatrix_;
	mutable ValueType w_;

}; // class Hadamard

template<typename T>
SizeType Hadamard<T>::numberOfBits_ = 0;
}

#endif // QUANTUMGATES_H
