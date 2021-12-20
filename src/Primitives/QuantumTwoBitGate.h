#ifndef QUANTUM_TWO_BIT_GATE_H
#define QUANTUM_TWO_BIT_GATE_H
#include "Node.h"
#include "Matrix.h"

namespace Gep {

template<typename ComplexOrRealType>
class TwoBitGateLibrary {

public:

	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

	static void fillCnot(MatrixType& gateMatrix)
	{
		gateMatrix.resize(4, 4);
		gateMatrix(0, 0) = gateMatrix(0, 1) = 1;
		gateMatrix(2, 3) = gateMatrix(3, 2) = 1;
	}
}; // class GateLibrary

template<typename VectorValueType>
class QuantumTwoBitGate : public Node<VectorValueType,
        typename PsimagLite::Real<typename VectorValueType::value_type::value_type>::Type> {

public:

	typedef typename VectorValueType::value_type ValueType;
	typedef typename ValueType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef Node<VectorValueType,
	typename PsimagLite::Real<typename VectorValueType::value_type::value_type>::Type> NodeType;
	typedef typename NodeType::VectorAnglesType VectorAnglesType;

	QuantumTwoBitGate(PsimagLite::String cr,
	                  SizeType bitNumber1,
	                  SizeType bitNumber2,
	                  SizeType numberOfBits,
	                  const MatrixType& gateMatrix)
	    : code_(cr),
	      bitNumber1_(bitNumber1),
	      bitNumber2_(bitNumber2),
	      gateMatrix_(gateMatrix),
	      w_(1 << numberOfBits)
	{
		code_ += ttos(bitNumber1);
		code_ += "_";
		code_ += ttos(bitNumber2);

		numberOfBits_ = numberOfBits;
	}

	virtual PsimagLite::String code() const { return code_; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v,
	                       const VectorAnglesType*,
	                       SizeType&) const
	{
		return exec(v);
	}

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);

		const ValueType& vv = v[0];
		const int n = vv.size();
		assert(n == (1 << numberOfBits_));  // 2^N

		std::fill(w_.begin(), w_.end(), 0);
		const SizeType mask2 = (1 << bitNumber2_);
		for (int i = 0; i < n; ++i) {
			const SizeType oldContent1 = getBitForIndex(i, bitNumber1_);
			assert(oldContent1 < 2);
			const SizeType oldContent2 = getBitForIndex(i, bitNumber2_);
			assert(oldContent2 < 2);
			const SizeType content2 = (oldContent1 + oldContent2) % 2;
			assert(content2 < 2);

			const SizeType j = (content2 == oldContent2) ? i : (i ^ mask2);

			w_[j] += vv[i];
		}

		return w_;
	}

private:

	SizeType findBasisState(SizeType ind, SizeType content) const
	{
		SizeType mask = 0;
		switch (content) {
		case 0:
			return ind;
		case 1:
			mask = 1 << bitNumber1_;
			return ind ^ mask;
		case 2:
			mask = 1 << bitNumber2_;
			return ind;
		case 3:
			mask = 1 << bitNumber1_;
			mask |= 1 << bitNumber2_;
			return ind ^ mask;
		default:
			err("findBasisState\n");
		}

		throw PsimagLite::RuntimeError("findBasisState\n");
	}

	SizeType getTwoBitForIndex(SizeType ind) const
	{
		const SizeType bit1 = getBitForIndex(ind, bitNumber1_);
		const SizeType bit2 = getBitForIndex(ind, bitNumber2_);
		return bit1 + bit2*2;
	}

	static SizeType getBitForIndex(SizeType ind, SizeType bitNumber)
	{
		const SizeType mask = (1 << bitNumber);
		const SizeType result = ind & mask;
		return (result > 0) ? 1 : 0;
	}

	static SizeType numberOfBits_;
	PsimagLite::String code_;
	SizeType bitNumber1_;
	SizeType bitNumber2_;
	MatrixType gateMatrix_;
	mutable ValueType w_;

}; // class QuantumTwoBitGate

template<typename T>
SizeType QuantumTwoBitGate<T>::numberOfBits_ = 0;
}

#endif // QUANTUM_TWO_BIT_GATE_H
