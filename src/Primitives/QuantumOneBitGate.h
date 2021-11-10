#ifndef QUANTUM_ONE_BIT_GATES_H
#define QUANTUM_ONE_BIT_GATES_H
#include "Node.h"
#include "Matrix.h"

namespace Gep {

template<typename ComplexOrRealType>
class OneBitGateLibrary {

public:

	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

	static void fillHadamard(MatrixType& gateMatrix)
	{
		static const ComplexOrRealType oneOverSqrt2 = 1/sqrt(2.);

		gateMatrix.resize(2, 2);
		gateMatrix(0, 0) = oneOverSqrt2;
		gateMatrix(0, 1) = oneOverSqrt2;
		gateMatrix(1, 0) = oneOverSqrt2;
		gateMatrix(1, 1) = -oneOverSqrt2;
	}

	static void fillPhase(MatrixType& gateMatrix)
	{
		gateMatrix.resize(2, 2);
		gateMatrix(0, 0) = 1;
		gateMatrix(1, 1) = ComplexOrRealType(0, 1);
	}

	// ind = 0 means rotation around x
	// ind = 1 means rotation around y
	// ind = 2 means rotation around z
	static void rotation(MatrixType& gateMatrix, SizeType ind, RealType angle)
	{
		const RealType cosine = cos(angle);
		const RealType sine = sin(angle);

		gateMatrix.resize(2, 2);
		if (ind == 0) {
			gateMatrix(0, 0) = cos(angle);
			gateMatrix(0, 1) = ComplexOrRealType(0, -sine);
			gateMatrix(1, 0) = ComplexOrRealType(0, -sine);
			gateMatrix(1, 1) = cos(angle);
			return;
		} else if (ind == 1) {
			gateMatrix(0, 0) = cosine;
			gateMatrix(0, 1) = -sine;
			gateMatrix(0, 1) = sine;
			gateMatrix(1, 1) = cosine;
			return;
		} else if (ind == 2) {
			gateMatrix(0, 0) = ComplexOrRealType(cosine, -sine);
			gateMatrix(0, 1) = 0;
			gateMatrix(0, 1) = 0;
			gateMatrix(1, 1) = ComplexOrRealType(cosine, sine);
			return;
		}
	}
}; // class GateLibrary

template<typename VectorValueType>
class QuantumOneBitGate : public Node<VectorValueType> {

public:

	typedef typename VectorValueType::value_type ValueType;
	typedef typename ValueType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

	QuantumOneBitGate(char c,
	           SizeType bitNumber,
	           SizeType numberOfBits,
	           const MatrixType& gateMatrix)
	    : code_("FF"),
	      bitNumber_(bitNumber),
	      gateMatrix_(gateMatrix),
	      w_(1 << numberOfBits)
	{
		code_[0] = c;
		assert(bitNumber < 10);
		PsimagLite::String c1 = ttos(bitNumber);
		assert(c1.length() == 1);
		code_[1] = c1[0];
		numberOfBits_ = numberOfBits;
	}

	virtual PsimagLite::String code() const { return code_; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);

		const ValueType& vv = v[0];
		const int n = vv.size();
		assert(n == (1 << numberOfBits_));  // 2^N

		std::fill(w_.begin(), w_.end(), 0);
		for (int i = 0; i < n; ++i) {
			SizeType j = findBasisState(i);
			SizeType bitI = getBitForIndex(i);
			SizeType bitJ = getBitForIndex(j);
			w_[i] += gateMatrix_(bitI, bitI)*vv[i];
			w_[j] += gateMatrix_(bitI, bitJ)*vv[i];
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
	PsimagLite::String code_;
	SizeType bitNumber_;
	MatrixType gateMatrix_;
	mutable ValueType w_;

}; // class QuantumOneBitGate

template<typename T>
SizeType QuantumOneBitGate<T>::numberOfBits_ = 0;
}

#endif // QUANTUM_ONE_BIT_GATES_H
