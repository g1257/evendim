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

	static SizeType directionCharToInteger(char c)
	{
		int val = c - 120;
		if (val >= 0 && val < 3) return val;

		throw PsimagLite::RuntimeError("findDirectionOfRotation\n");
	}

	static char directionIntegerToChar(SizeType ind)
	{
		if (ind < 3) {
			char c = ind + 120;
			return c;
		}

		throw PsimagLite::RuntimeError("findDirectionOfRotation\n");
	}
}; // class GateLibrary

template<typename VectorValueType>
class QuantumOneBitGate : public Node<VectorValueType,
        typename PsimagLite::Real<typename VectorValueType::value_type::value_type>::Type> {

public:

	typedef typename VectorValueType::value_type ValueType;
	typedef typename ValueType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef OneBitGateLibrary<ComplexOrRealType> OneBitGateLibraryType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	QuantumOneBitGate(PsimagLite::String cr,
	                  SizeType bitNumber,
	                  SizeType numberOfBits,
	                  const MatrixType& gateMatrix)
	    : code_(cr),
	      bitNumber_(bitNumber),
	      gateMatrix_(gateMatrix),
	      w_(1 << numberOfBits)
	{
		code_ += ttos(bitNumber);
		numberOfBits_ = numberOfBits;
	}

	virtual PsimagLite::String code() const { return code_; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v,
	                       const VectorRealType* angles,
	                       SizeType& currentIndex) const
	{
		if (!hasAngles()) return exec(v);

		assert(currentIndex < angles->size());
		SizeType angleToUse = angles->operator[](currentIndex++);
		assert(code_.size() > 1);
		SizeType ind = OneBitGateLibraryType::directionCharToInteger(code_[1]);
		OneBitGateLibraryType::rotation(gateMatrix_, ind, angleToUse);
		return exec(v);
	}

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

	bool hasAngles() const
	{
		assert(code_.size() > 0);
		return (code_[0] == 'R');
	}

	static SizeType numberOfBits_;
	PsimagLite::String code_;
	SizeType bitNumber_;
	mutable MatrixType gateMatrix_;
	mutable ValueType w_;

}; // class QuantumOneBitGate

template<typename T>
SizeType QuantumOneBitGate<T>::numberOfBits_ = 0;
}

#endif // QUANTUM_ONE_BIT_GATES_H
