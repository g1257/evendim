#ifndef QUANTUM_TWO_BIT_GATE_H
#define QUANTUM_TWO_BIT_GATE_H
#include "AST/Node.h"
#include "Matrix.h"

namespace Gep {

template <typename ComplexOrRealType>
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

template <typename VectorValueType>
class QuantumTwoBitGate : public PsimagLite::Node<VectorValueType,
                                                  typename PsimagLite::Real<typename VectorValueType::value_type::value_type>::Type> {

public:

	typedef typename VectorValueType::value_type ValueType;
	typedef typename ValueType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef PsimagLite::Node<VectorValueType,
	                         typename PsimagLite::Real<typename VectorValueType::value_type::value_type>::Type>
	    NodeType;
	typedef typename NodeType::VectorAnglesType VectorAnglesType;

	QuantumTwoBitGate(PsimagLite::String cr,
	                  SizeType bitNumber1,
	                  SizeType bitNumber2,
	                  SizeType numberOfBits,
	                  const MatrixType& gateMatrix)
	    : code_(cr)
	    , bitNumber1_(bitNumber1)
	    , bitNumber2_(bitNumber2)
	    , gateMatrix_(gateMatrix) // CNOT gate only has been implemented here
	{
		code_ += ttos(bitNumber1);
		code_ += "_";
		code_ += ttos(bitNumber2);

		numberOfBits_ = numberOfBits;
	}

	QuantumTwoBitGate* clone() const
	{
		return new QuantumTwoBitGate(*this);
	}

	virtual PsimagLite::String code() const { return code_; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v,
	                       const VectorAnglesType*,
	                       SizeType&) const
	{
		return exec(v);
	}

	// CNOT gate only has been implemented here
	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);

		const ValueType& vv = v[0];
		assert(vv.size() == static_cast<SizeType>(1 << numberOfBits_)); // 2^N
		return CNOT(vv, bitNumber1_, bitNumber2_);
	}

private:

	static SizeType numberOfBits_;
	PsimagLite::String code_;
	SizeType bitNumber1_;
	SizeType bitNumber2_;
	MatrixType gateMatrix_; // CNOT gate only has been implemented here
}; // class QuantumTwoBitGate

template <typename T>
SizeType QuantumTwoBitGate<T>::numberOfBits_ = 0;
}

#endif // QUANTUM_TWO_BIT_GATE_H
