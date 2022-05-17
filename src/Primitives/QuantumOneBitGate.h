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

	static void fillAnyGate(MatrixType& gateMatrix, PsimagLite::String name)
	{
		if (name == "H") {
			fillHadamard(gateMatrix);
			return;
		}

		if (name == "P") {
			fillPhaseOrT(gateMatrix, 0, 1);
			return;
		}

		if (name == "T") {
			RealType oneOverSqrt2 = 1.0/sqrt(2.0);
			fillPhaseOrT(gateMatrix, oneOverSqrt2, oneOverSqrt2);
			return;
		}

		if (name.length() == 2 && name[0] == 'S') {
			SizeType ind = directionCharToInteger(name[1]);
			fillPauli(gateMatrix, ind);
			return;
		}

		if (name.length() >= 2 && name[0] == 'R') {
			SizeType ind = directionCharToInteger(name[1]);
			RealType angle = 0;
			if (name.length() >= 4) {
				if (name[2] != ':')
					err("Expected : in rotation gate name " + name + "\n");
				PsimagLite::String angleStr = name.substr(3, name.length() - 3);
				angle = PsimagLite::atof(angleStr);
			}

			rotation(gateMatrix, ind, angle);
			return;
		}

		err("Gate with name " + name + " not implemented\n");
	}

	static void fillPauli(MatrixType& gateMatrix, SizeType dir)
	{
		gateMatrix.resize(2, 2);

		switch (dir) {
		case 0:
			gateMatrix(0, 1) = gateMatrix(1, 0) = 1;
			break;
		case 1:
			gateMatrix(0, 1) = ComplexOrRealType(0, -1);
			gateMatrix(1, 0) = ComplexOrRealType(0, 1);
			break;
		case 2:
			gateMatrix(0, 0) = 1;
			gateMatrix(1, 1) = -1;
			break;
		default:
			err("Direction can only be 0, 1, or 2\n");
		}
	}


	// ind = 0 means rotation around x
	// ind = 1 means rotation around y
	// ind = 2 means rotation around z
	static void rotation(MatrixType& gateMatrix, SizeType ind, RealType angle)
	{
		const RealType cosine = cos(0.5*angle);
		const RealType sine = sin(0.5*angle);

		gateMatrix.resize(2, 2);
		if (ind == 0) {
			gateMatrix(0, 0) = cosine;
			gateMatrix(0, 1) = ComplexOrRealType(0, -sine);
			gateMatrix(1, 0) = ComplexOrRealType(0, -sine);
			gateMatrix(1, 1) = cosine;
			return;
		} else if (ind == 1) {
			gateMatrix(0, 0) = cosine;
			gateMatrix(0, 1) = -sine;
			gateMatrix(1, 0) = sine;
			gateMatrix(1, 1) = cosine;
			return;
		} else if (ind == 2) {
			gateMatrix(0, 0) = ComplexOrRealType(cosine, -sine);
			gateMatrix(0, 1) = 0;
			gateMatrix(1, 0) = 0;
			gateMatrix(1, 1) = ComplexOrRealType(cosine, sine);
			return;
		}
	}

	// ind = 0 means rotation around x
	// ind = 1 means rotation around y
	// ind = 2 means rotation around z
	static void diffRotation(MatrixType& gateMatrix, SizeType ind, RealType angle)
	{
		const RealType cosine = cos(0.5*angle);
		const RealType sine = sin(0.5*angle);

		gateMatrix.resize(2, 2);
		if (ind == 0) {
			gateMatrix(0, 0) = -sine;
			gateMatrix(0, 1) = ComplexOrRealType(0, -cosine);
			gateMatrix(1, 0) = ComplexOrRealType(0, -cosine);
			gateMatrix(1, 1) = -sine;
			return;
		} else if (ind == 1) {
			gateMatrix(0, 0) = -sine;
			gateMatrix(0, 1) = -cosine;
			gateMatrix(1, 0) = cosine;
			gateMatrix(1, 1) = -sine;
			return;
		} else if (ind == 2) {
			gateMatrix(0, 0) = ComplexOrRealType(-sine, -cosine);
			gateMatrix(0, 1) = 0;
			gateMatrix(1, 0) = 0;
			gateMatrix(1, 1) = ComplexOrRealType(-sine, cosine);
			return;
		}
	}

	static char directionIntegerToChar(SizeType ind)
	{
		if (ind < 3) {
			char c = ind + 120;
			return c;
		}

		throw PsimagLite::RuntimeError("findDirectionOfRotation\n");
	}

	static SizeType directionCharToInteger(char c)
	{
		int val = c - 120;
		if (val >= 0 && val < 3) return val;

		throw PsimagLite::RuntimeError("findDirectionOfRotation\n");
	}

private:

	static void fillHadamard(MatrixType& gateMatrix)
	{
		static const ComplexOrRealType oneOverSqrt2 = 1/sqrt(2.);

		gateMatrix.resize(2, 2);
		gateMatrix(0, 0) = oneOverSqrt2;
		gateMatrix(0, 1) = oneOverSqrt2;
		gateMatrix(1, 0) = oneOverSqrt2;
		gateMatrix(1, 1) = -oneOverSqrt2;
	}

	static void fillPhaseOrT(MatrixType& gateMatrix, RealType a, RealType b)
	{
		gateMatrix.resize(2, 2);
		gateMatrix(0, 0) = 1;
		gateMatrix(1, 1) = ComplexOrRealType(a, b);
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
	      gateMatrix_(gateMatrix)
	{
		code_ += ttos(bitNumber);
		numberOfBits_ = numberOfBits;
	}

	QuantumOneBitGate* clone() const
	{
		return new QuantumOneBitGate(*this);
	}

	virtual PsimagLite::String code() const { return code_; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);

		const ValueType& vv = v[0];
		const int n = vv.size();
		assert(n == (1 << numberOfBits_));  // 2^N

		ValueType w(n);

		for (int i = 0; i < n; ++i) {
			SizeType j = findBasisState(i);
			SizeType bitI = getBitForIndex(i);
			SizeType bitJ = getBitForIndex(j);
			w[i] += gateMatrix_(bitI, bitI)*vv[i];
			w[j] += gateMatrix_(bitI, bitJ)*vv[i];
		}

		return w;
	}

	void setAngle(PsimagLite::String str) const
	{
		PsimagLite::String base;
		PsimagLite::String angleStr;
		extractAngle(base, angleStr, str);

#ifndef NDEBUG
		PsimagLite::String oldAngle;
		PsimagLite::String base2;
		extractAngle(base2, oldAngle, code_);
		assert(base == base2);
#endif

		const bool hA = hasAngles();

		assert(angleStr == "" || hA);

		if (!hA) return;

		code_ = str;

		if (angleStr == "") angleStr = "0";

		RealType angleToUse = std::stod(angleStr);
		assert(code_.size() > 1);
		char cDir = directionOfRotation();
		SizeType ind = OneBitGateLibraryType::directionCharToInteger(cDir);
		OneBitGateLibraryType::rotation(gateMatrix_, ind, angleToUse);
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
		return (code_[0] == 'R' || code_.substr(0, 2) == "_R");
	}

	char directionOfRotation() const
	{
		if (code_[0] == 'R') return code_[1];
		if (code_.substr(0, 2) == "_R") return code_[2];
		throw PsimagLite::RuntimeError("directionOfRotation\n");
	}

	static void extractAngle(PsimagLite::String& base,
	                         PsimagLite::String& angleStr,
	                         PsimagLite::String str)
	{
		typename PsimagLite::String::const_iterator it = std::find(str.begin(),
		                                                           str.end(),
		                                                           ':');

		if (it == str.end()) {
			angleStr = "";
			base = str;
			return;
		}

		base = str.substr(0, it - str.begin());
		angleStr = str.substr(it - str.begin() + 1, str.end() - it - 1);
	}

	static SizeType numberOfBits_;
	mutable PsimagLite::String code_;
	SizeType bitNumber_;
	mutable MatrixType gateMatrix_;
}; // class QuantumOneBitGate

template<typename T>
SizeType QuantumOneBitGate<T>::numberOfBits_ = 0;
}

#endif // QUANTUM_ONE_BIT_GATES_H
