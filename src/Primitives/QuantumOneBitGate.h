#ifndef QUANTUM_ONE_BIT_GATES_H
#define QUANTUM_ONE_BIT_GATES_H
#include "AST/Node.h"
#include "CustomQuantumGates.hh"
#include "Matrix.h"
#include "QuasiVector.hh"

namespace Gep {

template <typename ComplexOrRealType>
class OneBitGateLibrary {

public:

	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef CustomQuantumGates<ComplexOrRealType> CustomQuantumGatesType;

	static void fillAnyGate(MatrixType& gateMatrix,
	                        PsimagLite::String name,
	                        const CustomQuantumGatesType* customQuantumGates = nullptr)
	{
		if (name.substr(0, 2) == "CG" || name.substr(0, 2) == "PG") {
			assert(customQuantumGates);
			customQuantumGates->evaluate(gateMatrix, name);
		}
		else {
			fillAnyGateBasic(gateMatrix, name);
		}
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
		const RealType cosine = cos(0.5 * angle);
		const RealType sine = sin(0.5 * angle);

		gateMatrix.resize(2, 2);
		if (ind == 0) {
			gateMatrix(0, 0) = cosine;
			gateMatrix(0, 1) = ComplexOrRealType(0, -sine);
			gateMatrix(1, 0) = ComplexOrRealType(0, -sine);
			gateMatrix(1, 1) = cosine;
			return;
		}
		else if (ind == 1) {
			gateMatrix(0, 0) = cosine;
			gateMatrix(0, 1) = -sine;
			gateMatrix(1, 0) = sine;
			gateMatrix(1, 1) = cosine;
			return;
		}
		else if (ind == 2) {
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
		const RealType cosine = 0.5 * cos(0.5 * angle);
		const RealType sine = 0.5 * sin(0.5 * angle);

		gateMatrix.resize(2, 2);
		if (ind == 0) {
			gateMatrix(0, 0) = -sine;
			gateMatrix(0, 1) = ComplexOrRealType(0, -cosine);
			gateMatrix(1, 0) = ComplexOrRealType(0, -cosine);
			gateMatrix(1, 1) = -sine;
			return;
		}
		else if (ind == 1) {
			gateMatrix(0, 0) = -sine;
			gateMatrix(0, 1) = -cosine;
			gateMatrix(1, 0) = cosine;
			gateMatrix(1, 1) = -sine;
			return;
		}
		else if (ind == 2) {
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
		if (val >= 0 && val < 3)
			return val;

		throw PsimagLite::RuntimeError("findDirectionOfRotation\n");
	}

private:

	static void fillAnyGateBasic(MatrixType& gateMatrix, PsimagLite::String name)
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
			RealType oneOverSqrt2 = 1.0 / sqrt(2.0);
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

		if (name.length() > 2 && name.substr(0, 2) == "_R") {
			SizeType ind = directionCharToInteger(name[2]);
			RealType angle = 0;
			if (name.length() > 4) {
				if (name[3] != ':')
					err("Expected : in rotation gate name " + name + "\n");
				PsimagLite::String angleStr = name.substr(4, name.length() - 4);
				angle = PsimagLite::atof(angleStr);
			}

			diffRotation(gateMatrix, ind, angle);
			return;
		}

		err("Gate with name " + name + " not implemented\n");
	}

	static void fillHadamard(MatrixType& gateMatrix)
	{
		static const ComplexOrRealType oneOverSqrt2 = 1 / sqrt(2.);

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

template <typename VectorValueType>
class QuantumOneBitGate : public PsimagLite::Node<VectorValueType,
                                                  typename PsimagLite::Real<typename VectorValueType::value_type::value_type>::Type> {

public:

	typedef typename VectorValueType::value_type ValueType;
	typedef typename ValueType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef OneBitGateLibrary<ComplexOrRealType> OneBitGateLibraryType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename OneBitGateLibraryType::CustomQuantumGatesType CustomQuantumGatesType;

	QuantumOneBitGate(PsimagLite::String cr,
	                  SizeType bitNumber,
	                  SizeType numberOfBits,
	                  const MatrixType& gateMatrix)
	    : code_(cr)
	    , bitNumber_(bitNumber)
	    , gateMatrix_(gateMatrix)
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
		assert(vv.size() == static_cast<SizeType>(1 << numberOfBits_)); // 2^N

		return oneBitGate(vv, bitNumber_, gateMatrix_);
	}

	void setAngle(PsimagLite::String str) const
	{
		PsimagLite::String str2 = deleteSite(str);
		OneBitGateLibraryType::fillAnyGate(gateMatrix_, str2, customOneBitGate_);
	}

	static void setCustom(const CustomQuantumGatesType& customQuantumGates)
	{
		customOneBitGate_ = &customQuantumGates;
	}

private:

	static PsimagLite::String deleteSite(PsimagLite::String str)
	{
		long unsigned int ind = str.find(":");
		if (ind == PsimagLite::String::npos) {
			return deleteSiteNoColon(str);
		}

		SizeType l = str.length();
		PsimagLite::String str2 = str.substr(0, ind);
		PsimagLite::String str3 = deleteSiteNoColon(str2);
		PsimagLite::String str4 = str.substr(ind, l - ind);

		return str3 + str4;
	}

	static PsimagLite::String deleteSiteNoColon(PsimagLite::String str)
	{
		// delete site
		const SizeType l = str.length();
		SizeType counter = 0;
		for (SizeType i = 0; i < l; ++i) {
			SizeType j = l - i - 1;
			unsigned int num = str[j];
			if (num > 57 || num < 48)
				break;
			++counter;
		}

		SizeType ll = l - counter;
		assert(ll > 0);
		return str.substr(0, ll);
	}

	bool hasAngles() const
	{
		assert(code_.size() > 0);
		return (code_[0] == 'R' || code_.substr(0, 2) == "_R" || code_.substr(0, 2) == "PG");
	}

	char directionOfRotation() const
	{
		if (code_[0] == 'R')
			return code_[1];
		if (code_.substr(0, 2) == "_R")
			return code_[2];
		if (code_.substr(0, 2) == "PG")
			return 'x';
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
	static CustomQuantumGatesType const* customOneBitGate_;
	mutable PsimagLite::String code_;
	SizeType bitNumber_;
	mutable MatrixType gateMatrix_;
}; // class QuantumOneBitGate

template <typename T>
SizeType QuantumOneBitGate<T>::numberOfBits_ = 0;

template <typename T>
typename QuantumOneBitGate<T>::CustomQuantumGatesType const* QuantumOneBitGate<T>::customOneBitGate_ = nullptr;
}

#endif // QUANTUM_ONE_BIT_GATES_H
