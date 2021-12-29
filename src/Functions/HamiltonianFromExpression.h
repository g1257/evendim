#ifndef HAMILTONIANFROMEXPRESSION_H
#define HAMILTONIANFROMEXPRESSION_H
#include "Vector.h"
#include "CanonicalExpression.h"
#include "HamiltonianSpec.h"
#include "CrsMatrix.h"
#include "AuxForHamSpec.h"

namespace Gep {

template<typename ComplexType>
class HamiltonianFromExpression {

public:

	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef HamiltonianSpec<SparseMatrixType> HamiltonianSpecType;
	typedef PsimagLite::CanonicalExpression<HamiltonianSpecType> CanonicalExpressionType;
	typedef AuxForHamSpec AuxForHamSpecType;

	// 0.42*P3;H2+-4.2*Rx0:1.57*Sy1
	HamiltonianFromExpression(PsimagLite::String expression,
	                          SizeType numberOfBits)
	    : hamString_(expression), bits_(numberOfBits)
	{}

	void fillMatrix(SparseMatrixType& matrix)
	{
		HamiltonianSpecType hamSpec(bits_);
		CanonicalExpressionType canonicalExpression(hamSpec);
		hamString_ = killSpaces(hamString_);
		AuxForHamSpecType aux(bits_);
		typename HamiltonianSpecType::ResultType quasiMatrix(bits_);
		typename HamiltonianSpecType::ResultType emptyMatrix(bits_);
		canonicalExpression(quasiMatrix, hamString_, emptyMatrix, aux);
		matrix = quasiMatrix.getCRS();
	}

private:

	static PsimagLite::String killSpaces(PsimagLite::String str)
	{
		PsimagLite::String buffer;
		const SizeType n = str.length();
		for (SizeType i = 0; i < n; ++i)
			if (str[i] != ' ') buffer += str[i];
		return buffer;
	}

	PsimagLite::String hamString_;
	SizeType bits_;
};
}
#endif // HAMILTONIANFROMEXPRESSION_H
