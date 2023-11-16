#ifndef HAMILTONIANFROMEXPRESSION_H
#define HAMILTONIANFROMEXPRESSION_H
#include "AuxForHamSpec.h"
#include "CanonicalExpression.h"
#include "CrsMatrix.h"
#include "HamiltonianSpec.h"
#include "Vector.h"

namespace Gep {

template <typename ComplexType>
class HamiltonianFromExpression {

public:

	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef HamiltonianSpec<SparseMatrixType> HamiltonianSpecType;
	typedef PsimagLite::CanonicalExpression<HamiltonianSpecType> CanonicalExpressionType;
	typedef AuxForHamSpec AuxForHamSpecType;

	// 0.42*P3;H2+-4.2*Rx0:1.57*Sy1
	HamiltonianFromExpression(PsimagLite::String expression,
	                          SizeType numberOfBits)
	    : hamString_(expression)
	    , bits_(numberOfBits)
	{
		fillMatrix();
	}

	const SparseMatrixType& getMatrix() const { return matrix_; }

	// just for checking
	static void solveIt(const SparseMatrixType& matrix)
	{
		typedef typename PsimagLite::Real<ComplexType>::Type RealType;
		PsimagLite::Matrix<ComplexType> dense;
		crsMatrixToFullMatrix(dense, matrix);
		typename PsimagLite::Vector<RealType>::Type eigs(dense.rows());
		diag(dense, eigs, 'V');
		std::cout << "gs energy=" << eigs[0] << "\n";
	}

private:

	void fillMatrix()
	{
		HamiltonianSpecType hamSpec(bits_);
		CanonicalExpressionType canonicalExpression(hamSpec);
		hamString_ = killSpaces(hamString_);
		AuxForHamSpecType aux(bits_);
		typename HamiltonianSpecType::ResultType quasiMatrix(bits_);
		typename HamiltonianSpecType::ResultType emptyMatrix(bits_);
		canonicalExpression(quasiMatrix, hamString_, emptyMatrix, aux);
		matrix_ = quasiMatrix.getCRS();
	}

	static PsimagLite::String killSpaces(PsimagLite::String str)
	{
		PsimagLite::String buffer;
		const SizeType n = str.length();
		for (SizeType i = 0; i < n; ++i)
			if (str[i] != ' ')
				buffer += str[i];
		return buffer;
	}

	PsimagLite::String hamString_;
	SizeType bits_;
	SparseMatrixType matrix_;
};
}
#endif // HAMILTONIANFROMEXPRESSION_H
