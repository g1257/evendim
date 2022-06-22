#ifndef CUSTOMQUANTUMGATES_HH
#define CUSTOMQUANTUMGATES_HH
#include "Vector.h"
#include "InputCheck.h"
#include "InputNg.h"
#include <unordered_map>

namespace Gep {

template<typename ComplexType>
class CustomQuantumGates {

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::InputNg<InputCheck>::Readable InputNgReadableType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;

public:

	void push(const std::string& name, const PsimagLite::Matrix<PsimagLite::String>& symbolicMatrix)
	{
		if (indices_.count(name) > 0) {
			err("Matrix named " + name + " already defined\n");
		}

		indices_[name] = symbolicMatrices_.size();
		symbolicMatrices_.push_back(symbolicMatrix);
	}

	void push(const std::string& name, const MatrixType& matrix)
	{
		if (indices_.count(name) > 0) {
			err("Matrix named " + name + " already defined\n");
		}

		indices_[name] = numericMatrices_.size();
		numericMatrices_.push_back(matrix);
	}

	void evaluate(MatrixType& matrix, const std::string& name) const
	{
		if (name.substr(0, 2) == "CG") return;
		assert(name.substr(0, 2) == "PG");

		SizeType index = indices_.at(name);
		assert(index < symbolicMatrices_.size());
		const PsimagLite::Matrix<PsimagLite::String>& symbolicMatrix = symbolicMatrices_[index];
		evaluateSymbolicMatrix(matrix, symbolicMatrix);
	}

private:

	static void evaluateSymbolicMatrix(MatrixType& matrix,
	                                   const PsimagLite::Matrix<PsimagLite::String>& symbolicMatrix)
	{
		SizeType rows = symbolicMatrix.rows();
		SizeType cols = symbolicMatrix.cols();
		matrix.resize(rows, cols);
		for (SizeType i = 0; i < rows; ++i) {
			for (SizeType j = 0; j < cols; ++j) {
				matrix(i, j) = PsimagLite::atof(symbolicMatrix(i, j));
			}
		}
	}


	PsimagLite::Vector<PsimagLite::Matrix<PsimagLite::String> >::Type symbolicMatrices_;
	typename PsimagLite::Vector<MatrixType>::Type numericMatrices_;
	std::unordered_map<PsimagLite::String, SizeType> indices_;
};
}
#endif // CUSTOMQUANTUMGATES_HH
