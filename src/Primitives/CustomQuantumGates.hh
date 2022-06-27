#ifndef CUSTOMQUANTUMGATES_HH
#define CUSTOMQUANTUMGATES_HH
#include "Vector.h"
#include "InputCheck.h"
#include "InputNg.h"
#include "Sort.h"
#include "ExpressionCalculator.h"
#include "PsimagLite.h"
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

	void evaluate(MatrixType& matrix, const std::string& name2) const
	{
		if (name2.substr(0, 2) == "CG") return;
		assert(name2.substr(0, 2) == "PG");

		VectorStringType tokens;
		PsimagLite::split(tokens, name2, ":");
		double angle = 0;
		std::string name = name2;
		assert(tokens.size() == 1 || tokens.size() == 2);
		if (tokens.size() > 1) {
			name = tokens[0];
			angle = PsimagLite::atof(tokens[1]);
		}

		SizeType index = indices_.at(name);
		assert(index < symbolicMatrices_.size());
		const PsimagLite::Matrix<PsimagLite::String>& symbolicMatrix = symbolicMatrices_[index];
		evaluateSymbolicMatrix(matrix, symbolicMatrix, name, angle);
	}

private:

	static void evaluateSymbolicMatrix(MatrixType& matrix,
	                                   const PsimagLite::Matrix<PsimagLite::String>& symbolicMatrix,
	                                   const std::string& name,
	                                   double angle)
	{
		SizeType nparams = getNumberOfParams(symbolicMatrix, name);
		if (nparams != 1)
			err("Parametric gates can only have one parameter for now\n");

		std::vector<double> params(nparams, angle);
		SizeType rows = symbolicMatrix.rows();
		SizeType cols = symbolicMatrix.cols();
		matrix.resize(rows, cols);
		for (SizeType i = 0; i < rows; ++i) {
			for (SizeType j = 0; j < cols; ++j) {
				matrix(i, j) = maybeReplaceParam(symbolicMatrix(i, j), params);
			}
		}
	}

	static SizeType getNumberOfParams(const PsimagLite::Matrix<PsimagLite::String>& symbolicMatrix,
	                                  const std::string& name)
	{
		SizeType rows = symbolicMatrix.rows();
		SizeType cols = symbolicMatrix.cols();
		std::vector<SizeType> seen;
		for (SizeType i = 0; i < rows; ++i) {
			for (SizeType j = 0; j < cols; ++j) {
				pushParamsIfAny(seen, symbolicMatrix(i, j));
			}
		}

		SizeType maxParams = seen.size();
		if (maxParams == 0)
			err("Gate named " + name + " contains no parameters. Use CG instead.\n");

		PsimagLite::Sort<std::vector<SizeType> > sort;
		std::vector<SizeType> iperm(maxParams);
		sort.sort(seen, iperm);
		for (SizeType i = 0; i < maxParams; ++i) {
			if (seen[i] != i)
				err("Gate named " + name + " has parameter gap\n");
		}

		return maxParams;
	}

	// Up to 10 parameters from 0 to 9
	static void pushParamsIfAny(std::vector<SizeType>& seen, const std::string& value)
	{
		SizeType l = value.length();
		for (SizeType i = 0; i < l; ++i) {
			if (value[i] == 'p') {
				if (i + 1 == l)
					err("Symbolic value " + value + " syntax error. Nothing follows p\n");
				unsigned char c = value[i + 1];
				if ((c < 48) || (c > 57))
					err("Symbolic value " + value + " syntax error. p must be followed by a digit.\n");
				SizeType x = c - 48;
				if (std::find(seen.begin(), seen.end(), x) == seen.end())
					seen.push_back(x);
			}
		}
	}

	static std::complex<double> maybeReplaceParam(const std::string& value, const std::vector<double>& params)
	{
		typedef PsimagLite::ExpressionCalculator<std::complex<double> > ExpressionCalculatorType;
		typedef PsimagLite::PrepassData<double> PrepassDataType;

		PsimagLite::String paramsString = "p0";
		for (SizeType i = 1; i < params.size(); ++i) {
			paramsString += ",p" + ttos(i);
		}

		VectorStringType expression;

		PsimagLite::split(expression, value, ":");
		PrepassDataType pd(paramsString, params);
		PsimagLite::ExpressionPrepass<PrepassDataType>::prepass(expression, pd);
		ExpressionCalculatorType ec(expression);
		return ec();
	}

	PsimagLite::Vector<PsimagLite::Matrix<PsimagLite::String> >::Type symbolicMatrices_;
	typename PsimagLite::Vector<MatrixType>::Type numericMatrices_;
	std::unordered_map<PsimagLite::String, SizeType> indices_;
};
}
#endif // CUSTOMQUANTUMGATES_HH
