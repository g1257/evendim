#ifndef HAMILTONIANEXAMPLE_H
#define HAMILTONIANEXAMPLE_H
#include "InputNg.h"
#include "InputCheck.h"
#include "PsimagLite.h"
#include "CrsMatrix.h"
#include "HamiltonianFromExpression.h"

namespace Gep {

// H = coupling*\sum_{i} sigma^x_i sigma^x_{i + 1}
template<typename ComplexType>
class HamiltonianExample {

public:

	enum class TypeEnum {FILE, XX, ZZ, EXPRESSION};

	typedef PsimagLite::InputNg<InputCheck> InputNgType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef HamiltonianFromExpression<ComplexType> HamiltonianFromExpressionType;

	HamiltonianExample(typename InputNgType::Readable& io)
	    : hamTipo(TypeEnum::XX), periodic_(false), coupling_(1)
	{
		io.readline(bits_, "NumberOfBits="); // == number of "sites"

		PsimagLite::String ham;
		io.readline(ham, "Hamiltonian=");
		if (ham.substr(0, 5) == "file:") {
			fillFromFile(ham.substr(5, ham.length() - 5));
			const SizeType hilbert = 1<<bits_;
			hamTipo = TypeEnum::FILE;
			if (matrix_.rows() != hilbert)
				err("Matrix rows = " + ttos(matrix_.rows()) + " but " +
				    ttos(hilbert) + " expected.\n");
			return;
		}

		if (ham != "xx" && ham != "yy") {
			HamiltonianFromExpressionType hamExpression(ham, bits_);
			hamExpression.fillMatrix(matrix_);
			cacheVector_.resize(matrix_.rows());
			hamTipo = TypeEnum::EXPRESSION;
			return;
		}

		try {
			int tmp = 0;
			io.readline(tmp, "HamiltonianIsPeriodic=");
			periodic_ = (tmp > 0);
		} catch (std::exception&) {}

		try {
			io.readline(coupling_, "HamiltonianCoupling=");
		} catch (std::exception&) {}

		hamTipo = (ham == "xx") ? TypeEnum::XX : TypeEnum::ZZ;

		if (hamTipo == TypeEnum::XX)
			fillHxx();
		else if (hamTipo != TypeEnum::ZZ)
			err("Hamiltonian=xx or yy but not " + ham + "\n");
	}

	RealType energy(const VectorType& y) const
	{
		if (hamTipo == TypeEnum::ZZ) return energyZZ(y);

		std::fill(cacheVector_.begin(), cacheVector_.end(), 0);
		matrix_.matrixVectorProduct(cacheVector_, y);
		return PsimagLite::real(y*cacheVector_); // does conjugation of first vector
	}

private:

	void fillHxx()
	{
		SizeType hilbertSpace = (1 << bits_);
		matrix_.resize(hilbertSpace, hilbertSpace);
		cacheVector_.resize(hilbertSpace);

		VectorRealType v(hilbertSpace);
		VectorBoolType bcol(hilbertSpace);

		SizeType counter = 0;
		for (SizeType i = 0; i < hilbertSpace; ++i) {
			matrix_.setRow(i, counter);
			std::fill(v.begin(), v.end(), 0);
			std::fill(bcol.begin(), bcol.end(), false);
			for (SizeType site = 0; site < bits_; ++site) {
				// Flip bit at site
				SizeType maskSite = (1 << site);
				SizeType j = i ^ maskSite;
				SizeType site2 = site + 1;
				if (site2 >= bits_ && !periodic_) continue;
				assert(site2 <= bits_);
				if (site2 == bits_) site2 = 0;
				SizeType maskSite2 = (1 << site2);
				j ^= maskSite2;
				v[j] += coupling_;
				bcol[j] = true;
			}

			counter += fillThisRow(v, bcol);
		}

		matrix_.setRow(hilbertSpace, counter);
		matrix_.checkValidity();

		VectorRealType eigs(hilbertSpace);
		PsimagLite::Matrix<ComplexType> a = matrix_.toDense();
		diag(a, eigs, 'V');
		std::cout<<"Ground State Energy="<<eigs[0]<<"\n";
	}

	SizeType fillThisRow(VectorRealType& v, VectorBoolType& bcols)
	{
		const SizeType hilbertSpace = v.size();
		assert(hilbertSpace == bcols.size());
		SizeType counter = 0;
		for (SizeType i = 0; i < hilbertSpace; ++i) {
			if (!bcols[i]) continue;
			matrix_.pushCol(i);
			matrix_.pushValue(v[i]);
			++counter;
			bcols[i] = false;
			v[i] = 0;
		}

		return counter;
	}

	RealType energyZZ(const VectorType& v) const
	{
		const SizeType hilbertSpace = v.size();
		RealType e = 0;
		for (SizeType i = 0; i < hilbertSpace; ++i) {
			for (SizeType site = 0; site < bits_; ++site) {
				SizeType maskSite = (1 << site);
				SizeType j = i & maskSite;
				SizeType site2 = site + 1;
				if (site2 >= bits_ && !periodic_) continue;
				assert(site2 <= bits_);
				if (site2 == bits_) site2 = 0;
				SizeType maskSite2 = (1 << site2);
				SizeType k = i & maskSite2;
				SizeType jj = (j > 0);
				SizeType kk = (k > 0);
				RealType tmp = PsimagLite::real(PsimagLite::conj(v[i])*v[i]);
				RealType value = (jj == kk) ? tmp : -tmp;
				e += value;
			}
		}

		return e*coupling_;
	}

	void fillFromFile(PsimagLite::String filename)
	{
		SizeType hilbertSpace = (1 << bits_);
		cacheVector_.resize(hilbertSpace);
		std::ifstream fin(filename);
		if (!fin || !fin.good())
			err("Could not open file " + filename + "\n");

		SizeType rows = 0;
		fin>>rows;

		SizeType cols =0;
		fin>>cols;
		PsimagLite::Matrix<ComplexType> mat(rows, cols);
		for (SizeType i = 0; i < rows; ++i)
			for (SizeType j = 0; j < cols; ++j)
				fin >> mat(i, j);
		fullMatrixToCrsMatrix(matrix_, mat);
		VectorRealType eigs(hilbertSpace);
		diag(mat, eigs, 'V');
		std::cout<<"Ground State Energy="<<eigs[0]<<"\n";
	}

	TypeEnum hamTipo;
	bool periodic_;
	SizeType bits_;
	RealType coupling_;
	SparseMatrixType matrix_;
	mutable VectorType cacheVector_;
};
}
#endif // HAMILTONIANEXAMPLE_H
