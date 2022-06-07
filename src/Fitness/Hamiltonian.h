#ifndef EVENDIM_HAMILTONIAN_H
#define EVENDIM_HAMILTONIAN_H
#include "InputNg.h"
#include "InputCheck.h"
#include "PsimagLite.h"
#include "CrsMatrix.h"
#include "HamiltonianFromExpression.h"
#include "IsingGraph.hh"

namespace Gep {

// H = coupling*\sum_{i} sigma^x_i sigma^x_{i + 1}
template<typename ComplexType>
class Hamiltonian {

public:

	enum class TypeEnum {FILE, XX, EXPRESSION, ISING_GRAPH};

	typedef PsimagLite::InputNg<InputCheck> InputNgType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef typename PsimagLite::Vector<VectorType>::Type VectorVectorType;
	typedef typename PsimagLite::Vector<int>::Type VectorIntType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef HamiltonianFromExpression<ComplexType> HamiltonianFromExpressionType;
	typedef IsingGraph<ComplexType> IsingGraphType;

	Hamiltonian(typename InputNgType::Readable& io, SizeType numberOfThreads)
	    : hamTipo(TypeEnum::XX),
	      bits_(0),
	      periodic_(false),
	      isingGraph_(nullptr),
	      needsTransformAndTruncate_(false),
	      cacheVector_(numberOfThreads)
	{
		io.readline(bits_, "NumberOfBits="); // == number of "sites"

		PsimagLite::String ham;
		io.readline(ham, "Hamiltonian=");
		if (ham.substr(0, 5) == "file:") {

			fillFromFile(ham.substr(5, ham.length() - 5), io);

			hamTipo = TypeEnum::FILE;

			SizeType hilbert = (1 << bits_);

			if (matrix_.rows() != hilbert && !needsTransformAndTruncate_)
				err("Matrix rows = " + ttos(matrix_.rows()) + " but " +
				    ttos(hilbert) + " expected or a Basis= line needed\n");

			assert(isHermitian(matrix_, true));

			allocateCacheVector(hilbert);

			return;
		}

		RealType coupling = 1;
		try {
			io.readline(coupling, "HamiltonianCoupling=");
		} catch (std::exception&) {}

		bool hasPeriodic = false;
		try {
			int tmp = 0;
			io.readline(tmp, "HamiltonianIsPeriodic=");
			periodic_ = (tmp > 0);
			hasPeriodic = true;
		} catch (std::exception&) {}

		if (ham == "IsingGraph" || ham == "zz") {
			PsimagLite::String graphFile = "zz";
			if (ham == "IsingGraph") {

				if (hasPeriodic)
					err("IsingGraph does not support HamiltonianIsPeriodic= line\n");

				io.readline(graphFile, "GraphFile=");
				graphFile = "file:" + graphFile;
			}

			isingGraph_ = new IsingGraphType(bits_, coupling, periodic_, graphFile);
			hamTipo = TypeEnum::ISING_GRAPH;
			return;
		}

		if (ham != "xx") {
			HamiltonianFromExpressionType hamExpression(ham, bits_);
			hamExpression.fillMatrix(matrix_);
			assert(cacheVector_.size() > 0);
			allocateCacheVector(matrix_.rows());
			hamTipo = TypeEnum::EXPRESSION;
			return;
		}

		assert(ham == "xx");
		hamTipo = TypeEnum::XX;
		fillHxx(coupling);
	}

	RealType energy(const VectorType& y, SizeType threadNum) const
	{
		switch (hamTipo) {
		case  TypeEnum::ISING_GRAPH: {
			assert(isingGraph_);
			return isingGraph_->energyZZ(y);
			break;
		}

		default: {
			assert(cacheVector_.size() > threadNum);
			assert(cacheVector_[threadNum].size() == matrix_.rows());
			std::fill(cacheVector_[threadNum].begin(), cacheVector_[threadNum].end(), 0);

			matrix_.matrixVectorProduct(cacheVector_[threadNum], y);

			return PsimagLite::real(y*cacheVector_[threadNum]); // does conjugation of first vector
			break;
		}

		}
	}

	template<typename SomeChromosomeType>
	PsimagLite::String info(const SomeChromosomeType& chromosome) const
	{
		if (hamTipo != TypeEnum::ISING_GRAPH) return "";
		assert(isingGraph_);
		return isingGraph_->info(chromosome.exec(0));
	}

private:

	void fillHxx(RealType coupling)
	{
		SizeType hilbertSpace = (1 << bits_);
		matrix_.resize(hilbertSpace, hilbertSpace);
		allocateCacheVector(hilbertSpace);

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
				v[j] += coupling;
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

	void fillFromFile(PsimagLite::String filename, typename InputNgType::Readable& io)
	{
		std::ifstream fin(filename);
		if (!fin || !fin.good())
			err("Could not open file " + filename + "\n");

		SizeType rows = 0;
		fin>>rows;

		SizeType cols =0;
		fin>>cols;
		if (rows != cols)
			err("Hamiltonian must have rows==cols\n");

		PsimagLite::Matrix<ComplexType> mat(rows, cols);
		for (SizeType i = 0; i < rows; ++i)
			for (SizeType j = 0; j < cols; ++j)
				fin >> mat(i, j);

		try {
			io.read(basis_, "Basis=");
			needsTransformAndTruncate_ = true;
			transformAndTruncate(mat);
		} catch (std::exception&) {
			fullMatrixToCrsMatrix(matrix_, mat);
			printGs(mat);
		}
	}

	static void printGs(PsimagLite::Matrix<ComplexType>& mat)
	{
		assert(mat.rows() == mat.cols());
		VectorRealType eigs(mat.rows());
		diag(mat, eigs, 'V');
		std::cout<<"Ground State Energy="<<eigs[0]<<"\n";
		std::cout<<mat;
	}

	void allocateCacheVector(SizeType hilbertSpace)
	{
		for (SizeType thread = 0; thread < cacheVector_.size(); ++thread)
			cacheVector_[thread].resize(hilbertSpace);
	}

	void transformAndTruncate(PsimagLite::Matrix<ComplexType>& mat)
	{
		SizeType n = mat.rows();
		assert(n == mat.cols());
		PsimagLite::Matrix<ComplexType> dense(n, n);

		for (SizeType i = 0; i < n; ++i) {
			int ii = basis_[i];
			if (ii < 0)
				continue;
			for (SizeType j = 0; j < n; ++j) {
				int jj = basis_[j];
				if (j < 0)
					continue;

				dense(i, j) = mat(ii, jj);
			}
		}

		fullMatrixToCrsMatrix(matrix_, dense);
		printGs(dense);
	}

	TypeEnum hamTipo;
	SizeType bits_;
	bool periodic_;
	IsingGraphType* isingGraph_;
	SparseMatrixType matrix_;
	VectorIntType basis_;
	bool needsTransformAndTruncate_;
	mutable VectorVectorType cacheVector_;
};
}
#endif // EVENDIM_HAMILTONIAN_H
