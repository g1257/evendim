#ifndef EVENDIM_HAMILTONIAN_CPU_H
#define EVENDIM_HAMILTONIAN_CPU_H
#include "../Engine/InputCheck.h"
#include "../Primitives/QuasiVector.hh"
#include "CrsMatrix.h"
#include "HamiltonianFromExpression.h"
#include "InputNg.h"
#include "IsingGraph.hh"
#include "PsimagLite.h"
#include "SchwingerModel.hh"

namespace Gep {

// H = coupling*\sum_{i} sigma^x_i sigma^x_{i + 1}
template <typename ComplexType>
class Hamiltonian {

public:

	enum class TypeEnum { FILE,
		              XX,
		              EXPRESSION,
		              ISING_GRAPH };

	typedef PsimagLite::InputNg<InputCheck> InputNgType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef QuasiVector<ComplexType> QuasiVectorType;
	typedef typename PsimagLite::Vector<QuasiVectorType>::Type VectorQuasiVectorType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef HamiltonianFromExpression<ComplexType> HamiltonianFromExpressionType;
	typedef IsingGraph<ComplexType> IsingGraphType;

	Hamiltonian(typename InputNgType::Readable& io, SizeType numberOfThreads)
	    : hamTipo(TypeEnum::XX)
	    , bits_(0)
	    , periodic_(false)
	    , isingGraph_(nullptr)
	    , needsTransformAndTruncate_(false)
	{
		io.readline(bits_, "NumberOfBits="); // == number of "sites"

		PsimagLite::String ham;
		io.readline(ham, "Hamiltonian=");
		if (ham.substr(0, 5) == "file:") {

			fillFromFile(ham.substr(5, ham.length() - 5), io);

			hamTipo = TypeEnum::FILE;

			SizeType hilbert = (1 << bits_);

			if (matrix_.rows() != hilbert && !needsTransformAndTruncate_)
				err("Matrix rows = " + ttos(matrix_.rows()) + " but " + ttos(hilbert) + " expected or a Basis= line needed\n");

			assert(isHermitian(matrix_, true));

			return;
		}

		RealType coupling = 1;
		try {
			io.readline(coupling, "HamiltonianCoupling=");
		}
		catch (std::exception&) {
		}

		bool hasPeriodic = false;
		try {
			int tmp = 0;
			io.readline(tmp, "HamiltonianIsPeriodic=");
			periodic_ = (tmp > 0);
			hasPeriodic = true;
		}
		catch (std::exception&) {
		}

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
			isingGraph_->solve();
			return;
		}

		if (ham == "xx") {
			hamTipo = TypeEnum::XX;
			fillHxx(coupling);
			return; // EARLY EXIT HERE
		}

		if (ham == "zxz") {
			PsimagLite::String strZxZ = createNnn("Sz", "Sx", "Sz", bits_);
			HamiltonianFromExpressionType hamZxZ(strZxZ, bits_);
			const SparseMatrixType& matrixZxZ = hamZxZ.getMatrix();

			PsimagLite::String strX = createLocal("Sx", bits_);
			HamiltonianFromExpressionType hamX(strX, bits_);
			const SparseMatrixType& matrixX = hamX.getMatrix();

			PsimagLite::String strXx = createNn("Sx", "Sx", bits_);
			HamiltonianFromExpressionType hamXx(strXx, bits_);
			const SparseMatrixType& matrixXx = hamXx.getMatrix();

			RealType hJ = 0;
			io.readline(hJ, "HamiltonianJ=");
			matrix_ = hJ * matrixZxZ;

			RealType h1 = 0;
			io.readline(h1, "Hamiltonianh1=");
			addMatrixWithWeight(matrix_, h1, matrixX);

			RealType h2 = 0;
			io.readline(h2, "Hamiltonianh2=");
			addMatrixWithWeight(matrix_, h2, matrixXx);
		}
		else if (ham == "schwinger") {
			RealType m_param = 0;
			io.readline(m_param, "Hamiltonianm=");

			RealType g_param = 0;
			io.readline(g_param, "Hamiltoniang=");

			SchwingerModel<ComplexType> schwinger_model(bits_, periodic_, m_param, g_param);
			matrix_ = schwinger_model.matrix();

			PsimagLite::Matrix<ComplexType> a = matrix_.toDense();
			printGs(a);
		}
		else {
			std::cerr << "Asumming Hamiltonian Expression\n";
			HamiltonianFromExpressionType hamExpression(ham, bits_);
			matrix_ = hamExpression.getMatrix();
		}

		hamTipo = TypeEnum::EXPRESSION;
		HamiltonianFromExpressionType::solveIt(matrix_);
	}

	RealType energy(const QuasiVectorType& y, SizeType /* threadNum */) const
	{
		switch (hamTipo) {
		case TypeEnum::ISING_GRAPH: {
			assert(isingGraph_);
			return isingGraph_->energyZZ(y.toVector());
			break;
		}

		default: {
			return tensorEnergy(y, matrix_, y);
			break;
		}
		}
	}

	SizeType numberOfSites() const { return bits_; }

private:

	void fillHxx(RealType coupling)
	{
		SizeType hilbertSpace = (1 << bits_);
		matrix_.resize(hilbertSpace, hilbertSpace);

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
				if (site2 >= bits_ && !periodic_)
					continue;
				assert(site2 <= bits_);
				if (site2 == bits_)
					site2 = 0;
				SizeType maskSite2 = (1 << site2);
				j ^= maskSite2;
				v[j] += coupling;
				bcol[j] = true;
			}

			counter += fillThisRow(v, bcol);
		}

		matrix_.setRow(hilbertSpace, counter);
		matrix_.checkValidity();

		PsimagLite::Matrix<ComplexType> a = matrix_.toDense();
		printGs(a);
	}

	SizeType fillThisRow(VectorRealType& v, VectorBoolType& bcols)
	{
		const SizeType hilbertSpace = v.size();
		assert(hilbertSpace == bcols.size());
		SizeType counter = 0;
		for (SizeType i = 0; i < hilbertSpace; ++i) {
			if (!bcols[i])
				continue;
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
		fin >> rows;

		SizeType cols = 0;
		fin >> cols;
		if (rows != cols)
			err("Hamiltonian must have rows==cols\n");

		PsimagLite::Matrix<ComplexType> mat(rows, cols);
		for (SizeType i = 0; i < rows; ++i)
			for (SizeType j = 0; j < cols; ++j)
				fin >> mat(i, j);

		bool hasScale = false;
		VectorType scale;
		try {
			io.read(scale, "ScaleHamiltonian");
			hasScale = true;
		}
		catch (std::exception&) {
		}

		scaleHamiltonian(mat, scale, hasScale);

		try {
			io.read(basis_, "Basis");
			needsTransformAndTruncate_ = true;
			transformAndTruncate(mat);
			std::cerr << "Has basis of size " << basis_.size() << "\n";
		}
		catch (std::exception&) {
			fullMatrixToCrsMatrix(matrix_, mat);
			printGs(mat);
		}
	}

	static void scaleHamiltonian(PsimagLite::Matrix<ComplexType>& mat,
	                             const VectorType& scale,
	                             bool hasScale)
	{
		if (!hasScale)
			return;

		if (scale.size() != 2)
			err("Expecting ScaleHamiltonian a vector of two entries\n");

		for (SizeType i = 0; i < mat.rows(); ++i) {
			for (SizeType j = 0; j < mat.cols(); ++j) {
				ComplexType val = scale[0] * mat(i, j);
				if (i == j)
					val += scale[1];
				mat(i, j) = val;
			}
		}
	}

	static void printGs(PsimagLite::Matrix<ComplexType>& mat)
	{
		assert(mat.rows() == mat.cols());
		VectorRealType eigs(mat.rows());
		diag(mat, eigs, 'V');
		std::cout << "Ground State Energy=" << eigs[0] << "\n";

		if (eigs.size() > 100)
			return; // <<--- EARLY EXIT HERE
			        //
		std::cout << "Eigenvector------------\n";
		ComplexType sum = 0;
		for (SizeType i = 0; i < mat.rows(); ++i) {
			ComplexType val = mat(i, 0);
			sum += val * PsimagLite::conj(val);
			if (std::norm(val) < 1e-8)
				continue;

			std::cout << i << " " << mat(i, 0) << "\n";
		}

		std::cout << "-------- End eigenvector=" << sum << "\n\n";
	}

	void transformAndTruncate(PsimagLite::Matrix<ComplexType>& mat)
	{
		SizeType n = mat.rows();
		assert(n == mat.cols());
		SizeType hilbert = (1 << bits_);
		PsimagLite::Matrix<ComplexType> dense(hilbert, hilbert);
		assert(basis_.size() == n);

		for (SizeType i = 0; i < n; ++i) {
			SizeType ii = basis_[i];
			for (SizeType j = 0; j < n; ++j) {
				SizeType jj = basis_[j];

				dense(ii, jj) = mat(i, j);
			}
		}
		// dense(8+4+2+1=15,) = mat(0, 0);

		fullMatrixToCrsMatrix(matrix_, dense);
		printGs(dense);
	}

	static PsimagLite::String createNnn(const PsimagLite::String& A,
	                                    const PsimagLite::String& B,
	                                    const PsimagLite::String& C,
	                                    SizeType n)
	{
		if (n < 3)
			err("createNnn needs at least three sites\n");
		SizeType nMinusTwo = n - 2;
		assert(nMinusTwo < n);
		PsimagLite::String buffer;
		for (SizeType i = 0; i < nMinusTwo; ++i) {
			SizeType j = i + 1;
			SizeType k = i + 2;
			if (i > 0)
				buffer += "+";
			buffer += A + ttos(i) + "*" + B + ttos(j) + "*" + C + ttos(k);
		}

		return buffer;
	}

	static PsimagLite::String createNn(const PsimagLite::String& A,
	                                   const PsimagLite::String& B,
	                                   SizeType n)
	{
		if (n < 2)
			err("createNn needs at least two sites\n");
		SizeType nMinusOne = n - 1;
		assert(nMinusOne < n);
		PsimagLite::String buffer;
		for (SizeType i = 0; i < nMinusOne; ++i) {
			SizeType j = i + 1;
			if (i > 0)
				buffer += "+";
			buffer += A + ttos(i) + "*" + B + ttos(j);
		}

		return buffer;
	}

	static PsimagLite::String createLocal(const PsimagLite::String& A,
	                                      SizeType n)
	{
		if (n < 1)
			err("createLocal needs at least one site\n");
		PsimagLite::String buffer;
		for (SizeType i = 0; i < n; ++i) {
			if (i > 0)
				buffer += "+";
			buffer += A + ttos(i);
		}

		return buffer;
	}

	static void addMatrixWithWeight(SparseMatrixType& m, RealType weight, const SparseMatrixType& a)
	{
		m += weight * a;
	}

	TypeEnum hamTipo;
	SizeType bits_;
	bool periodic_;
	IsingGraphType* isingGraph_;
	SparseMatrixType matrix_;
	VectorSizeType basis_;
	bool needsTransformAndTruncate_;
};
}
#endif // EVENDIM_HAMILTONIAN_CPU_H
