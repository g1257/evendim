#include "../../PsimagLite/src/LanczosSolver.h"
#include "Fitness/SchwingerModel.hh"

using ComplexType = std::complex<double>;
using VectorRealType = std::vector<double>;
using VectorType = std::vector<ComplexType>;

double observable(const VectorType& v, const VectorRealType& obs)
{
	double sum = 0.;
	SizeType n = v.size();
	assert(n == obs.size());
	for (SizeType i = 0; i < n; ++i) {
		sum += std::real(std::conj(v[i]) * obs[i] * v[i]);
	}

	return sum;
}

double observable(const PsimagLite::Matrix<ComplexType>& a, const VectorRealType& obs)
{
	double sum = 0.;
	SizeType n = a.rows();
	assert(n == obs.size());
	for (SizeType i = 0; i < n; ++i) {
		sum += std::real(std::conj(a(i, 0)) * obs[i] * a(i, 0));
	}

	return sum;
}

void solveSparse(const PsimagLite::CrsMatrix<ComplexType>& matrix,
                 const VectorRealType& chi,
                 double factor_for_density)
{
	using SparseMatrixType = PsimagLite::CrsMatrix<ComplexType>;
	using VectorType = std::vector<ComplexType>;
	using SolverParametersType = PsimagLite::ParametersForSolver<double>;
	using LanczosSolverType = PsimagLite::LanczosSolver<SolverParametersType, SparseMatrixType, VectorType>;

	int n = matrix.rows();

	SolverParametersType params;
	params.lotaMemory = true;
	params.minSteps = std::min(n, 100);

	LanczosSolverType lanczos_solver(matrix, params);

	double energy = 0.;
	VectorType z(n);
	VectorType init_vector(n);
	PsimagLite::fillRandom(init_vector);
	lanczos_solver.computeOneState(energy, z, init_vector, 0);
	std::cout << "Energy= " << energy << " density= " << energy * factor_for_density << "\n";
	std::cout << "Chi = " << observable(z, chi) << "\n";
}

void solveDense(const PsimagLite::CrsMatrix<ComplexType>& matrix,
                const VectorRealType& chi,
                double factor_for_density)
{
	PsimagLite::Matrix<ComplexType> a = matrix.toDense();
	if (!isHermitian(a, true)) {
		throw std::runtime_error("Not Hermitian\n");
	}

	using VectorRealType = std::vector<double>;
	VectorRealType eigs(a.rows());
	diag(a, eigs, 'V');
	double energy = eigs[0];
	std::cout << "Energy= " << energy << " density= " << energy * factor_for_density << "\n";
	std::cout << "Chi = " << observable(a, chi) << "\n";
}

int main(int argc, char* argv[])
{

	if (argc != 2) {
		std::cerr << "USAGE: " << argv[0] << " bits\n";
		return 1;
	}

	SizeType bits = std::stoi(argv[1]);
	bool periodic = false;
	double m = 0.5;
	double g = 0.3;

	Gep::SchwingerModel<ComplexType> schwinger(bits, periodic, m, g);
	double factor_for_density = 2. / bits;
	if (bits > 6) {
		solveSparse(schwinger.matrix(), schwinger.chi(), factor_for_density);
	}

	if (bits < 12) {
		solveDense(schwinger.matrix(), schwinger.chi(), factor_for_density);
	}
}
