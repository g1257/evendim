#include "../../PsimagLite/src/LanczosSolver.h"
#include "Fitness/SchwingerModel.hh"

using ComplexType = std::complex<double>;

double solveSparse(const PsimagLite::CrsMatrix<ComplexType>& matrix)
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
	return energy;
}

double solveDense(const PsimagLite::CrsMatrix<ComplexType>& matrix)
{
	PsimagLite::Matrix<ComplexType> a = matrix.toDense();
	if (!isHermitian(a, true)) {
		throw std::runtime_error("Not Hermitian\n");
	}

	using VectorRealType = std::vector<double>;
	VectorRealType eigs(a.rows());
	diag(a, eigs, 'N');
	return eigs[0];
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
		double energy = solveSparse(schwinger.matrix());
		std::cout << "Energy= " << energy << " density= " << energy * factor_for_density << "\n";
	}

	if (bits < 12) {
		double energy = solveDense(schwinger.matrix());
		std::cout << "Energy= " << energy << " density= " << energy * factor_for_density << "\n";
	}
}
