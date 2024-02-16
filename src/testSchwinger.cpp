#include "Fitness/SchwingerModel.hh"

template <typename ComplexType>
void printGs(PsimagLite::Matrix<ComplexType>& mat)
{
	assert(mat.rows() == mat.cols());
	std::vector<double> eigs(mat.rows());
	diag(mat, eigs, 'V');
	std::cout << "Ground State Energy=" << eigs[0] << "\n";
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

int main(int argc, char* argv[])
{
	using ComplexType = std::complex<double>;

	SizeType bits = 8;
	bool periodic = false;
	double m = 0.5;
	double g = 0.3;

	Gep::SchwingerModel<ComplexType> schwinger(bits, periodic, m, g);
	PsimagLite::Matrix<ComplexType> a = schwinger.matrix().toDense();
	printGs(a);
}
