#ifndef HAMILTONIANEXAMPLE_H
#define HAMILTONIANEXAMPLE_H
#include "InputNg.h"
#include "InputCheck.h"
#include "PsimagLite.h"
#include "CrsMatrix.h"

namespace Gep {

// H = coupling*\sum_{i} sigma^x_i sigma^x_{i + 1}
template<typename ComplexType>
class HamiltonianExample {

public:

	typedef PsimagLite::InputNg<InputCheck> InputNgType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;

	HamiltonianExample(typename InputNgType::Readable& io)
	{
		SizeType bits = 0;
		io.readline(bits, "NumberOfBits="); // == number of "sites"

		bool periodic = false;
		try {
			int tmp = 0;
			io.readline(tmp, "HamiltonianIsPeriodic=");
			periodic = (tmp > 0);
		} catch (std::exception&) {}

		RealType coupling = 0;
		io.readline(coupling, "HamiltonianCoupling=");

		SizeType hilbertSpace = (1 << bits);
		matrix_.resize(hilbertSpace, hilbertSpace);
		cacheVector_.resize(hilbertSpace);

		VectorRealType v(hilbertSpace);
		VectorBoolType bcol(hilbertSpace);

		SizeType counter = 0;
		for (SizeType i = 0; i < hilbertSpace; ++i) {
			matrix_.setRow(i, counter);
			std::fill(v.begin(), v.end(), 0);
			std::fill(bcol.begin(), bcol.end(), false);
			for (SizeType site = 0; site < bits; ++site) {
				// Flip bit at site
				SizeType maskSite = (1 << site);
				SizeType j = i ^ maskSite;
				SizeType site2 = site + 1;
				if (site2 >= bits && !periodic) continue;
				assert(site2 <= bits);
				if (site2 == bits) site2 = 0;
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
		PsimagLite::Matrix<RealType> a = matrix_.toDense();
		diag(a, eigs, 'V');
		std::cout<<"Ground State Energy="<<eigs[0]<<"\n";

	}

	RealType energy(const VectorType& y) const
	{
		std::fill(cacheVector_.begin(), cacheVector_.end(), 0);
		matrix_.matrixVectorProduct(cacheVector_, y);
		return PsimagLite::real(y*cacheVector_); // does conjugation of first vector
	}

private:

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
	SparseMatrixType matrix_;
	mutable VectorType cacheVector_;
};
}
#endif // HAMILTONIANEXAMPLE_H
