#ifndef SCHWINGERMODEL_HH
#define SCHWINGERMODEL_HH
#include "CrsMatrix.h"

// 2308.04481 Eq. (2)
namespace Gep {

template <typename ComplexType>
class SchwingerModel {
public:

	using RealType = typename PsimagLite::Real<ComplexType>::Type;
	using SparseMatrixType = PsimagLite::CrsMatrix<ComplexType>;

	enum class Spin { UP,
		          DOWN };

	class State {
	public:

		State(SizeType ind)
		    : ind_(ind)
		{
		}

		Spin operator[](SizeType i) const
		{
			SizeType mask = (1 << i);
			return (mask & ind_) ? Spin::DOWN : Spin::UP;
		}

	private:

		SizeType ind_;
	};

	SchwingerModel(SizeType bits, RealType param_m, RealType param_g)
	    : bits_(bits)
	{
		SizeType hilbertSpace = (1 << bits);
		matrix_.resize(hilbertSpace, hilbertSpace);

		SizeType counter = 0;
		for (SizeType i = 0; i < hilbertSpace; ++i) {
			matrix_.setRow(i, counter);

			State state(i);
			ComplexType val = getMassTerm(state, param_m) + getGterm(state, param_g);
			if (std::abs(val) == 0.)
				continue;
			matrix_.pushCol(i);
			matrix_.pushValue(val);
			++counter;
		}

		matrix_.setRow(hilbertSpace, counter);
		matrix_.checkValidity();
	}

	const SparseMatrixType& matrix() const { return matrix_; }

private:

	RealType getMassTerm(const State& state, RealType param_m) const
	{
		SizeType twoL = bits_;
		double mass_term = twoL; // identity operator
		for (SizeType i = 0; i < twoL; ++i) {
			int sign = (i & 1) ? -1 : 1;
			int z = (state[i] == Spin::UP) ? 1 : -1;
			mass_term += sign * z;
		}

		mass_term *= param_m * 0.5;
		return mass_term;
	}

	RealType getGterm(const State& state, RealType param_g) const
	{
		SizeType twoL = bits_;
		double g_term = 0.;
		// FIXME: CHECK LIMIT OF THIS FOR LOOP
		for (SizeType i = 0; i < twoL - 1; ++i) {
			double qterm = sumOfQs(state, i);
			g_term += qterm * qterm;
		}

		g_term *= param_g * param_g;
		return g_term;
	}

	static RealType sumOfQs(const State& state, SizeType jnd)
	{
		RealType sum = 0.;
		for (SizeType i = 0; i <= jnd; ++i) { // note the <=
			int sign = (i & 1) ? -1 : 1;
			int z = (state[i] == Spin::UP) ? 1 : -1;
			sum += (sign + z);
		}

		return -sum * 0.5;
	}

	SizeType bits_;
	SparseMatrixType matrix_;
};
}

#endif // SCHWINGERMODEL_HH
