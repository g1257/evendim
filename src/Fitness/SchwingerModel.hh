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
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;

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

	SchwingerModel(SizeType bits, bool periodic, RealType param_m, RealType param_g)
	    : bits_(bits)
	    , periodic_(periodic)
	{
		constexpr double coupling = 1. / 2.; // s+ s- coupling constant

		SizeType hilbertSpace = (1 << bits);
		matrix_.resize(hilbertSpace, hilbertSpace);

		VectorRealType v(hilbertSpace);
		VectorBoolType bcol(hilbertSpace);

		SizeType counter = 0;

		for (SizeType i = 0; i < hilbertSpace; ++i) {
			matrix_.setRow(i, counter);

			State state(i);

			// diagonal terms
			ComplexType val = getMassTerm(state, param_m) + getGterm(state, param_g);
			if (std::abs(val) != 0.) {

				matrix_.pushCol(i);
				matrix_.pushValue(val);
				++counter;
			}

			// off-diagonal term
			SizeType total = bits_;
			for (SizeType site = 0; site < total; ++site) {
				SizeType site2 = site + 1;
				if (site2 >= total && !periodic_)
					continue;
				assert(site2 <= total);
				if (site2 == total)
					site2 = 0;

				// up up and down down states do not contribute
				if (state[site2] == state[site])
					continue;

				// Flip bit at site
				SizeType maskSite = (1 << site);
				SizeType j = i ^ maskSite;

				// Flip bit at site2
				SizeType maskSite2 = (1 << site2);
				j ^= maskSite2;
				v[j] += coupling;
				bcol[j] = true;
			}

			counter += fillThisRow(v, bcol);
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
		for (SizeType i = 0; i < twoL; ++i) {
			double qterm = sumOfQs(state, i);
			g_term += qterm * qterm;
		}

		g_term *= param_g * param_g * 0.5;
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

	SizeType bits_;
	bool periodic_;
	SparseMatrixType matrix_;
};
}

#endif // SCHWINGERMODEL_HH
