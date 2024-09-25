#ifndef ISINGGRAPH_HH
#define ISINGGRAPH_HH
#include "../Engine/ProgramGlobals.h"
#include "CrsMatrix.h"
#include "Graph.hh"
#include "Vector.h"

namespace Gep {

template <typename ComplexType>
class IsingGraph {

public:

	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef Graph GraphType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef
	    typename PsimagLite::Vector<ComplexType>::Type VectorType;

	IsingGraph(SizeType bits, RealType coupling, bool periodic, PsimagLite::String graphFile)
	    : bits_(bits)
	    , coupling_(coupling)
	    , graph_(graphFile, bits, periodic)
	{
		if (graph_.vertices() != bits)
			err("Graph vertices != bits\n");
	}

	RealType energyZZ(const VectorType& v) const
	{
		const SizeType hilbertSpace = v.size();
		RealType e = 0;
		assert(bits_ > 1);
		for (SizeType i = 0; i < hilbertSpace; ++i) {
			for (SizeType site = 0; site < bits_ - 1;
			     ++site) {
				SizeType maskSite = (1 << site);
				SizeType j = i & maskSite;
				for (SizeType site2 = site + 1;
				     site2 < bits_;
				     ++site2) {
					if (!graph_.connected(site,
					                      site2))
						continue;
					SizeType maskSite2 = (1 << site2);
					SizeType k = i & maskSite2;
					SizeType jj = (j > 0);
					SizeType kk = (k > 0);
					RealType tmp = PsimagLite::real(
					    PsimagLite::conj(v[i]) * v[i]);
					RealType value = (jj == kk) ? tmp : -tmp;
					e += value;
				}
			}
		}

		return e * coupling_;
	}

	// Needed by XACC
	std::string buildExpression() const
	{
		const std::string coupling_str = (coupling_ == 1) ? "" : ttos(coupling_) + "*";
		std::string str;
		bool firstCall = true;
		assert(bits_ > 1);
		for (SizeType site = 0; site < bits_ - 1; ++site) {
			for (SizeType site2 = site + 1; site2 < bits_; ++site2) {
				if (!graph_.connected(site, site2)) {
					continue;
				}

				if (!firstCall) {
					str += " + ";
				}
				else {
					firstCall = false;
				}

				str += coupling_str + "Sz" + ttos(site) + "*Sz" + ttos(site2);
			}
		}

		return str;
	}

	// Use only to obtain the exact solution
	void solve()
	{
		SizeType total = (1 << bits_);
		VectorType v(total);
		RealType emin = 0;
		std::vector<SizeType> ind;
		for (SizeType i = 0; i < total; ++i) {
			v[i] = 1;
			RealType e = this->energyZZ(v);
			if (e < emin || i == 0) {
				ind.resize(1, i);
				emin = e;
			}
			else if (e == emin) {
				ind.push_back(i);
			}
			v[i] = 0;
		}

		std::cout << "IsingGraph::emin=" << emin << "\n";
		for (SizeType i = 0; i < ind.size(); ++i) {
			std::cout << ind[i] << " ";
		}

		std::cout << "\n";
	}

private:

	SizeType bits_;
	RealType coupling_;
	GraphType graph_;
};
} // namespace Gep
#endif // ISINGGRAPH_HH
