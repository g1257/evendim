#ifndef ISINGGRAPH_HH
#define ISINGGRAPH_HH
#include "Vector.h"
#include "Graph.hh"
#include "CrsMatrix.h"
#include "../Engine/ProgramGlobals.h"

namespace Gep {

template<typename ComplexType>
class IsingGraph {

public:

	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef Graph GraphType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;

	IsingGraph(SizeType bits, RealType coupling, bool periodic, PsimagLite::String graphFile)
	    : bits_(bits),
	      coupling_(coupling),
	      graph_(graphFile, bits, periodic)
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
			for (SizeType site = 0; site < bits_ - 1; ++site) {
				SizeType maskSite = (1 << site);
				SizeType j = i & maskSite;
				for (SizeType site2 = site + 1; site2 < bits_; ++site2) {
					if (!graph_.connected(site, site2)) continue;
					SizeType maskSite2 = (1 << site2);
					SizeType k = i & maskSite2;
					SizeType jj = (j > 0);
					SizeType kk = (k > 0);
					RealType tmp = PsimagLite::real(PsimagLite::conj(v[i])*v[i]);
					RealType value = (jj == kk) ? tmp : -tmp;
					e += value;
				}
			}
		}

		return e*coupling_;
	}

	PsimagLite::String info(const VectorType& v) const
	{
		const SizeType n = v.size();
		PsimagLite::String buffer;
		for (SizeType i = 0; i < n; ++i) {
			if (std::norm(v[i]) > 1e-4) buffer += ttos(i) + " ";
		}

		return buffer;
	}

private:

	SizeType bits_;
	RealType coupling_;
	GraphType graph_;
};
}
#endif // ISINGGRAPH_HH
