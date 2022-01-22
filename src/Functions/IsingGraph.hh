#ifndef ISINGGRAPH_HH
#define ISINGGRAPH_HH
#include "Vector.h"
#include "Graph.hh"
#include "CrsMatrix.h"

namespace Gep {

template<typename ComplexType>
class IsingGraph {

public:

	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef Graph GraphType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;

	IsingGraph(SizeType bits, RealType coupling, PsimagLite::String graphFile)
	    : FILE_PREFIX("file:"),
	      bits_(bits),
	      coupling_(coupling),
	      graph_(graphFile)
	{
		if (graph_.vertices() != bits)
			err("Graph vertices != bits\n");

		const bool fromFile = (graphFile.substr(0, 5) == FILE_PREFIX);
		if (fromFile)
			err("Unimplemented\n");
	}

	RealType energyZZ(const VectorType& v, bool periodic) const
	{
		const SizeType hilbertSpace = v.size();
		RealType e = 0;
		for (SizeType i = 0; i < hilbertSpace; ++i) {
			for (SizeType site = 0; site < bits_; ++site) {
				SizeType maskSite = (1 << site);
				SizeType j = i & maskSite;
				SizeType site2 = site + 1;
				if (site2 >= bits_ && !periodic) continue;
				assert(site2 <= bits_);
				if (site2 == bits_) site2 = 0;
				SizeType maskSite2 = (1 << site2);
				SizeType k = i & maskSite2;
				SizeType jj = (j > 0);
				SizeType kk = (k > 0);
				RealType tmp = PsimagLite::real(PsimagLite::conj(v[i])*v[i]);
				RealType value = (jj == kk) ? tmp : -tmp;
				e += value;
			}
		}

		return e*coupling_;
	}

private:

	const PsimagLite::String FILE_PREFIX;
	SizeType bits_;
	RealType coupling_;
	GraphType graph_;
};
}
#endif // ISINGGRAPH_HH
