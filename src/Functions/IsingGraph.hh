#ifndef ISINGGRAPH_HH
#define ISINGGRAPH_HH
#include "Vector.h"
#include "Graph.hh"
#include "CrsMatrix.h"

namespace Gep {

template<typename ComplexType>
class IsingGraph {

public:

	typedef Graph GraphType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;

	IsingGraph(SizeType bits, PsimagLite::String graphFile)
	    : graph_(graphFile)
	{
		if (graph_.vertices() != bits)
			err("Graph vertices != bits\n");
	}

	void fillMatrix(SparseMatrixType& matrix)
	{
		err("IsingGraph::fillMatrix(): unimplemented\n");
	}

private:

	GraphType graph_;
};
}
#endif // ISINGGRAPH_HH
