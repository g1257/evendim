#ifndef GRAPH_HH
#define GRAPH_HH
#include "PsimagLite.h"

namespace Gep {

class Graph {

public:

	Graph(PsimagLite::String graphFile)
	    : graphFile_(graphFile)
	{
		err("Graph::ctor() unimplemented\n");
	}

	SizeType vertices() const
	{
		throw PsimagLite::RuntimeError("Graph::vertices(): unimplemented\n");
	}

private:

	PsimagLite::String graphFile_;
};
}
#endif // GRAPH_HH
