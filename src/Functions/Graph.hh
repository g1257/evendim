#ifndef GRAPH_HH
#define GRAPH_HH
#include "PsimagLite.h"

namespace Gep {

class Graph {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;

	Graph(PsimagLite::String graphFile, SizeType vertices = 0, bool periodic = false)
	    : graphFile_(graphFile), vertices_(vertices)
	{
		if (graphFile == "zz") {
			createChain(periodic);
			return;
		}

		err("Graph::ctor() unimplemented\n");
	}

	SizeType vertices() const
	{
		return vertices_;
	}

	const VectorSizeType& neighbors(SizeType site) const
	{
		assert(site < allNeighbors_.size());
		return allNeighbors_[site];
	}

private:

	void createChain(bool periodic)
	{
		for (SizeType vertex = 0; vertex < vertices_; ++vertex) {
			SizeType nextSite = vertex + 1;
			if (nextSite == vertices_) {
				if (!periodic) {
					allNeighbors_.push_back(VectorSizeType());
					break;
				}

				nextSite = 0;
			}

			VectorSizeType tmpVector(1, nextSite);
			allNeighbors_.push_back(tmpVector);
		}
	}

	PsimagLite::String graphFile_;
	SizeType vertices_;
	VectorVectorSizeType allNeighbors_;
};
}
#endif // GRAPH_HH
