#include "Fitness/Graph.hh"
#include "PsimagLite.h"

int main(int argc, char* argv[])
{
	if (argc != 2)
		err("USAGE: " + PsimagLite::String(argv[0]) + " order\n");

	const SizeType vertices = PsimagLite::atoi(argv[1]);

	if (vertices < 2)
		err(PsimagLite::String(argv[0]) + ": Expected at least two vertices\n");

	SizeType n = vertices * (vertices - 1) / 2;
	Gep::Graph::LongUintType nstates = (1 << n);
	std::cout << "Graph, order " << n << "\n";
	for (Gep::Graph::LongUintType state = 0; state < nstates; ++state) {
		Gep::Graph graph(state, vertices);
		if (!graph.isConnected())
			continue;
		std::cout << graph;
		std::cout << "----------------------------------------\n";
	}
}
