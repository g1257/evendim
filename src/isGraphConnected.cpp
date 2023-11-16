#include "Fitness/Graph.hh"

int main(int argc, char* argv[])
{
	if (argc != 2)
		err("USAGE: " + PsimagLite::String(argv[0]) + " filename\n");

	Gep::Graph graph("file:" + PsimagLite::String(argv[1]));
	std::cout << graph.isConnected() << "\n";
}
