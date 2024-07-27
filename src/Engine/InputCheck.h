#ifndef INPUTCHECK_H
#define INPUTCHECK_H
#include "PsimagLite.h"

namespace Gep {

class InputCheck {

public:

	PsimagLite::String import() const
	{
		PsimagLite::String str;

		str += "integer HeadSize;\n";
		str += "integer Population;\n";
		str += "integer Generations;\n";
		str += "integer NumberOfBits;\n";
		str += "integer Samples;\n";
		str += "integer MinimizerVerbose;\n";
		str += "integer ProgressBar;\n";
		str += "integer HamiltonianIsPeriodic;\n";

		str += "string Primitives;\n";
		str += "string MinimizerAlgorithm;\n";
		str += "string RunType;\n";
		str += "string Hamiltonian;\n";
		str += "string InVectorFile;\n";

		str += "real MinimizerTolerance;\n";
		str += "real MinimizerDelta;\n";
		str += "real MinimizerDelta2;\n";
		str += "real HamiltonianCoupling;\n";
		str += "string GraphFile;\n";
		str += "integer Threads;\n";
		str += "integer Genes;\n";
		str += "integer Chead;\n";
		str += "string EngineOptions;\n";
		str += "vector Basis;\n";

		str += "integer Descendants;\n";
		str += "integer Mutations;\n";
		str += "integer Inversions;\n";
		str += "string UseXaccOptimizer;";
		str += "string XaccVerbosityLevel;";
		return str;
	}

	// POD inputs
	bool checkSimpleLabel(const PsimagLite::String& label,
	                      SizeType line) const
	{
		return true;
	}

	// SolverOptions= check
	void check(const PsimagLite::String& label,
	           const PsimagLite::String& val,
	           SizeType)
	{
	}

	bool check(const PsimagLite::String& label,
	           const PsimagLite::Vector<PsimagLite::String>::Type& vec,
	           SizeType line) const
	{
		return true;
	}
};
}
#endif // INPUTCHECK_H
