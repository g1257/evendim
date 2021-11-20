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
		str += "real MinimizerTolerance;\n";
		str += "string Primitives;\n";
		str += "string MinimizerAlgorithm;\n";
		str += "real MinimizerDelta;\n";
		str += "real MinimizerDelta2;\n";
		str += "integer Samples;\n";
		str += "integer MinimizerVerbose;\n";
		str += "integer ProgressBar;\n";
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
	{}

	bool check(const PsimagLite::String& label,
	           const PsimagLite::Vector<PsimagLite::String>::Type& vec,
	           SizeType line) const
	{
		return true;
	}


};
}
#endif // INPUTCHECK_H
