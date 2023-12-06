#ifndef PROGRAMGLOBALS_H
#define PROGRAMGLOBALS_H
#include "PsimagLite.h"
#include "Vector.h"

namespace Gep {

namespace ProgramGlobals {
	template <typename SomeType, typename SomeRngType>
	void randomVector(std::vector<SomeType>& outVector,
	                  SomeRngType& rng,
	                  const SomeType& a,
	                  const SomeType& b)
	{
		typedef typename PsimagLite::Real<SomeType>::Type RealType;

		const SizeType n = outVector.size();
		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i) {
			SomeType val = a * rng() + b;
			outVector[i] = val;
			sum += PsimagLite::real(val * PsimagLite::conj(val));
		}

		assert(sum > 0);
		RealType factor = 1 / sqrt(sum);
		for (SizeType i = 0; i < n; ++i)
			outVector[i] *= factor;
	}

	template <typename SomeType>
	void writeVector(std::ostream& os, const std::vector<SomeType>& outVector)
	{
		const SizeType n = outVector.size();
		os << n << "\n";
		for (SizeType i = 0; i < n; ++i) {
			SomeType val = (PsimagLite::norm(outVector[i]) < 1e-6) ? 0 : outVector[i];
			os << val << " ";
		}

		os << "\n";
	}

	void pushVector(PsimagLite::Vector<PsimagLite::String>::Type& dest,
	                const PsimagLite::Vector<PsimagLite::String>::Type& src,
	                SizeType upTo = 0)
	{
		const SizeType total = src.size();
		if (upTo == 0)
			upTo = total;
		if (upTo > total)
			throw std::runtime_error("pushVector\n");

		for (SizeType j = 0; j < upTo; ++j)
			dest.push_back(src[j]);
	}

	PsimagLite::String vecStrToStr(const PsimagLite::Vector<PsimagLite::String>::Type& vecStr,
	                               PsimagLite::String sep)
	{
		PsimagLite::String ret;
		const SizeType n = vecStr.size();
		for (SizeType i = 0; i < n; ++i)
			ret += vecStr[i] + sep;
		return ret;
	}

	template <typename SomeType>
	static void readVector(std::vector<SomeType>& inVector, PsimagLite::String vectorFilename)
	{
		std::ifstream fin(vectorFilename);
		if (!fin || !fin.good() || fin.bad())
			err("Could not open file " + vectorFilename + "\n");
		int x = 0;
		fin >> x;
		if (x <= 0) {
			fin.close();
			err("First entry of file " + vectorFilename + " should be vector size\n");
		}

		inVector.resize(x);
		int i = 0;
		for (; i < x; ++i) {
			fin >> inVector[i];
			if (fin.eof())
				break;
		}

		if (i == x)
			return;

		fin.close();
		err("File " + vectorFilename + " should contain " + ttos(x) + " vector entries.\n");
	}
}
}
#endif // PROGRAMGLOBALS_H
