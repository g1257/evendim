#ifndef PROGRAMGLOBALS_H
#define PROGRAMGLOBALS_H
#include "Vector.h"

namespace Gep {

static void pushVector(PsimagLite::Vector<PsimagLite::String>::Type& dest,
                       const PsimagLite::Vector<PsimagLite::String>::Type& src,
                       SizeType upTo = 0)
{
	const SizeType total = src.size();
	if (upTo == 0) upTo = total;
	if (upTo > total)
		err("pushVector\n");

	for (SizeType j = 0; j < upTo; ++j)
		dest.push_back(src[j]);
}

static PsimagLite::String vecStrToStr(const PsimagLite::Vector<PsimagLite::String>::Type& vecStr,
                                      PsimagLite::String sep)
{
	PsimagLite::String ret;
	const SizeType n = vecStr.size();
	for (SizeType i = 0; i < n; ++i)
		ret += vecStr[i] + sep;
	return ret;
}

}
#endif // PROGRAMGLOBALS_H
