#ifndef OPTIONS_H
#define OPTIONS_H
#include "InputNg.h"
#include "PsimagLite.h"
#include <algorithm>
#include <cctype>
#include <numeric>

namespace Gep {

class Options {

public:

	typedef typename PsimagLite::String::value_type CharType;
	typedef typename PsimagLite::String::const_iterator StringConstIterator;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	Options(PsimagLite::String str)
	{
		PsimagLite::split(vdata_, str, ",");
		lowerAll();
	}

	void operator+=(PsimagLite::String moreData)
	{
		VectorStringType vmore;
		PsimagLite::split(vmore, moreData, ",");
		vdata_.insert(vdata_.end(), vmore.begin(), vmore.end());
		lowerAll();
	}

	bool isSet(PsimagLite::String what) const
	{
		what = toLower(what);
		VectorStringType::const_iterator it = std::find(vdata_.begin(),
		                                                vdata_.end(),
		                                                what);
		return (it != vdata_.end());
	}

private:

	static PsimagLite::String toLower(PsimagLite::String data)
	{
		std::transform(data.begin(), data.end(), data.begin(), [](unsigned char c) { return std::tolower(c); });
		return data;
	}

	void lowerAll()
	{
		std::transform(vdata_.begin(),
		               vdata_.end(),
		               vdata_.begin(),
		               [](PsimagLite::String s) { return toLower(s); });
	}

	VectorStringType vdata_;
};
}
#endif // OPTIONS_H
