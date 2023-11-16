#ifndef AUXFORHAMSPEC_H
#define AUXFORHAMSPEC_H
#include "PsimagLite.h"

namespace Gep {

class AuxForHamSpec {

public:

	AuxForHamSpec(SizeType numberOfBits)
	    : bits_(numberOfBits)
	{
	}

	SizeType numberOfBits() const { return bits_; }

private:

	SizeType bits_;
};
}
#endif // AUXFORHAMSPEC_H
