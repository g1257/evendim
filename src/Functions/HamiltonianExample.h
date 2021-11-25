#ifndef HAMILTONIANEXAMPLE_H
#define HAMILTONIANEXAMPLE_H
#include "InputNg.h"
#include "InputCheck.h"

namespace Gep {

class HamiltonianExample {

public:

	typedef PsimagLite::InputNg<InputCheck> InputNgType;

	HamiltonianExample(typename InputNgType::Readable& io)
	{}

};
}
#endif // HAMILTONIANEXAMPLE_H
