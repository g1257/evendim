#ifndef BASEFITNESS_H
#define BASEFITNESS_H

namespace Gep {

class NullClass {};

template<typename EvolutionType>
class BaseFitness {
public:

	typedef NullClass FitnessParamsType;

	const SizeType status() const { return 0; }

};
}
#endif // BASEFITNESS_H
