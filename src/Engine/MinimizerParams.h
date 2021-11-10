#ifndef MINIMIZERPARAMS_H
#define MINIMIZERPARAMS_H
#include "Vector.h"

namespace Gep {

template<typename RealType_>
struct MinimizerParams {

	typedef RealType_ RealType;

	enum EnumAlgo {SIMPLEX, CONJUGATE_GRADIENT};

	MinimizerParams(EnumAlgo algo_,
	                SizeType maxIter_,
	                RealType delta_,
	                RealType delta2_,
	                RealType tol_,
	                SizeType saveEvery_,
	                bool verbose_)
	    : verbose(verbose_),
	      algo(algo_),
	      maxIter(maxIter_),
	      saveEvery(saveEvery_),
	      delta(delta_),
	      delta2(delta2_),
	      tol(tol_)
	{}

	bool verbose;
	EnumAlgo algo;
	SizeType maxIter;
	SizeType saveEvery;
	RealType delta;
	RealType delta2;
	RealType tol;
};


template<typename RealType>
std::ostream& operator<<(std::ostream& os, const MinimizerParams<RealType>& m)
{
	os<<"algo= "<<m.algo<<"\n";
	os<<"maxIter= "<<m.maxIter<<"\n";
	os<<"delta= "<<m.delta<<"\n";
	os<<"delta2= "<<m.delta2<<"\n";
	os<<"tolerance= "<<m.tol<<"\n";
	os<<"verbose= "<<m.verbose<<"\n";
	return os;
}
} // namespace Gep

#endif // MINIMIZERPARAMS_H

