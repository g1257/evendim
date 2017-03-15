#ifndef PARAMETERSENGINE_H
#define PARAMETERSENGINE_H

namespace Gep {

struct Options
{
	Options(SizeType p = 0,
	        SizeType h = 0,
	        SizeType g = 1,
	        SizeType ch = 0,
	        SizeType adfs1 = 0,
	        bool se = false)
	    : population(p),head(h),genes(g),chead(ch),adfs(adfs1),stopEarly(se)
	{}

	SizeType population;
	SizeType head;
	SizeType genes;
	SizeType chead;
	SizeType adfs;
	bool stopEarly;
};

template<typename RealType>
class ParametersEngine {

public:

	ParametersEngine(const Options& op,
	                 RealType d = 2.0,
	                 RealType m = 0.5,
	                 RealType i = 0.5,
	                 SizeType s = 50)
	    : population(op.population),
	      head(op.head),
	      genes(op.genes),
	      chead(op.chead),
	      adfs(op.adfs),
	      descendants(static_cast<SizeType>(op.population*d)),
	      mutation(static_cast<SizeType>(op.population*m)),
	      inversion(static_cast<SizeType>(op.population*i)),
	      samples(s),
	      stopEarly(op.stopEarly)
	{}

	SizeType population;
	SizeType head;
	SizeType genes;
	SizeType chead;
	SizeType adfs;
	SizeType descendants;
	SizeType mutation;
	SizeType inversion;
	SizeType samples;
	bool stopEarly;

}; // class ParametersEngine

} // namespace Gep

#endif // PARAMETERSENGINE_H
