/*
Copyright (c) 2017, UT-Battelle, LLC

evendim, Version 0.

This file is part of evendim.
evendim is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
evendim is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with evendim. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef PARAMETERSENGINE_H
#define PARAMETERSENGINE_H
#include "PsimagLite.h"
#include "InputNg.h"
#include "InputCheck.h"

namespace Gep {

struct Options {

	typedef PsimagLite::InputNg<InputCheck> InputNgType;

	Options(SizeType p = 0,
	        SizeType h = 0,
	        SizeType g = 1,
	        SizeType ch = 0,
	        SizeType adfs1 = 0,
	        SizeType samples1 = 50,
	        bool se = false,
	        bool pb = false,
	        PsimagLite::String prim = "")
	    : population(p),
	      head(h),
	      genes(g),
	      chead(ch),
	      adfs(adfs1),
	      samples(samples1),
	      stopEarly(se),
	      progressBar(pb),
	      primitives(prim)
	{}

	Options(InputNgType::Readable& io)
	    : population(0),
	      head(0),
	      genes(1),
	      chead(0),
	      adfs(0),
	      samples(50),
	      stopEarly(false),
	      progressBar(false),
	      primitives("")
	{
		io.readline(population, "Population=");

		io.readline(head, "HeadSize=");

		try {
			io.readline(samples, "Samples=");
		} catch (std::exception&) {}

		try {
			int tmp = 0;
			io.readline(tmp, "StopEarly=");
			stopEarly = (tmp > 0);
		} catch (std::exception&) {}

		try {
			int tmp = 0;
			io.readline(tmp, "ProgressBar=");
			progressBar = (tmp > 0);
		} catch (std::exception&) {}

		try {
			io.readline(primitives, "Primitives=");
		} catch (std::exception&) {}
	}

	SizeType population;
	SizeType head;
	SizeType genes;
	SizeType chead;
	SizeType adfs;
	SizeType samples;
	bool stopEarly;
	bool progressBar;
	PsimagLite::String primitives; // comma-separated list of primitives
};

template<typename RealType>
class ParametersEngine {

public:

	ParametersEngine(const Options& op,
	                 RealType d = 2.0,
	                 RealType m = 0.5,
	                 RealType i = 0.5)
	    : population(op.population),
	      head(op.head),
	      genes(op.genes),
	      chead(op.chead),
	      adfs(op.adfs),
	      descendants(static_cast<SizeType>(op.population*d)),
	      mutation(static_cast<SizeType>(op.population*m)),
	      inversion(static_cast<SizeType>(op.population*i)),
	      samples(op.samples),
	      stopEarly(op.stopEarly),
	      progressBar(op.progressBar)
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
	bool progressBar;
}; // class ParametersEngine

} // namespace Gep

#endif // PARAMETERSENGINE_H
