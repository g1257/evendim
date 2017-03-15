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
