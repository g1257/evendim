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
#include "Options.h"

namespace Gep {

struct ParametersInput {

	typedef PsimagLite::InputNg<InputCheck> InputNgType;

	ParametersInput(SizeType gen = 0,
                    SizeType p = 0,
	                SizeType h = 0,
	                SizeType g = 1,
	                SizeType ch = 0,
	                SizeType adfs1 = 0,
	                SizeType samples1 = 50,
	                SizeType threads1 = 1,
	                PsimagLite::String options1 = "",
	                PsimagLite::String prim = "")
	    : generations(gen),
	      population(p),
	      head(h),
	      genes(g),
	      chead(ch),
	      adfs(adfs1),
	      samples(samples1),
	      threads(threads1),
	      options(new Options(options1)),
	      primitives(prim)
	{}

	ParametersInput(InputNgType::Readable& io)
	    : generations(0),
	      population(0),
	      head(0),
	      genes(1),
	      chead(0),
	      adfs(0),
	      samples(50),
	      threads(1),
	      options(nullptr),
	      primitives("")
	{
/* PSIDOC ParamtersEngineInFunction
The engine parameters can be specified with
\verb!Generations=100;! in the input file,
and similarly for the others, which are as follows.
\begin{itemize}
\item[Generations] The number of GEP generations. Integer. Mandatory.
\item[Population] The number of GEP individuals. Integer. Mandatory.
\item[HeadSize] The size of the head (that is, the maximum effective gene size). Integer. Mandatory.
\item[Samples] The samples to be cached. Optional. Defaults to 50 and is unused in quantumGEP.
\item[Threads] The number of shared memory threads to use. Optional. Defaults to 1.
Not all fitness classes support paralellization, that is, a number greater than one here.
\item[Primitives] A comma-separated list of quantum gates to consider by GEP. String. Optional.
Defaults to "C,H,P".
\item[EngineOptions] A comma-separated list of options. String. Optional. Default to the empty string.
\end{itemize}

The EngineOptions are case-insensitive and can be none or more of the following.
\begin{itemize}
\item[stopEarly] Stops quantumGEP as soon as a perfect individual (that is, circuit) is found.
\item[noncanonical] Disables the canonicalization step.
\item[progressBar] Prints a progress bar for each generation.
\item[printCompact] Prints individuals in compact form.
\end{itemize}
*/
		io.readline(generations, "Generations=");

		io.readline(population, "Population=");

		io.readline(head, "HeadSize=");

		try {
			io.readline(samples, "Samples=");
		} catch (std::exception&) {}

		try {
			io.readline(threads, "Threads=");
		} catch (std::exception&) {}

		try {
			io.readline(primitives, "Primitives=");
		} catch (std::exception&) {}

		PsimagLite::String str;
		try {
			io.readline(str, "EngineOptions=");
		} catch (std::exception&) {}

		options = new Options(str);
	}

	~ParametersInput()
	{
		delete options;
		options = nullptr;
	}

	ParametersInput(const ParametersInput&) = delete;

	ParametersInput& operator=(const ParametersInput&) = delete;

	SizeType generations;
	SizeType population;
	SizeType head;
	SizeType genes;
	SizeType chead;
	SizeType adfs;
	SizeType samples;
	SizeType threads;
	Options* options;
	PsimagLite::String primitives; // comma-separated list of primitives
};

template<typename RealType>
class ParametersEngine {

public:

	ParametersEngine(const ParametersInput& op,
	                 RealType d = 2.0,
	                 RealType m = 0.5,
	                 RealType i = 0.5)
	    : generations(op.generations),
	      population(op.population),
	      head(op.head),
	      genes(op.genes),
	      chead(op.chead),
	      adfs(op.adfs),
	      descendants(static_cast<SizeType>(op.population*d)),
	      mutation(static_cast<SizeType>(op.population*m)),
	      inversion(static_cast<SizeType>(op.population*i)),
	      samples(op.samples),
	      threads(op.threads),
	      options(*op.options)
	{}

	SizeType generations;
	SizeType population;
	SizeType head;
	SizeType genes;
	SizeType chead;
	SizeType adfs;
	SizeType descendants;
	SizeType mutation;
	SizeType inversion;
	SizeType samples;
	SizeType threads;
	const Options& options;
}; // class ParametersEngine

} // namespace Gep

#endif // PARAMETERSENGINE_H
