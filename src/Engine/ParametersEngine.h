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
#include "InputCheck.h"
#include "InputNg.h"
#include "Options.h"
#include "PsimagLite.h"

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
	    : generations(gen)
	    , population(p)
	    , head(h)
	    , genes(g)
	    , chead(ch)
	    , adfs(adfs1)
	    , descendants(2 * p)
	    , mutation(0.5 * p)
	    , inversion(0.5 * p)
	    , samples(samples1)
	    , threads(threads1)
	    , options(new Options(options1))
	    , primitives(prim)
	{
	}

	ParametersInput(InputNgType::Readable& io)
	    : generations(0)
	    , population(0)
	    , head(0)
	    , genes(1)
	    , chead(0)
	    , adfs(0)
	    , samples(50)
	    , threads(1)
	    , options(nullptr)
	    , primitives("")
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

		descendants = 2 * population;
		inversion = 0.5 * population;
		mutation = 0.5 * population;

		io.readline(head, "HeadSize=");

		tryToRead(genes, "Genes=", io);
		if (genes > 1) {
			PsimagLite::String str("Automatically Setting ADFS to 1\n");
			std::cout << str;
			std::cerr << str;
			adfs = 1;
		}

		tryToRead(descendants, "Descendants=", io);
		tryToRead(inversion, "Inversions=", io);
		tryToRead(mutation, "Mutations=", io);

		tryToRead(chead, "Chead=", io);
		if (genes > 1 && chead == 0) {
			throw PsimagLite::RuntimeError("genes > 1 but chead == 0\n");
		}

		tryToRead(samples, "Samples=", io);
		tryToRead(threads, "Threads=", io);

		try {
			io.readline(primitives, "Primitives=");
		}
		catch (std::exception&) {
		}

		PsimagLite::String str;
		try {
			io.readline(str, "EngineOptions=");
		}
		catch (std::exception&) {
		}

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
	SizeType descendants;
	SizeType mutation;
	SizeType inversion;
	SizeType samples;
	SizeType threads;
	Options* options;
	PsimagLite::String primitives; // comma-separated list of primitives

private:

	static void tryToRead(SizeType& what, const std::string& label, InputNgType::Readable& io)
	{
		try {
			io.readline(what, label);
		}
		catch (std::exception&) {
		}
	}
};

template <typename RealType>
class ParametersEngine {

public:

	ParametersEngine(const ParametersInput& op)
	    : generations(op.generations)
	    , population(op.population)
	    , head(op.head)
	    , genes(op.genes)
	    , chead(op.chead)
	    , adfs(op.adfs)
	    , descendants(op.descendants)
	    , mutation(op.mutation)
	    , inversion(op.inversion)
	    , samples(op.samples)
	    , threads(op.threads)
	    , options(*op.options)
	{
	}

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
