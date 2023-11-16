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
#include "Engine.h"
#include "Evolution.h"
#include "Fitness/Example1Fitness.h"
#include "Fitness/Example2Fitness.h"
#include "Fitness/Example3Fitness.h"
#include "Primitives/PlusMinusMultiplyDivide.h"
#include <unistd.h>

/* PSIDOC EngineOverviewFunction
 The main loop in gep2 is
\begin{lstlisting}
  // total = number of generations
  for (SizeType i = 0; i < total; ++i)
    engine.evolve(i);
\end{lstlisting}

 The \verb!Engine! class is templated on a \verb!Fitness! template that represents
 the training class, and determines how fit a GEP individual is.
 It is also templated on the Evolution type.
 The \verb!Engine! constructor takes an input parameters object, and an evolution object.
 */
template <template <typename> class FitnessTemplate,
          typename EvolutionType>
void main1(EvolutionType& evolution, const Gep::ParametersInput& gepOptions)
{
	typedef Gep::Engine<FitnessTemplate, EvolutionType> EngineType;

	typename EngineType::ParametersEngineType params(gepOptions);
	EngineType engine(params, evolution);

	for (SizeType i = 0; i < gepOptions.generations; ++i) {
		engine.evolve(i);
	}
}

/* PSIDOC Gep2main
        This driver program named gep2 runs different ``example'' cases consisting
        of using GEP to find a function knowning only some inputs and outputs.

        The primitives are under Primitives/PlusMinusMultiplyDivide.h, and
        the functions under Example1.h, Example2.h and Example3.h

        The following command line arguments to gep2 are mandatory.
        \begin{itemize}
        \item[-i] inputs. The number of inputs to the function.
        \item[-h] headSize. The maximum number of the head or effective gene size.
        \item[-p] population. The number of GEP individuals to consider in each generation.
        \item[-t] generations. The number of generations to run GEP.
        \end{itemize}

        The following command line arguments to gep2 are optional.
        \begin{itemize}
        \item[-e] example. The example number to run: 1, 2 or 3. Defaults to 1.
        \item[-g] genes. The number of genes to be used. Defaults to 1.
        \item[-s] seed. The seed for the random number generator.  Defaults to 1234.
        \item[-c] constants. The number of GEP constants to use. Default to 0.
        \item[-H] maximum head size for ADF. ADF stands for automatic defined funtions. Defaults to 0.
        \item[-a] adfs. The number of ADFs to use. Defaults to 0.
        \item[-n] samples. The number of training samples to cache. Defaults to 100.
        \item[-v] indicates that GEP should be verbose. Defaults to false.
        \item[-S] indicates that GEP should stop when a perfect individual is found. Defaults to false.
        \end{itemize}

        The options -v and -S take no arguments.
 */
int main(int argc, char* argv[])
{
	SizeType inputs = 0;
	SizeType seed = 1234;
	SizeType constants = 0;
	SizeType example = 0;
	bool verbose = false;
	Gep::ParametersInput gepOptions;

	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -i inputs -h head [-p population -t total -g genes -H chead]\n";
	while ((opt = getopt(argc, argv, "i:h:g:s:p:t:c:H:a:e:n:Sv")) != -1) {
		switch (opt) {
		case 'i':
			inputs = atoi(optarg);
			break;
		case 'h':
			gepOptions.head = atoi(optarg);
			break;
		case 'g':
			gepOptions.genes = atoi(optarg);
			break;
		case 's':
			seed = atoi(optarg);
			break;
		case 'p':
			gepOptions.population = atoi(optarg);
			break;
		case 't':
			gepOptions.generations = atoi(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		case 'c':
			constants = atoi(optarg);
			break;
		case 'H':
			gepOptions.chead = atoi(optarg);
			break;
		case 'a':
			gepOptions.adfs = atoi(optarg);
			break;
		case 'e':
			example = atoi(optarg);
			break;
		case 'n':
			gepOptions.samples = atoi(optarg);
			break;
		case 'S':
			*gepOptions.options += "stopEarly";
			break;
		default:
			throw PsimagLite::RuntimeError(strUsage);
			return 1;
		}
	}

	// sanity checks here
	if (inputs == 0
	    || gepOptions.head == 0
	    || gepOptions.population == 0
	    || gepOptions.generations == 0) {
		throw PsimagLite::RuntimeError(strUsage);
		return 1;
	}

	if (gepOptions.chead > 0 && gepOptions.adfs == 0)
		throw PsimagLite::RuntimeError("FATAL: You selected ADF head size H > 0 but ADF number a == 0\n");

	if (gepOptions.chead > 0 && gepOptions.adfs == 0)
		throw PsimagLite::RuntimeError("FATAL: You selected ADF number a > 0 but ADF head size H == 0\n");

	bool hasAdfs = (gepOptions.chead > 0 && gepOptions.adfs > 0);
	if (gepOptions.genes > 1 && !hasAdfs)
		throw PsimagLite::RuntimeError("FATAL: genes > 1 but no ADF\n");

	/* PSIDOC EvolutionInFunction
 Evolution is templated on Primitives, which represents the GEP primitives or ``operators''
 to be considered. Evolution's constructor takes a primitives object, a seed, and an verbose
 boolean. In this file, gep2.cpp, Primitives is set to the class PlusMinusMultiplyDivide so that
 the primitives are plus, minus, multiply and divide.
	 */
	typedef Gep::PlusMinusMultiplyDivide<double> PrimitivesType;
	typedef Gep::Evolution<PrimitivesType> EvolutionType;

	PrimitivesType primitives(inputs, gepOptions.genes, constants);
	EvolutionType evolution(primitives, seed, verbose);

	if (example < 2) {
		main1<Gep::Example1Fitness, EvolutionType>(evolution, gepOptions);
		return 0;
	}
	else if (example == 2) {
		main1<Gep::Example2Fitness, EvolutionType>(evolution, gepOptions);
		return 0;
	}

	main1<Gep::Example3Fitness, EvolutionType>(evolution, gepOptions);
}
