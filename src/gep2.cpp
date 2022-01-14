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
#include "Evolution.h"
#include "Primitives/PlusMinusMultiplyDivide.h"
#include "Engine.h"
#include <unistd.h>
#include "Functions/Example1.h"
#include "Functions/Example2.h"
#include "Functions/Example3.h"

template<template<typename> class FitnessTemplate,
         typename EvolutionType>
void main1(EvolutionType& evolution,
           const Gep::ParametersInput& gepOptions,
           SizeType total)
{
	typedef FitnessTemplate<EvolutionType> FitnessType;
	typedef Gep::Engine<FitnessType> EngineType;

	typename EngineType::ParametersEngineType params(gepOptions);
	EngineType engine(params, evolution);

	for (SizeType i = 0; i < total; i++)
		if (engine.evolve(i) && params.options.isSet("stopEarly")) break;
}

int main(int argc, char* argv[])
{
	SizeType inputs = 0;
	SizeType total = 0;
	SizeType seed = 1000;
	SizeType constants = 0;
	SizeType example = 0;
	bool verbose = false;
	Gep::ParametersInput gepOptions;

	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -i inputs -h head [-p population -t total -g genes -H chead]\n";
	while ((opt = getopt(argc, argv,"i:h:g:s:p:t:c:H:a:e:n:Sv")) != -1) {
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
			total = atoi(optarg);
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
	if (inputs == 0 || gepOptions.head == 0 || gepOptions.population == 0 || total == 0) {
		throw PsimagLite::RuntimeError(strUsage);
		return 1;
	}

	if (gepOptions.genes > 1 && (gepOptions.chead == 0 || gepOptions.adfs == 0))
		throw PsimagLite::RuntimeError(strUsage);

	typedef Gep::PlusMinusMultiplyDivide<double> PrimitivesType;
	typedef Gep::Evolution<PrimitivesType> EvolutionType;

	PrimitivesType primitives(inputs,gepOptions.genes,constants);
	EvolutionType evolution(primitives,seed,verbose);

	if (example < 2) {
		main1<Gep::Example1,EvolutionType>(evolution,gepOptions,total);
		return 0;
	} else if (example == 2) {
		main1<Gep::Example2,EvolutionType>(evolution,gepOptions,total);
		return 0;
	}

	main1<Gep::Example3,EvolutionType>(evolution,gepOptions,total);
}
