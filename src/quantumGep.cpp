/*
Copyright (c) 2017-2021, UT-Battelle, LLC

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
#include "Primitives/QuantumCircuit.h"
#include "Engine.h"
#include <unistd.h>
#include "Functions/QuantumOracle.h"

template<template<typename> class FitnessTemplate,
         typename EvolutionType>
void main1(EvolutionType& evolution,
           const Gep::Options& gepOptions,
           SizeType total)
{
	typedef FitnessTemplate<EvolutionType> FitnessType;
	typedef Gep::Engine<FitnessType> EngineType;
	typedef typename FitnessType::MinimizerParamsType MinimizerParamsType;
	typedef typename MinimizerParamsType::RealType RealType;

	// TODO FIXME Read from input file
	typename MinimizerParamsType::EnumAlgo algo = MinimizerParamsType::EnumAlgo::CONJUGATE_GRADIENT;
	SizeType maxIter = 100;
	SizeType saveEvery = 0;
	RealType delta = 0.01;
	RealType delta2 = 0.01;
	RealType tol = 1e-3;
	bool verbose = true;

	MinimizerParamsType minParams(algo, maxIter, delta, delta2, tol, saveEvery, verbose);
	typename EngineType::ParametersEngineType params(gepOptions);
	EngineType engine(params, evolution, &minParams);

	for (SizeType i = 0; i < total; i++)
		if (engine.evolve() && gepOptions.stopEarly) break;
}

int main(int argc, char* argv[])
{
	SizeType inputs = 0;
	SizeType total = 0;
	SizeType seed = 1000;
	SizeType numberOfBits = 0;
	bool verbose = false;
	Gep::Options gepOptions;

	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -i inputs -h head -b numberOfBits -p population -t total [-g genes -H chead]\n";
	while ((opt = getopt(argc, argv,"i:h:g:s:p:t:b:H:a:Sv")) != -1) {
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
		case 'b':
			numberOfBits = atoi(optarg);
			break;
		case 'H':
			gepOptions.chead = atoi(optarg);
			break;
		case 'a':
			gepOptions.adfs = atoi(optarg);
			break;
		case 'S':
			gepOptions.stopEarly = true;
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

	if (numberOfBits == 0)
		err("You need -b numberOfBits\n");

	if (gepOptions.genes > 1 && (gepOptions.chead == 0 || gepOptions.adfs == 0))
		throw PsimagLite::RuntimeError(strUsage);

	typedef std::complex<double> ComplexType;
	typedef PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef Gep::QuantumCircuit<VectorType> PrimitivesType;
	typedef Gep::Evolution<PrimitivesType> EvolutionType;

	PrimitivesType primitives(inputs, gepOptions.genes, numberOfBits);
	EvolutionType evolution(primitives,seed,verbose);

	main1<Gep::QuantumOracle,EvolutionType>(evolution,gepOptions,total);
}
