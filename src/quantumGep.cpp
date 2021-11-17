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
#include "InputNg.h"
#include "InputCheck.h"

template<template<typename> class FitnessTemplate,
         typename EvolutionType>
void main1(EvolutionType& evolution,
           const Gep::Options& gepOptions,
           PsimagLite::InputNg<Gep::InputCheck>::Readable& io)
{
	typedef FitnessTemplate<EvolutionType> FitnessType;
	typedef Gep::Engine<FitnessType> EngineType;
	typedef typename FitnessType::MinimizerParamsType MinimizerParamsType;

	SizeType total = 0;
	io.readline(total, "Generations=");
	if (total == 0)
		err("Generations must be greater than zero\n");

	MinimizerParamsType minParams(io);
	typename EngineType::ParametersEngineType params(gepOptions);
	EngineType engine(params, evolution, &minParams);

	for (SizeType i = 0; i < total; i++)
		if (engine.evolve() && gepOptions.stopEarly) break;
}

int main(int argc, char* argv[])
{
	PsimagLite::String filename;
	bool verbose = false;

	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename [-v]\n";
	while ((opt = getopt(argc, argv,"f:v")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'v':
			verbose = true;
			break;
		default:
			throw PsimagLite::RuntimeError(strUsage);
			return 1;
		}
	}

	if (filename == "")
		throw PsimagLite::RuntimeError(strUsage);

	Gep::InputCheck inputCheck;
	PsimagLite::InputNg<Gep::InputCheck>::Writeable input(filename, inputCheck);
	PsimagLite::InputNg<Gep::InputCheck>::Readable io(input);

	Gep::Options gepOptions(io);

	PsimagLite::String gates;
	if (gepOptions.primitives == "" || gepOptions.primitives == "?")
		gates = "C,H,P";
	else
		gates = gepOptions.primitives;

	if (gepOptions.primitives == "?") {
		std::cout<<"Default gates are: "<<gates<<"\n";
		return 0;
	}

	// sanity checks here
	if (gepOptions.head == 0 || gepOptions.population == 0) {
		throw PsimagLite::RuntimeError(strUsage);
		return 1;
	}

	SizeType numberOfBits = 0;
	io.readline(numberOfBits, "NumberOfBits=");

	if (numberOfBits == 0)
		err("You need numberOfBits > 0\n");

	SizeType seed = 12345;
	try {
		io.readline(seed, "RngSeed=");
	} catch (std::exception&) {}

	if (gepOptions.genes > 1 && (gepOptions.chead == 0 || gepOptions.adfs == 0))
		throw PsimagLite::RuntimeError(strUsage);

	typedef std::complex<double> ComplexType;
	typedef PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef Gep::QuantumCircuit<VectorType> PrimitivesType;
	typedef Gep::Evolution<PrimitivesType> EvolutionType;

	PrimitivesType primitives(numberOfBits, gates);
	EvolutionType evolution(primitives, seed, verbose);

	main1<Gep::QuantumOracle,EvolutionType>(evolution, gepOptions, io);
}
