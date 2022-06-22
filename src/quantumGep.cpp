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
#include "Fitness/QuantumFitness.h"
#include "Fitness/GroundStateFitness.h"
#include "Fitness/Hamiltonian.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "FloatingPoint.h"

template<template<typename> class FitnessTemplate, typename EvolutionType>
void main2(EvolutionType& evolution,
           const Gep::ParametersEngine<double>& params,
           PsimagLite::InputNg<Gep::InputCheck>::Readable& io)
{
	typedef Gep::Engine<FitnessTemplate, EvolutionType> EngineType;
	typedef typename EngineType::FitnessType FitnessType;
	typedef typename FitnessType::FitnessParamsType FitnessParamsType;

	FitnessParamsType fitParams(io, params.threads);

	EngineType engine(params, evolution, &fitParams);

	for (SizeType i = 0; i < params.generations; i++)
		if (engine.evolve(i) && params.options.isSet("stopEarly")) break;
}

/* PSIDOC quantumGepMain
This driver program named quantumGep uses GEP to find a quantum circuit.
There are two usages: (i) the quantum circuit to be found implements
a function known only by some of its input and outputs, and (ii) the quantum
circuit to be found yields the ground state of a known Hamiltonian when applied
to an initial quantum state.
The primitives are under Primitives/QuantumCircuit.h, and
consist of one-bit and two-bit gates.

quantumGep takes one mandatory argument: -f filename, with the name of the input file.
It takes the following optional arguments.
\begin{itemize}
\item[-S] threads. The number of threads for shared memory parallelization.
\item[-p] precision. The precision for printing numbers.
\item[-v] indicates that quantumGep be verbose.
\end{itemize}
*/
int main(int argc, char* argv[])
{
	PsimagLite::String filename;
	SizeType threads = 0;
	bool verbose = false;

	PsimagLite::FloatingPoint::enableExcept();

	int opt = 0;
	int precision = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename [-S threads] [-p precision] [-v]\n";
	while ((opt = getopt(argc, argv,"f:S:p:v")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'v':
			verbose = true;
			break;
		case 'p':
			precision = atoi(optarg);
			break;
		case 'S':
			threads = PsimagLite::atoi(optarg);
			break;
		default:
			throw PsimagLite::RuntimeError(strUsage);
			return 1;
		}
	}

	if (filename == "")
		throw PsimagLite::RuntimeError(strUsage);

	if (precision > 0) {
		std::cout.precision(precision);
		std::cerr.precision(precision);
	}

	Gep::InputCheck inputCheck;
	PsimagLite::InputNg<Gep::InputCheck>::Writeable input(filename, inputCheck);
	PsimagLite::InputNg<Gep::InputCheck>::Readable io(input);

	Gep::ParametersInput gepOptions(io);

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

	PsimagLite::String runType;
	io.readline(runType, "RunType=");

	if (runType == "GroundState") gepOptions.samples = 1;

	typedef std::complex<double> ComplexType;
	typedef PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef Gep::QuantumCircuit<VectorType> PrimitivesType;
	typedef Gep::Evolution<PrimitivesType> EvolutionType;
	Gep::ParametersEngine<double> params(gepOptions);
	if (threads > 0) params.threads = threads;
	PsimagLite::CodeSectionParams codeSection(params.threads,
	                                          1, // threads2
	                                          false, // setAffinities,
	                                          0); // threadsStackSize;
	PsimagLite::Concurrency::setOptions(codeSection);

	PrimitivesType primitives(numberOfBits, gates, io);
	EvolutionType evolution(primitives, seed, verbose);

	if (runType == "FunctionFit") {
		main2<Gep::QuantumFitness, EvolutionType>(evolution, params, io);
	} else if (runType == "GroundState") {
		main2<Gep::GroundStateFitness, EvolutionType>(evolution, params, io);
	} else {
		err("RunType=FunctionFit or GroundState, but not " + runType + "\n");
	}
}
