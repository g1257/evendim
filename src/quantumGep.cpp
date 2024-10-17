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
#include "Engine.h"
#include "Evolution.h"
#include "Fitness/GroundStateFitness.h"
#include "Fitness/QuantumFitness.h"
#include "FloatingPoint.h"
#include "InputCheck.h"
#include "InputNg.h"
#include "MpiShim.hh"
#include "Primitives/QuantumCircuit.h"
#include "Primitives/QuasiVector.hh"
#include "RedirectOutput.hh"
#include "XaccBackend.hh"
#include <unistd.h>

template <template <typename> class FitnessTemplate, typename EvolutionType>
void main2(EvolutionType& evolution,
           const Gep::ParametersEngine<double>& params,
           PsimagLite::InputNg<Gep::InputCheck>::Readable& io)
{
	typedef Gep::Engine<FitnessTemplate, EvolutionType> EngineType;
	typedef typename EngineType::FitnessType FitnessType;
	typedef typename FitnessType::FitnessParamsType FitnessParamsType;

	FitnessParamsType fitParams(io, params.threads);

	EngineType engine(params, evolution, &fitParams);

	for (SizeType i = 0; i < params.generations; i++) {
		engine.evolve(i);
	}
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
	Gep::MpiShim mpi_shim(argc, argv);
	PsimagLite::String filename;
	SizeType threads = 0;
	bool verbose = false;

	PsimagLite::FloatingPoint::enableExcept();

	int opt = 0;
	int precision = 0;
	PsimagLite::String label;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename [-S threads] [-p precision] [-v]\n";
	while ((opt = getopt(argc, argv, "f:S:p:l:v")) != -1) {
		switch (opt) {
		case 'f':
			filename = mpi_shim.buildInput(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		case 'p':
			precision = atoi(optarg);
			break;
		case 'l':
			label = optarg;
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

	PsimagLite::RedirectOutput::setAppName(argv[0]);

	if (!label.empty()) {
		if (mpi_shim.isMPI()) {
			err("Cannot use -l command line option with MPI\n");
		}
	}

	if (label != "-") {
		constexpr bool unbuffered = true;
		std::string output = (label.empty()) ? mpi_shim.buildOutput(filename) : label;
		PsimagLite::RedirectOutput::doIt(output, std::ofstream::out, unbuffered);
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
		std::cout << "Default gates are: " << gates << "\n";
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
	}
	catch (std::exception&) {
	}

	if (gepOptions.chead > 0 && gepOptions.adfs == 0)
		throw PsimagLite::RuntimeError("FATAL: You selected ADF head size H > 0 but ADF number a == 0\n");

	if (gepOptions.chead > 0 && gepOptions.adfs == 0)
		throw PsimagLite::RuntimeError("FATAL: You selected ADF number a > 0 but ADF head size H == 0\n");

	bool hasAdfs = (gepOptions.chead > 0 && gepOptions.adfs > 0);
	if (gepOptions.genes > 1 && !hasAdfs)
		throw PsimagLite::RuntimeError("FATAL: genes > 1 but no ADF\n");

	PsimagLite::String runType;
	io.readline(runType, "RunType=");

	if (runType == "GroundState")
		gepOptions.samples = 1;

	typedef std::complex<double> ComplexType;
	typedef Gep::QuasiVector<ComplexType> QuasiVectorType;
	typedef Gep::QuantumCircuit<QuasiVectorType> PrimitivesType;
	typedef Gep::Evolution<PrimitivesType> EvolutionType;
	Gep::ParametersEngine<double> params(gepOptions);
	if (threads > 0)
		params.threads = threads;
	PsimagLite::CodeSectionParams codeSection(params.threads,
	                                          1, // threads2
	                                          false, // setAffinities,
	                                          0); // threadsStackSize;
	PsimagLite::Concurrency::setOptions(codeSection);

	// Xacc backend if needed
	Gep::XaccBackend xaccBackend(argc, argv);

	PrimitivesType primitives(numberOfBits, gates, io);
	EvolutionType evolution(primitives, seed, verbose);

	if (runType == "FunctionFit") {
		main2<Gep::QuantumFitness, EvolutionType>(evolution, params, io);
	}
	else if (runType == "GroundState") {
		main2<Gep::GroundStateFitness, EvolutionType>(evolution, params, io);
	}
	else {
		err("RunType=FunctionFit or GroundState, but not " + runType + "\n");
	}
}
