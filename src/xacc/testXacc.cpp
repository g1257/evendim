#include "LinearTreeExecXacc.hh"
#include "XaccBackendActual.hh"
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

std::string implodeVecString(const std::vector<std::string>& vecStr)
{
	const char* const delim = ", ";

	std::ostringstream imploded;
	std::copy(vecStr.begin(), vecStr.end(), std::ostream_iterator<std::string>(imploded, delim));
	return imploded.str();
}

int main(int argc, char** argv)
{
	// This call xacc::init in its ctor and xacc:fin in its dtor
	Gep::XaccBackend xaccBackend(argc, argv);

	using LinearTreeExecType = Gep::LinearTreeExec<std::complex<double>>;

	typename LinearTreeExecType::VecStringType mycircuit { "Sx0" };

	constexpr int threadNum = 0; // no parallelization for now
	LinearTreeExecType linearTreeXacc(mycircuit, threadNum);

	// Does energy = <0|C H C |0>, with H = Z_0
	double energy = linearTreeXacc.energy();
	std::cout << "Circuit is " << implodeVecString(mycircuit) << "\n";
	std::cout << "Energy is " << energy << "\n";
}
