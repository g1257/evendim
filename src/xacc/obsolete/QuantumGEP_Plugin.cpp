#include "QuantumGEP_Plugin.hh"

#include "PauliOperator.hpp"
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace xacc;
using namespace xacc::quantum;

namespace xacc {
namespace algorithm {
	bool QuantumGEP::initialize(const HeterogeneousMap& parameters)
	{
		/** Checks for the required parameters and other optional keywords
		 * @param[in] HeterogeneousMap A map of strings to keys
		 */

		std::cout << "QuantumGEP::initialize\n";

		if (!parameters.pointerLikeExists<Accelerator>("accelerator")) {
			xacc::error("Acc was false");
			return false;
		}

		if (!parameters.pointerLikeExists<Optimizer>("optimizer")) {
			xacc::error("Opt was false");
			return false;
		}

		optimizer = parameters.getPointerLike<Optimizer>("optimizer");
		accelerator = parameters.getPointerLike<Accelerator>("accelerator");

		return true;
	}

	const std::vector<std::string> QuantumGEP::requiredParameters() const
	{
		return { "accelerator" };
	}

	void QuantumGEP::execute(const std::shared_ptr<AcceleratorBuffer> buffer) const
	{

		std::cout << "execute\n";
		return;
	}

} // namespace algorithm
} // namespace xacc
