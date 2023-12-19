#include "xacc/ToPauliMatrices.hh"
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
	if (argc < 2)
		return -1;
	std::string name(argv[1]);

	Gep::ToPauliMatrices toPauliMatrices(name);

	std::cout << name << " expands to " << toPauliMatrices() << "\n";
}
