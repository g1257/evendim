#ifndef EVENDIM_MPI_SHIM_SERIAL_HH
#define EVENDIM_MPI_SHIM_SERIAL_HH
#include "PsimagLite.h"
#include <map>

namespace Gep {

class MpiShim {

public:

	MpiShim(int argc, char** argv)
	{
	}

	unsigned int rank() const { return 0; }

	unsigned int size() const { return 1; }

	bool isMPI() const { return false; }

	static std::string buildInput(const std::string& name)
	{
		bool is_open = false;
		std::string buffer;
		std::string var;

		for (char c : name) {
			if (c == ';' && is_open) {
				is_open = false;
				buffer += getVarValue(var);
				var = "";
				continue;
			}

			if (c == '%' && !is_open) {
				is_open = true;
				var = "";
				continue;
			}

			if (is_open) {
				var += c;
				continue;
			}

			buffer += c;
		}

		return buffer;
	}

	static PsimagLite::String buildOutput(PsimagLite::String filename)
	{
		PsimagLite::String rootname = PsimagLite::basename(filename);
		size_t index = rootname.find(".", 0);
		if (index != PsimagLite::String::npos) {
			rootname.erase(index, filename.length());
		}

		return "runFor" + rootname + ".cout";
	}

private:

	static std::string getVarValue(const std::string& name)
	{
		if (name == "i") {
			return ttos(0);
		}
		else if (name == "n") {
			return ttos(1);
		}
		else {
			return "";
		}
	}
};

}

#endif
