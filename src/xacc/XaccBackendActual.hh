#ifndef XACC_BACKEND_ACTUAL_HH
#define XACC_BACKEND_ACTUAL_HH

namespace Gep {

class XaccBackend {

public:

	XaccBackend(int argc, char* argv[])
	{
		xacc::Initialize(argc, argv);
	}

	~XaccBackend()
	{
		xacc::Finalize();
	}
};
}
#endif
