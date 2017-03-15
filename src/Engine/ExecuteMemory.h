#ifndef EXECUTEMEMORY_H
#define EXECUTEMEMORY_H

#include "Vector.h"
#include <sys/mman.h>
#include <unistd.h>
#include <stdio.h>

class ExecuteMemory {

public:

	ExecuteMemory(SizeType length)
	    : length_(length),ptr_(0)
	{
		checkLength();
		int prot = PROT_READ | PROT_WRITE | PROT_EXEC;
		ptr_ = (char *)mmap(0,length,prot,MAP_PRIVATE | MAP_ANONYMOUS,-1,0);
		if ((void*)ptr_ == MAP_FAILED) {
			perror("mmap");
			throw PsimagLite::RuntimeError("ExecuteMemory::ctor()\n");
		}
	}

	~ExecuteMemory()
	{
		int x = munmap(ptr_,length_);
		if (x == -1) {
			perror("munmap");
		}
	}

	template<typename FunctionType>
	int operator()(const PsimagLite::String& str,int x, int y) const
	{
		if (str.length()>=length_) {
			throw PsimagLite::RuntimeError("ExecuteMemory: string too long\n");
		}

		for (SizeType i = 0; i < str.length(); i++) {
			ptr_[i] = str[i];
		}
		FunctionType f = (FunctionType) ptr_;
		return f(x,y);
	}

private:

	void checkLength() const
	{
		long sz = sysconf(_SC_PAGESIZE);
		if (length_ % sz !=0)
			throw PsimagLite::RuntimeError("ExecuteMemory::checkLength()\n");
	}

	SizeType length_;
	char *ptr_;

}; // class ExecuteMemory

#endif // EXECUTEMEMORY_H
