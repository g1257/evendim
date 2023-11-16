/*
Copyright (c) 2017, UT-Battelle, LLC

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
#ifndef EXECUTEMEMORY_H
#define EXECUTEMEMORY_H

#include "Vector.h"
#include <stdio.h>
#include <sys/mman.h>
#include <unistd.h>

class ExecuteMemory {

public:

	ExecuteMemory(SizeType length)
	    : length_(length)
	    , ptr_(0)
	{
		checkLength();
		int prot = PROT_READ | PROT_WRITE | PROT_EXEC;
		ptr_ = (char*)mmap(0, length, prot, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
		if ((void*)ptr_ == MAP_FAILED) {
			perror("mmap");
			throw PsimagLite::RuntimeError("ExecuteMemory::ctor()\n");
		}
	}

	~ExecuteMemory()
	{
		int x = munmap(ptr_, length_);
		if (x == -1) {
			perror("munmap");
		}
	}

	template <typename FunctionType>
	int operator()(const PsimagLite::String& str, int x, int y) const
	{
		if (str.length() >= length_) {
			throw PsimagLite::RuntimeError("ExecuteMemory: string too long\n");
		}

		for (SizeType i = 0; i < str.length(); i++) {
			ptr_[i] = str[i];
		}
		FunctionType f = (FunctionType)ptr_;
		return f(x, y);
	}

private:

	void checkLength() const
	{
		long sz = sysconf(_SC_PAGESIZE);
		if (length_ % sz != 0)
			throw PsimagLite::RuntimeError("ExecuteMemory::checkLength()\n");
	}

	SizeType length_;
	char* ptr_;

}; // class ExecuteMemory

#endif // EXECUTEMEMORY_H
