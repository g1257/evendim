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
#ifndef TREE_H
#define TREE_H
#include "PsimagLite.h"

namespace Gep {

template<typename PrimitivesType>
class Tree {

public:

	typedef Tree<PrimitivesType> TreeType;
	typedef typename PrimitivesType::ValueType ValueType;
	typedef typename PrimitivesType::NodeType NodeType;
	typedef typename PsimagLite::Vector<const TreeType*>::Type VectorTreeType;
	typedef typename NodeType::AnglesType AnglesType;
	typedef typename PsimagLite::Vector<AnglesType>::Type VectorAnglesType;
	typedef typename PsimagLite::Vector<ValueType>::Type VectorValueType;

	Tree(const PrimitivesType& primitives,const NodeType& node, bool verbose)
	    : primitives_(primitives),
	      node_(node),
	      verbose_(verbose)
	{}

	~Tree()
	{}

	ValueType exec() const
	{
		if (verbose_) std::cout<<" type= "<<node_.code()<<"\n";
		VectorValueType values(descendants_.size());

		for (SizeType i = 0; i < descendants_.size(); i++) {
			values[i] = descendants_[i]->exec();
		}

//		for (SizeType i = 0; i < values.size(); i++) {
//			if (values[i]<0 || values[i]==0 || values[i] >0) continue;
//			throw std::runtime_error("exec\n");
//		}

		ValueType tmp = node_.exec(values);
		if (verbose_) {
			std::cout<<"tmp= "<<tmp<<" type= "<<node_.code();
			for (SizeType i = 0; i < values.size(); i++)
				std::cout<<values[i]<<" ";
			std::cout<<"\n";
		}

		 return tmp;
	}

	void set(const VectorValueType& values) const
	{
		for (SizeType i = 0; i < descendants_.size(); i++)
			descendants_[i]->set(values);

		const PsimagLite::String str = node_.code();
		if (str.length() == 0)
			err("Node with empty code!?\n");

		if (str.length() > 1) return;
		const char c = str[0];
		if (c < 48 || c > 57) return;
		SizeType index = c - 48;
		assert(index < values.size());
		node_.set(values[index]);
	}

	void setDescendants(const TreeType& n0)
	{
		descendants_.push_back(&n0);
	}

	void setDescendants(const TreeType& n0,const TreeType& n1)
	{
		descendants_.push_back(&n0);
		descendants_.push_back(&n1);
	}

private:

	Tree& operator=(const Tree&) = delete;

	const PrimitivesType& primitives_;
	const NodeType& node_;
	bool verbose_;
	VectorTreeType descendants_;
};

} // namespace Gep

#endif // TREE_H
