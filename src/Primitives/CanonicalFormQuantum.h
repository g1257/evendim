#ifndef CANONICALFORMQUANTUM_H
#define CANONICALFORMQUANTUM_H
#include "AST/Node.h"
#include "NodeFactory.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"
#include "Sort.h"
#include "Vector.h"
#include <cassert>
#include <queue>

namespace Gep {

template <typename ValueType_, typename AnglesType>
class CanonicalFormQuantum {

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Vector<ValueType_>::Type VectorValueType;
	typedef PsimagLite::Node<VectorValueType, AnglesType> NodeType;
	typedef typename PsimagLite::Vector<NodeType*>::Type VectorNodeType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::queue<PsimagLite::String> QueueStringType;
	typedef PsimagLite::Vector<QueueStringType>::Type VectorQueueStringType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef NodeFactory<NodeType> NodeFactorType;

	enum class RotationEnum { INVALID,
		                  X,
		                  Y,
		                  Z };

	CanonicalFormQuantum(const VectorStringType& data,
	                     const NodeFactorType& nodeFactory)
	    : data_(data)
	    , nodeFactory_(nodeFactory)
	    , needsChange_(false)
	{
		getEffectiveAndJunk();
		needsChange_ |= orderGatesByBitNg(effective_);
		needsChange_ |= compactifyRotations(effective_, junkDna_);
	}

	void changeIfNeeded(VectorStringType& vstr) const
	{
		if (!needsChange_)
			return;
		vstr.resize(effective_.size() + junkDna_.size());
		for (SizeType i = 0; i < effective_.size(); ++i)
			vstr[i] = effective_[i];

		for (SizeType i = 0; i < junkDna_.size(); ++i)
			vstr[i + effective_.size()] = junkDna_[i];
	}

private:

	void getEffectiveAndJunk()
	{
		bool flag = true;
		for (auto it = data_.begin(); it != data_.end(); ++it) {
			if (flag)
				effective_.push_back(*it);
			else
				junkDna_.push_back(*it);

			if (*it == "0")
				flag = false;
		}
	}

	bool orderGatesByBitNg(VectorStringType& effective) const
	{
		bool flag = false;
		while (orderGatesByBitOneRound(effective)) {
			flag = true;
		}

		return flag;
	}

	bool orderGatesByBitOneRound(VectorStringType& effective) const
	{
		bool flag = false;
		const SizeType n = effective.size();
		for (SizeType i = 1; i < n; ++i) {
			bool gateMoved = moveThisGateIfPossible(effective, i);
			flag |= gateMoved;
		}

		return flag;
	}

	bool moveThisGateIfPossible(VectorStringType& effective, SizeType ind) const
	{
		if (ind == 0)
			return false;

		static const ValueType_ value;
		constexpr bool isCell = false;
		constexpr SizeType threadNum = 0;

		const NodeType& node = nodeFactory_.findNodeFromCode(effective[ind],
		                                                     value,
		                                                     isCell,
		                                                     threadNum);
		if (node.isInput())
			return false;
		VectorSizeType bits;
		getBits(bits, node.code());

		int jnd = ind - 1;
		int location = -1;
		for (; jnd >= 0; --jnd) {
			const NodeType& nodePrev = nodeFactory_.findNodeFromCode(effective[jnd],
			                                                         value,
			                                                         isCell,
			                                                         threadNum);
			if (nodePrev.isInput())
				err("moveThisGateIfPossible: input found before gate!?\n");
			VectorSizeType bitsPrev;

			getBits(bitsPrev, nodePrev.code());

			if (!isBefore(bitsPrev, bits))
				break;
			location = jnd;
		}

		if (location < 0)
			return false;

		std::swap(effective[ind], effective[location]);
		return true;
	}

	static bool isBefore(const VectorSizeType& bits1, const VectorSizeType& bits2)
	{
		const SizeType n1 = bits1.size();
		const SizeType n2 = bits2.size();
		for (SizeType i = 0; i < n2; ++i)
			for (SizeType j = 0; j < n1; ++j)
				if (bits1[j] <= bits2[i])
					return false;

		return true;
	}

	static void getBits(VectorSizeType& bits, PsimagLite::String code)
	{
		VectorStringType tokens;
		PsimagLite::split(tokens, code, ":");
		if (tokens.size() > 2)
			err("getBit: code " + code + " with two or more colons\n");

		if (tokens.size() == 2)
			code = tokens[0]; // ignore angles

		bits.clear();

		while (true) {
			const SizeType n = code.size();
			const SizeType bp = getBreakpoint(code);
			if (bp == n)
				break;
			if (n == bp + 1)
				err("getBits: Internal error\n");
			const PsimagLite::String buffer = code.substr(bp + 1, n - bp - 1);
			const SizeType bit = PsimagLite::atoi(buffer);
			bits.push_back(bit);
			if (code[bp] != '_')
				break;
			code = code.substr(0, bp);
		}
	}

	static SizeType getBreakpoint(PsimagLite::String code)
	{
		const SizeType n = code.size();
		SizeType ind = 0;
		for (; ind < n; ++ind) {
			const SizeType j = n - ind - 1;
			unsigned char c = code[j];
			if (!std::isdigit(c))
				break;
		}

		return n - ind - 1;
	}

	class Track {

	public:

		enum class StateEnum { COPY,
			               IGNORE,
			               NEW };

		typedef typename PsimagLite::Vector<RotationEnum>::Type VectorRotationEnumType;
		typedef typename PsimagLite::Vector<AnglesType>::Type VectorAnglesType;
		typedef typename PsimagLite::Vector<StateEnum>::Type VectorStateEnumType;

		Track(SizeType n)
		    : data_(n, StateEnum::COPY)
		    , locations_(0)
		    , dirs_(n, RotationEnum::INVALID)
		    , angles_(n)
		    , bits_(n)
		{
		}

		StateEnum state(SizeType i) const
		{
			assert(i < data_.size());
			return data_[i];
		}

		void increaseLocations()
		{
			++locations_;
		}

		void compactify(RotationEnum prevDir,
		                SizeType bit,
		                AnglesType prevAngle,
		                SizeType ind)
		{
			if (locations_ == 0)
				return;

			const SizeType start = ind - locations_ - 1;
			assert(start >= 0 && ind > start);
			for (SizeType i = start; i < ind; ++i)
				data_[i] = StateEnum::IGNORE;

			assert(ind > 0);
			data_[ind - 1] = StateEnum::NEW;

			dirs_[ind - 1] = prevDir;
			bits_[ind - 1] = bit;
			angles_[ind - 1] = prevAngle;

			locations_ = 0;
		}

		PsimagLite::String code(SizeType ind)
		{
			assert(ind < dirs_.size());

			RotationEnum dir = dirs_[ind];
			PsimagLite::String str = "R";
			switch (dir) {
			case RotationEnum::X:
				str += "x";
				break;
			case RotationEnum::Y:
				str += "y";
				break;
			case RotationEnum::Z:
				str += "z";
				break;
			default:
				err("Canonicalization: Invalid direction for rotation\n");
			}

			assert(bits_.size() > ind);
			SizeType bit = bits_[ind];
			str += ttos(bit);

			str += ":";

			assert(angles_.size() > ind);
			AnglesType angle = angles_[ind];
			str += ttos(angle);
			return str;
		}

	private:

		VectorStateEnumType data_;
		SizeType locations_;
		VectorRotationEnumType dirs_;
		VectorAnglesType angles_;
		VectorSizeType bits_;
	};

	bool compactifyRotations(VectorStringType& effective, VectorStringType& junkDna) const
	{
		constexpr SizeType threadNum = 0;
		static const ValueType_ value;
		constexpr bool isCell = false;
		RotationEnum prevDir = RotationEnum::INVALID;
		SizeType prevBit = 1e6;
		AnglesType prevAngle = 0;
		VectorSizeType bits;
		const SizeType n = effective.size();
		Track track(n);

		for (SizeType i = 0; i < n; ++i) {
			const NodeType& node = nodeFactory_.findNodeFromCode(effective[i],
			                                                     value,
			                                                     isCell,
			                                                     threadNum);

			if (node.isInput()) {
				track.compactify(prevDir, prevBit, prevAngle, i);
				prevDir = RotationEnum::INVALID;
				continue;
			}

			RotationEnum dir = getRotationDirection(node.code());
			if (dir == RotationEnum::INVALID) {
				track.compactify(prevDir, prevBit, prevAngle, i);
				prevDir = RotationEnum::INVALID;
				continue; // not a rotation gate
			}

			AnglesType angle = getAngle(node.code());

			bool mayCompactify = true;
			if (prevDir != dir) {

				track.compactify(prevDir, prevBit, prevAngle, i);

				mayCompactify = false;
				prevDir = dir;
				prevAngle = angle;
			}

			getBits(bits, node.code());
			if (bits.size() != 1)
				err("Rotation gate must have one bit!?\n");

			assert(bits.size() == 1);
			const SizeType bit = bits[0];
			if (prevBit != bit) {

				track.compactify(prevDir, prevBit, prevAngle, i);

				prevBit = bit;
				prevAngle = angle;

				continue;
			}

			if (!mayCompactify)
				continue; // dir doesn't match previous

			track.increaseLocations();
			prevAngle += angle;
		}

		bool needsChange = false;
		VectorStringType newData;
		SizeType ignored = 0;
		for (SizeType i = 0; i < n; ++i) {
			typename Track::StateEnum state = track.state(i);
			switch (state) {
			case Track::StateEnum::COPY:
				newData.push_back(effective[i]);
				break;
			case Track::StateEnum::IGNORE:
				++ignored;
				needsChange = true;
				break;
			case Track::StateEnum::NEW:
				PsimagLite::String newCode = track.code(i);
				newData.push_back(newCode);
				needsChange = true;
			}
		}

		if (!needsChange)
			return false;

		// keep total size constant
		for (SizeType i = 0; i < ignored; ++i)
			junkDna.push_back("0");

		effective.swap(newData);
		return true;
	}

	static RotationEnum getRotationDirection(PsimagLite::String code)
	{
		if (code.length() < 2)
			err("getRotationDirection: code " + code + "\n");
		if (code[0] != 'R')
			return RotationEnum::INVALID;

		switch (code[1]) {
		case 'x':
			return RotationEnum::X;
		case 'y':
			return RotationEnum::Y;
		case 'z':
			return RotationEnum::Z;
		default:
			return RotationEnum::INVALID;
		}
	}

	static AnglesType getAngle(PsimagLite::String code)
	{
		VectorStringType tokens;
		PsimagLite::split(tokens, code, ":");
		if (tokens.size() > 2)
			err("getBit: code " + code + " with two or more colons\n");

		PsimagLite::String angle = "0";
		if (tokens.size() == 2)
			angle = tokens[1]; // ignore gate name and bit
		return PsimagLite::atof(angle);
	}

	static bool isEqual(const VectorStringType& v1, const VectorStringType& v2)
	{
		const SizeType n = v1.size();
		if (n != v2.size())
			return false;
		for (SizeType i = 0; i < n; ++i)
			if (v1[i] != v2[i])
				return false;
		return true;
	}

	const VectorStringType& data_;
	const NodeFactorType& nodeFactory_;
	bool needsChange_;
	VectorStringType effective_;
	VectorStringType junkDna_;
};
}
#endif // CANONICALFORMQUANTUM_H
