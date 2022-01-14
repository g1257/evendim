#ifndef CANONICALFORMQUANTUM_H
#define CANONICALFORMQUANTUM_H
#include "Vector.h"
#include "Node.h"
#include "ProgramGlobals.h"
#include "Sort.h"
#include <cassert>
#include "PsimagLite.h"
#include <queue>

namespace Gep {

template<typename ValueType_, typename AnglesType>
class CanonicalFormQuantum {

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Vector<ValueType_>::Type VectorValueType;
	typedef Node<VectorValueType, AnglesType> NodeType;
	typedef typename PsimagLite::Vector<NodeType*>::Type VectorNodeType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::queue<PsimagLite::String> QueueStringType;
	typedef PsimagLite::Vector<QueueStringType>::Type VectorQueueStringType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;

	enum class RotationEnum {INVALID, X, Y, Z};

	CanonicalFormQuantum(const VectorStringType& data,
	                     const VectorNodeType& nodes)
	    : data_(data), needsChange_(false)
	{
		getEffectiveAndJunk();
		needsChange_ |= orderGatesByBit(effective_, nodes);
		needsChange_ |= compactifyRotations(effective_, nodes);
	}

	void changeIfNeeded(VectorStringType& vstr) const
	{
		if (!needsChange_) return;
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

			if (*it == "0") flag = false;
		}
	}

	static bool orderGatesByBit(VectorStringType& effective, const VectorNodeType& nodes)
	{
		constexpr SizeType maxConstant = 12345;
		static const ValueType_ value;
		constexpr bool isCell = false;
		VectorQueueStringType v(100);
		VectorSizeType order(effective.size(), maxConstant);
		VectorSizeType bits;
		bool flag = false;

		for (auto it = effective.begin(); it != effective.end(); ++it) {
			const NodeType& node = ProgramGlobals::findNodeFromCode<NodeType>(*it,
			                                                                  nodes,
			                                                                  value,
			                                                                  isCell);

			if (node.isInput()) continue;

			getBits(bits, node.code());
			if (bits.size() > 1) return false;// FIXME TODO: consider more than one bit gates

			assert(bits.size() == 1);
			const SizeType bit = bits[0];
			assert(bit < v.size());
			v[bit].push(*it);
			order[it - effective.begin()] = bit;
			flag = true;
		}

		if (!flag) return false;

		VectorSizeType iperm(order.size());
		PsimagLite::Sort<VectorSizeType> sort;
		sort.sort(order, iperm);

		VectorStringType newData(effective.size());

		const SizeType n = effective.size();
		for (SizeType i = 0; i < n; ++i) {
			const SizeType bit = order[i]; // order is ordered
			if (bit >= maxConstant) {
				newData[i] = effective[i];
				continue;
			}

			assert(bit < v.size());
			assert(v[bit].size() > 0);
			newData[i] = v[bit].front();

			v[bit].pop();
		}

		const bool b = isEqual(newData, effective);
		if (b) return false;

		effective = newData;
		return true;
	}

	static void getBits(VectorSizeType& bits, PsimagLite::String code)
	{
		VectorStringType tokens;
		PsimagLite::split(tokens, code, ":");
		if (tokens.size() > 2)
			err("getBit: code " + code + " with two or more colons\n");

		if (tokens.size() == 2) code = tokens[0]; // ignore angles

		bits.clear();

		while (true) {
			const SizeType n = code.size();
			const SizeType bp = getBreakpoint(code);
			if (bp == n) break;
			if (n  == bp + 1)
				err("getBits: Internal error\n");
			const PsimagLite::String buffer = code.substr(bp + 1, n - bp - 1);
			const SizeType bit = PsimagLite::atoi(buffer);
			bits.push_back(bit);
			if (code[bp] != '_') break;
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
			if (!std::isdigit(c)) break;
		}

		return n - ind - 1;
	}

	class Track {

	public:

		enum class StateEnum { COPY, IGNORE, NEW};

		typedef typename PsimagLite::Vector<RotationEnum>::Type VectorRotationEnumType;
		typedef typename PsimagLite::Vector<AnglesType>::Type VectorAnglesType;
		typedef typename PsimagLite::Vector<StateEnum>::Type VectorStateEnumType;

		Track(SizeType n)
		    : data_(n, StateEnum::COPY),
		      locations_(0),
		      dirs_(n, RotationEnum::INVALID),
		      angles_(n),
		      bits_(n)
		{}

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
			if (locations_ == 0) return;

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

	static bool compactifyRotations(VectorStringType& effective,
	                                const VectorNodeType& nodes)
	{
		static const ValueType_ value;
		constexpr bool isCell = false;
		RotationEnum prevDir = RotationEnum::INVALID;
		SizeType prevBit = 1e6;
		AnglesType prevAngle = 0;
		VectorSizeType bits;
		const SizeType n = effective.size();
		Track track(n);

		for (SizeType i = 0; i < n; ++i) {
			const NodeType& node = ProgramGlobals::findNodeFromCode<NodeType>(effective[i],
			                                                                  nodes,
			                                                                  value,
			                                                                  isCell);

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
			if (bits.size() != 1) err("Rotation gate must have one bit!?\n");

			assert(bits.size() == 1);
			const SizeType bit = bits[0];
			if (prevBit != bit) {

				track.compactify(prevDir, prevBit, prevAngle, i);

				prevBit = bit;
				prevAngle = angle;

				continue;
			}

			if (!mayCompactify) continue; // dir doesn't match previous

			track.increaseLocations();
			prevAngle += angle;
		}

		bool needsChange = false;
		VectorStringType newData;
		for (SizeType i = 0; i < n; ++i) {
			typename Track::StateEnum state = track.state(i);
			switch (state) {
			case Track::StateEnum::COPY:
				newData.push_back(effective[i]);
				break;
			case Track::StateEnum::IGNORE:
				break;
			case Track::StateEnum::NEW:
				PsimagLite::String newCode = track.code(i);
				newData.push_back(newCode);
				needsChange = true;
			}
		}

		if (!needsChange) return false;

		effective.swap(newData);
		return true;
	}

	static RotationEnum getRotationDirection(PsimagLite::String code)
	{
		if (code.length() < 2)
			err("getRotationDirection: code " + code + "\n");
		if (code[0] != 'R') return RotationEnum::INVALID;

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
		if (tokens.size() == 2) angle = tokens[1]; // ignore gate name and bit
		return PsimagLite::atof(angle);
	}

	static bool isEqual(const VectorStringType& v1, const VectorStringType& v2)
	{
		const SizeType n = v1.size();
		if (n != v2.size()) return false;
		for (SizeType i = 0; i < n; ++i)
			if (v1[i] != v2[i]) return false;
		return true;
	}

	const VectorStringType& data_;
	bool needsChange_;
	VectorStringType effective_;
	VectorStringType junkDna_;
};
}
#endif // CANONICALFORMQUANTUM_H
