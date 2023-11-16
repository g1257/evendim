#ifndef HAMILTONIANSPEC_H
#define HAMILTONIANSPEC_H
#include "../Primitives/QuantumOneBitGate.h"
#include "AuxForHamSpec.h"
#include "Vector.h"

namespace Gep {

template <typename SparseMatrixType>
class HamiltonianSpec {

public:

	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef AuxForHamSpec AuxiliaryType;
	typedef OneBitGateLibrary<ComplexOrRealType> OneBitGateLibraryType;
	typedef typename OneBitGateLibraryType::MatrixType MatrixType;

	class MyResult {

	public:

		MyResult(SizeType numberOfBits)
		    : bits_(numberOfBits)
		    , data_(1 << numberOfBits, 1 << numberOfBits)
		{
		}

		MyResult(PsimagLite::String str, SizeType numberOfBits)
		    : bits_(numberOfBits)
		    , data_(1 << numberOfBits, 1 << numberOfBits)
		{
			MatrixType gateMatrix;
			std::pair<PsimagLite::String, SizeType> nameBitPair = extractNameAndBit(str);
			if (nameBitPair.second >= numberOfBits)
				err("Bit too large in expression " + ttos(nameBitPair.second));

			OneBitGateLibraryType::fillAnyGate(gateMatrix, nameBitPair.first);
			blowUp(gateMatrix, nameBitPair.second);
		}

		MyResult& operator+=(const MyResult& other)
		{
			data_ += other.data_;
			return *this;
		}

		MyResult& operator*=(const MyResult& other)
		{
			SparseMatrixType r;
			multiply(r, data_, other.data_); // order FIXME TODO
			data_ = r;
			return *this;
		}

		MyResult& operator*=(const ComplexOrRealType& scalar)
		{
			SparseMatrixType r;
			r = scalar * data_;
			data_ = r;
			return *this;
		}

		const SparseMatrixType& getCRS() const { return data_; }

	private:

		void blowUp(const MatrixType& gateMatrix, SizeType bit)
		{
			const SizeType smallCols = gateMatrix.cols();
			assert(smallCols == gateMatrix.rows());
			const SizeType rows = data_.rows();
			assert(rows == data_.cols());

			SizeType counter = 0;
			for (SizeType i = 0; i < rows; ++i) {
				data_.setRow(i, counter);
				const SizeType ip = extractIndexAtBit(i, bit);
				for (SizeType jp = 0; jp < smallCols; ++jp) {
					const ComplexOrRealType val = gateMatrix(ip, jp);
					if (std::norm(val) == 0)
						continue;
					SizeType j = replaceIndexAtBit(i, bit, jp);
					data_.pushValue(val);
					data_.pushCol(j);
					++counter;
				}
			}

			data_.setRow(rows, counter);
		}

		// Ry0:1.57
		static std::pair<PsimagLite::String, SizeType> extractNameAndBit(PsimagLite::String str)
		{
			const SizeType n = str.length();
			PsimagLite::String bufferName;
			PsimagLite::String bufferBit;

			SizeType ind = 0;
			for (; ind < n; ++ind) {
				const char c = str[ind];
				if (isADigit(c)) {
					bufferBit += c;
					++ind;
					break;
				}
				else {
					bufferName += c;
				}
			}

			bool hasAngle = false;
			for (; ind < n; ++ind) {
				const char c = str[ind];
				if (isADigit(c)) {
					bufferBit += c;
				}
				else {
					if (c != ':')
						err("Expected digit or : in " + str + "\n");
					hasAngle = true;
					break;
				}
			}

			PsimagLite::String angle = (hasAngle) ? str.substr(ind, str.length() - ind)
			                                      : "";
			return std::pair<PsimagLite::String, SizeType>(bufferName + angle,
			                                               PsimagLite::atoi(bufferBit));
		}

		static bool isADigit(const char c)
		{
			return (c >= 48 && c <= 57);
		}

		SizeType extractIndexAtBit(SizeType ind, SizeType bit) const
		{
			checkBits(ind, bit);
			SizeType mask = (1 << bit);
			return (mask & ind) ? 1 : 0;
		}

		SizeType replaceIndexAtBit(SizeType ind, SizeType bit, SizeType jp) const
		{
			checkBits(ind, bit);
			assert(jp < 2);
			SizeType kp = extractIndexAtBit(ind, bit);
			if (kp == jp)
				return ind; // early exit here

			SizeType mask = (1 << bit);
			return mask ^ ind;
		}

		void checkBits(SizeType ind, SizeType bit) const
		{
			if (ind < data_.rows())
				return;
			if (bit < bits_)
				return;
			err("HamiltonianSpec: Bits out of range in expression\n");
		}

		SizeType bits_;
		SparseMatrixType data_;
	};

	typedef MyResult ResultType;

	HamiltonianSpec(SizeType numberOfBits)
	    : bits_(numberOfBits)
	    , fullMatrix_(1 << bits_, 1 << bits_)
	{
	}

	ResultType operator()(PsimagLite::String str, AuxiliaryType& aux) const
	{
		return ResultType(str, aux.numberOfBits());
	}

	static bool isEmpty(const MyResult& quasiMatrix)
	{
		return isZero(quasiMatrix.getCRS());
	}

	static bool metaEqual(const MyResult&, const MyResult&)
	{
		return true;
	}

private:

	SizeType bits_;
	SparseMatrixType fullMatrix_;
};
}
#endif // HAMILTONIANSPEC_H
