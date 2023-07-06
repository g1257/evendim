#ifndef QUASIVECTOR_HH
#define QUASIVECTOR_HH
#include "ProgramGlobals.h"
#include "Vector.h"
#include <string>

namespace Gep
{

template <typename ComplexOrRealType>
class QuasiVector
{

public:

	using VectorType =
	    typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using RealType =
	    typename PsimagLite::Real<ComplexOrRealType>::Type;
	using value_type = ComplexOrRealType;

	QuasiVector()
	    : size_(0)
	    , isExp_(false)
	{
	}

	QuasiVector(SizeType size)
	    : size_(size)
	    , isExp_(false)
	{
	}

	QuasiVector(const std::string& filename) { fromFile(filename); }

	void fromFile(const std::string& filename)
	{
		isExp_ = true;
		Gep::ProgramGlobals::readVector(data_, filename);
		size_ = data_.size();
	}

	template <typename SomeRngType>
	void randomize(SizeType size, SomeRngType& rng, const ComplexOrRealType& a, const ComplexOrRealType& b)
	{
		blowUp(size);
		needsExp("randomize");
		ProgramGlobals::randomVector(data_, rng, a, b);
	}

	void flipABit(const QuasiVector& src, SizeType bit)
	{
		assert(size_ == src.size());
		SizeType mask = (1 << bit);
		for (SizeType i = 0; i < size_; ++i) {
			SizeType j = i ^ mask;
			data_[j] = src.data_[i];
		}
	}

	// PUBLIC CONST FUNCTIONS BELOW

	// Use toVector() only for IsingGraph
	const VectorType& toVector() const
	{
		// cop out for now; remove later
		needsExp("toVector");
		return data_;
	}

	bool hasWeight(SizeType ind, const RealType& epsilon) const
	{
		assert(ind < data_.size());
		return (std::norm(data_[ind]) > epsilon);
	}

	SizeType size() const { return size_; }

	RealType norm() const
	{
		needsExp("norm");
		return PsimagLite::norm(data_);
	}

	void print(std::ostream& os) const
	{
		needsExp("print");
		ProgramGlobals::writeVector(os, data_);
	}

	friend std::ostream& operator<<(std::ostream& os,
	    const QuasiVector& qv)
	{
		qv.needsExp("operator<<");
		for (SizeType i = 0; i < qv.data_.size(); ++i) {
			os << qv.data_[i] << " ";
		}

		return os;
	}

	friend RealType vectorDiff2(const QuasiVector& v1,
	    const QuasiVector& v2)
	{
		return vectorDiff2_(v1.toVector(), v2.toVector());
	}

	friend RealType diffVectorDiff2(const QuasiVector& v1,
	    const QuasiVector& v2,
	    const QuasiVector& v3)
	{
		return diffVectorDiff2_(v1.toVector(), v2.toVector(), v3.toVector());
	}

	friend QuasiVector
	oneBitGate(const QuasiVector& src, SizeType bit, const PsimagLite::Matrix<ComplexOrRealType>& gate)
	{
		SizeType n = src.size();
		QuasiVector w(n);

		w.blowUp(n);
		for (SizeType i = 0; i < n; ++i) {
			SizeType j = findBasisState(i, bit);
			SizeType bitI = getBitForIndex(i, bit);
			SizeType bitJ = getBitForIndex(j, bit);
			w.data_[i] += gate(bitI, bitI) * src.data_[i];
			w.data_[j] += gate(bitI, bitJ) * src.data_[i];
		}

		return w;
	}

	friend QuasiVector CNOT(const QuasiVector& src, SizeType bit1, SizeType bit2)
	{
		const int n = src.size(); // 2^Nbits

		QuasiVector w(n);
		w.blowUp(n);
		const SizeType mask2 = (1 << bit2);
		for (int i = 0; i < n; ++i) {
			const SizeType oldContent1 = getBitForIndex(i, bit1);
			assert(oldContent1 < 2);
			const SizeType oldContent2 = getBitForIndex(i, bit2);
			assert(oldContent2 < 2);
			const SizeType content2 = (oldContent1 + oldContent2) % 2;
			assert(content2 < 2);

			const SizeType j = (content2 == oldContent2) ? i : (i ^ mask2);

			w.data_[j] += src.data_[i];
		}

		return w;
	}

	// <v1|H|v2>
	// caching has been disabled here!
	template <typename SomeMatrixType>
	friend RealType tensorEnergy(const QuasiVector& v1,
	    const SomeMatrixType& H,
	    const QuasiVector& v2)
	{
		assert(v1.size() == v2.size());
		assert(H.cols() == v2.size());
		assert(H.rows() == v1.size());
		VectorType tmpVector(v1.size());
		H.matrixVectorProduct(tmpVector, v2.toVector());
		return PsimagLite::real(v1.toVector() * tmpVector);
	}

private:

	void blowUp(SizeType size)
	{
		data_.resize(size);
		isExp_ = true;
		size_ = size;
	}

	void needsExp(const std::string& info) const
	{
		if (isExp_)
			return;
		err(info + " unimplemented or non-working unless exponential "
			   "representation\n");
	}

	static RealType vectorDiff2_(const VectorType& v1,
	    const VectorType& v2)
	{
		const SizeType n = v1.size();
		assert(n == v2.size());
		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i)
			sum += std::abs(v2[i] - v1[i]);

		return sum / n;
	}

	static RealType diffVectorDiff2_(const VectorType& v1,
	    const VectorType& v2,
	    const VectorType& v3)
	{
		const SizeType n = v1.size();
		assert(n == v2.size());
		assert(n == v3.size());
		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i) {
			RealType denom = std::abs(v2[i] - v1[i]);
			if (denom == 0)
				denom = 1;
			const RealType re1 = PsimagLite::real(v2[i] - v1[i]);
			const RealType im1 = PsimagLite::imag(v2[i] - v1[i]);
			sum += (PsimagLite::real(v3[i]) * re1 + PsimagLite::imag(v3[i]) * im1) / denom;
		}

		return sum / n;
	}

	static SizeType findBasisState(SizeType ind, SizeType bit)
	{
		const SizeType mask = (1 << bit);
		return ind ^ mask;
	}

	static SizeType getBitForIndex(SizeType ind, SizeType bitNumber)
	{
		const SizeType mask = (1 << bitNumber);
		const SizeType result = ind & mask;
		return (result > 0) ? 1 : 0;
	}

	SizeType size_;
	bool isExp_;
	VectorType data_;
};
} // namespace Gep
#endif // QUASIVECTOR_HH
