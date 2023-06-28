#ifndef QUASIVECTOR_HH
#define QUASIVECTOR_HH
#include "Vector.h"
#include "ProgramGlobals.h"
#include <string>

namespace Gep {

template<typename ComplexOrRealType>
class QuasiVector {

public:

    using VectorType = typename PsimagLite::Vector<ComplexOrRealType>::Type;
    using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;
    using value_type = ComplexOrRealType;

    QuasiVector() : size_(0), isExp_(false) {}

    QuasiVector(SizeType size) : size_(size), isExp_(false) {}

    QuasiVector(const std::string& filename)
    {
        fromFile(filename);
    }

    void fromFile(const std::string& filename)
    {
        isExp_ = true;
        Gep::ProgramGlobals::readVector(data_, filename);
        size_ = data_.size();
    }

    template<typename SomeRngType>
    void randomize(SomeRngType& rng)
    {
        needsExp("randomize");
        ProgramGlobals::randomVector(data_, rng);
    }

    void blowUp(SizeType size)
    {
        data_.resize(size);
        isExp_ = true;
        size_ = size;
    }

    void setTo(const ComplexOrRealType& val)
    {
        needsExp("setTo");
        std::fill(data_.begin(), data_.end(), val);
    }

    ComplexOrRealType& operator[](SizeType ind)
    {
        // needsExp obviously, but disabled for performance here
        assert(ind < data_.size());
        return data_[ind];
    }

    // PUBLIC CONST FUNCTIONS BELOW

    const ComplexOrRealType& operator[](SizeType ind) const
    {
        // needsExp obviously, but disabled for performance here
        assert(ind < data_.size());
        return data_[ind];
    }

    const VectorType& toVector() const
    {
        // cop out for now; remove later
        needsExp("toVector");
        return data_;
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

    friend ComplexOrRealType operator*(const QuasiVector& a, const QuasiVector& b)
    {
        return a.toVector()*b.toVector();
    }

    friend std::ostream& operator<<(std::ostream& os, const QuasiVector& qv)
    {
        qv.needsExp("operator<<");
        for (SizeType i = 0; i < qv.data_.size(); ++i) {
            os<<qv.data_[i]<<" ";
        }

        return os;
    }

    friend RealType diffVectorDiff2(const QuasiVector& v1,
                                    const QuasiVector& v2,
                                    const QuasiVector& v3)
    {
        return diffVectorDiff2_(v1.toVector(), v2.toVector(), v3.toVector());
    }

private:

    void needsExp(const std::string& info) const
    {
        if (isExp_) return;
        err(info + " unimplemented or non-working unless exponential representation\n");
    }

    static RealType vectorDiff2(const VectorType& v1, const VectorType& v2)
    {
        const SizeType n = v1.size();
        assert(n == v2.size());
        RealType sum = 0;
        for (SizeType i = 0; i < n; ++i)
            sum += std::abs(v2[i] - v1[i]);

        return sum/n;
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
            if (denom == 0) denom = 1;
            const RealType re1 = PsimagLite::real(v2[i] - v1[i]);
            const RealType im1 = PsimagLite::imag(v2[i] - v1[i]);
            sum += (PsimagLite::real(v3[i])*re1 + PsimagLite::imag(v3[i])*im1)/denom;
        }

        return sum/n;
    }

    SizeType size_;
    bool isExp_;
    VectorType data_;
};
}
#endif // QUASIVECTOR_HH
