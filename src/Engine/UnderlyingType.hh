#ifndef UNDERLYINGTYPE_HH
#define UNDERLYINGTYPE_HH
#include "../Primitives/QuasiVector.hh"
#include <vector>

namespace Gep {

template <typename T>
struct UnderlyingType {
	using Type = T;
};

template <typename T>
struct UnderlyingType<std::vector<T>> {
	using Type = T;
};

template <typename T>
struct UnderlyingType<QuasiVector<T>> {
	using Type = T;
};

}
#endif // UNDERLYINGTYPE_HH
