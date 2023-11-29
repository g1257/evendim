#ifndef UNDERLYINGTYPE_HH
#define UNDERLYINGTYPE_HH

namespace Gep {

template <typename T>
struct UnderlyingType {
	using Type = T;
};

template <typename T>
struct UnderlyingType<QuasiVector<T>> {
	using Type = T;
};

}
#endif // UNDERLYINGTYPE_HH
