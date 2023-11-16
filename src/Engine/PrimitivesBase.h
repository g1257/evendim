#ifndef PRIMITIVESBASE_H
#define PRIMITIVESBASE_H

namespace Gep {

/* PSIDOC PrimitivesBase
This is the interface for Primitives, that is, these
are the functions that a programmer needs to write
to implement a new Primitives class to use with EVENDIM.

PSIDOCCOPY PrimitivesBase::nodes
PSIDOCCOPY PrimitivesBase::dcValues
PSIDOCCOPY PrimitivesBase::dcArray
*/
template <typename ValueType_>
class PrimitivesBase {
public:

	/* PSIDOC PrimitivesBase::nodes
PSIDOCCOPY $FirstProtoBelow
Returns a vector of nodes containing all different nodes.
	 */
	virtual const VectorNodeType& nodes() const = 0;

	/* PSIDOC PrimitivesBase::dcValues
PSIDOCCOPY $FirstProtoBelow
Returns a vector of defined constant values
	 */
	virtual const VectorValueType& dcValues() const = 0;

	/* PSIDOC PrimitivesBase::dcArray
PSIDOCCOPY $FirstProtoBelow
Returns a vector of defined constant names
	 */
	virtual const VectorStringType& dcArray() const = 0;
};
}
#endif // PRIMITIVESBASE_H
