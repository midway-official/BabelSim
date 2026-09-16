#pragma once

#include "babelsim/field.h"

// Read-only geometric data exposed as ordinary Fields. These functions only
// materialize values already owned by the Mesh; they do not load a case file,
// select output, or introduce a physical algorithm.
namespace babelsim::geometry {

// Cell-centred geometry.
ScalarField cellVolumes(const Mesh& mesh);
VectorField cellCentres(const Mesh& mesh);

// Face-centred geometry. Face area vectors use the Mesh orientation: outward
// on a boundary and owner-to-neighbour on an internal face.
ScalarField faceAreas(const Mesh& mesh);
VectorField faceAreaVectors(const Mesh& mesh);
VectorField faceCentres(const Mesh& mesh);
VectorField faceUnitNormals(const Mesh& mesh);

}  // namespace babelsim::geometry
