// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_IRL_POLY_H_
#define SRC_GEOMETRY_POLYHEDRONS_IRL_POLY_H_

#include "irl/geometry/general/geometry_type_traits.h"
#include "irl/geometry/general/stored_vertex_access.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/tet.h"
#include "irl/parameters/defined_types.h"

namespace IRL{

template <class Derived, class VertexType>
class IRLPolySpecialization: public BasePolyhedron<Derived, VertexType, ProxyTet<Derived>> {
public: 

HalfEdgePolyhedron<VertexType> generateHalfEdgeVersion(void) const; 

template<class HalfEdgePolyhedronType>
void setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const;

static constexpr UnsignedIndex_t getNumberOfSimplicesInDecomposition(void);

static constexpr std::array<UnsignedIndex_t, 4> getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet);

ProxyTet<Derived> getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const;

}; 

template <class VertexType>
class StoredIRLPoly: public StoredVertexAccess<StoredIRLPoly<VertexType>,VertexType,32>, public IRLPolySpecialization<StoredIRLPoly<VertexType>, VertexType>{
friend StoredVertexAccess<StoredIRLPoly<VertexType>,VertexType,32>;

public: 

using StoredVertexAccess<StoredIRLPoly<VertexType>,VertexType,32>::StoredVertexAccess;

StoredIRLPoly(void) = default;
}; 


// Predefined types 
using IRLPoly= StoredIRLPoly<Pt>; 


template<class VertexType>
struct is_polyhedron<StoredIRLPoly<VertexType>> : std::true_type {};
} // namespace IRL 

#include "irl_poly.tpp"
#endif //SRC_GEOMETRY_POLYHEDRONS_IRL_POLY_H_
