// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS__H_
#define SRC_GEOMETRY_POLYHEDRONS__H_

#include "irl/geometry/general/geometry_type_traits.h"
#include "irl/geometry/general/stored_vertex_access.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/tet.h"
#include "irl/parameters/defined_types.h"

namespace IRL{

template <class Derived, class VertexType>
class Specialization: public BasePolyhedron<Derived, VertexType, ProxyTet<Derived>> {
public: 

HalfEdgePolyhedron<VertexType> generateHalfEdgeVersion(void) const; 

template<class HalfEdgePolyhedronType>
void setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const;

static constexpr UnsignedIndex_t getNumberOfSimplicesInDecomposition(void);

static constexpr std::array<UnsignedIndex_t, 4> getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet);

ProxyTet<Derived> getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const;

}; 

template <class VertexType>
class Stored: public StoredVertexAccess<Stored<VertexType>,VertexType,0>, public Specialization<Stored<VertexType>, VertexType>{
friend StoredVertexAccess<Stored<VertexType>,VertexType,0>;

public: 

using StoredVertexAccess<Stored<VertexType>,VertexType,0>::StoredVertexAccess;

Stored(void) = default;
}; 


// Predefined types 
using = Stored<Pt>; 


template<class VertexType>
struct is_polyhedron<Stored<VertexType>> : std::true_type {};
} // namespace IRL 

#include ".tpp"
#endif //SRC_GEOMETRY_POLYHEDRONS__H_
