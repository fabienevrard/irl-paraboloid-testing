// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_CUBE_HOLE_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_CUBE_HOLE_TPP_

#include <cassert>

#include "irl/geometry/general/moment_calculation_through_simplices.h"
#include "irl/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace cube_hole_triangulation {
static constexpr UnsignedIndex_t datum_index = 12;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 14>
    face_triangle_decomposition{{{3, 2, 1},
                                 {3, 1, 0},
                                 {0, 1, 5},
                                 {0, 5, 4},
                                 {4, 5, 6},
                                 {4, 6, 7},
                                 {6, 2, 3},
                                 {6, 3, 7},
                                 {1, 2, 6},
                                 {1, 6, 5},
                                 {0, 4, 7},
                                 {0, 7, 3},
                                 {9, 8, 10},
                                 {9, 10, 11}}};
}  // namespace cube_hole_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType>
CubeHoleSpecialization<Derived, VertexType>::generateHalfEdgeVersion(
    void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void CubeHoleSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(40, 13, 11);

  for (UnsignedIndex_t v = 0; v < 13; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 40> ending_vertex_mapping{
      {2, 1, 0, 3, 1, 5,  4,  0, 5, 6,  7, 4,  2,  3, 7,  6,  2,  6, 5,  1,
       4, 7, 3, 0, 8, 10, 11, 9, 9, 12, 8, 11, 12, 9, 10, 12, 11, 8, 12, 10}};
  static constexpr std::array<UnsignedIndex_t, 40> previous_half_edge_mapping{
      {3,  0,  1,  2,  7,  4,  5,  6,  11, 8,  9,  10, 15, 12,
       13, 14, 19, 16, 17, 18, 23, 20, 21, 22, 27, 24, 25, 26,
       30, 28, 29, 33, 31, 32, 36, 34, 35, 39, 37, 38}};
  static constexpr std::array<UnsignedIndex_t, 40> next_half_edge_mapping{
      {1,  2,  3,  0,  5,  6,  7,  4,  9,  10, 11, 8,  13, 14,
       15, 12, 17, 18, 19, 16, 21, 22, 23, 20, 25, 26, 27, 24,
       29, 30, 28, 32, 33, 31, 35, 36, 34, 38, 39, 37}};
  static constexpr std::array<UnsignedIndex_t, 40> face_mapping{
      {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4,  4,  4,
       5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10}};
  static constexpr std::array<UnsignedIndex_t, 40> opposite_half_edge_mapping{
      {13, 16, 4,  23, 2,  19, 8,  20, 6,  18, 15, 21, 17, 0,
       22, 10, 1,  12, 9,  5,  7,  11, 14, 3,  28, 37, 34, 31,
       24, 33, 38, 27, 36, 29, 26, 39, 32, 25, 30, 35}};
  for (UnsignedIndex_t n = 0;
       n < static_cast<UnsignedIndex_t>(ending_vertex_mapping.size()); ++n) {
    HalfEdgeType& current_half_edge = a_half_edge_version->getHalfEdge(n);
    current_half_edge = HalfEdgeType(
        &a_half_edge_version->getVertex(ending_vertex_mapping[n]),
        &a_half_edge_version->getHalfEdge(previous_half_edge_mapping[n]),
        &a_half_edge_version->getHalfEdge(next_half_edge_mapping[n]),
        &a_half_edge_version->getFace(face_mapping[n]));
    current_half_edge.setOppositeHalfEdge(
        &a_half_edge_version->getHalfEdge(opposite_half_edge_mapping[n]));
    current_half_edge.getFace()->setStartingHalfEdge(&current_half_edge);
    current_half_edge.getVertex()->setHalfEdge(&current_half_edge);
  }
}

template <class Derived, class VertexType>
constexpr UnsignedIndex_t CubeHoleSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      cube_hole_triangulation::face_triangle_decomposition.size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
CubeHoleSpecialization<Derived, VertexType>::getSimplexIndicesFromDecomposition(
    const UnsignedIndex_t a_tet) {
  assert(a_tet < CubeHoleSpecialization::getNumberOfSimplicesInDecomposition());
  return {cube_hole_triangulation::face_triangle_decomposition[a_tet][0],
          cube_hole_triangulation::face_triangle_decomposition[a_tet][1],
          cube_hole_triangulation::face_triangle_decomposition[a_tet][2],
          cube_hole_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived>
CubeHoleSpecialization<Derived, VertexType>::getSimplexFromDecomposition(
    const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

}  // namespace IRL
#endif  // SRC_GEOMETRY_POLYHEDRONS_CUBE_HOLE_TPP_
