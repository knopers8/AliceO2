// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   MIDBase/GeometryTransformer.h
/// \brief  Geometry transformer for MID
/// \author Diego Stocco <Diego.Stocco at cern.ch>
/// \date   07 July 2017

#ifndef O2_MID_GEOMETRYTRANSFORMER_H
#define O2_MID_GEOMETRYTRANSFORMER_H

#include <array>
#include "MathUtils/Cartesian3D.h"

namespace o2
{
namespace mid
{
/// Class to convert the local MID RPC coordinates into global ones
/// and viceversa
class GeometryTransformer
{
 public:
  GeometryTransformer();
  virtual ~GeometryTransformer() = default;

  GeometryTransformer(const GeometryTransformer&) = delete;
  GeometryTransformer& operator=(const GeometryTransformer&) = delete;
  GeometryTransformer(GeometryTransformer&&) = delete;
  GeometryTransformer& operator=(GeometryTransformer&&) = delete;

  template <typename T>
  Point3D<T> localToGlobal(int deId, const Point3D<T>& position) const
  {
    /// Converts local coordinates into global ones
    return mTransformations[deId](position);
  }
  template <typename T>
  Point3D<T> globalToLocal(int deId, const Point3D<T>& position) const
  {
    /// Converts global coordinates into local ones
    return mTransformations[deId].ApplyInverse(position);
  }
  template <typename T>
  Point3D<T> localToGlobal(int deId, T xPos, T yPos) const
  {
    /// Converts local coordinates into global ones
    return localToGlobal(deId, Point3D<T>(xPos, yPos, 0.));
  }
  template <typename T>
  Point3D<T> globalToLocal(int deId, T xPos, T yPos, T zPos) const
  {
    /// Converts global coordinates into local ones
    return globalToLocal(deId, Point3D<T>(xPos, yPos, zPos));
  }
  template <typename T>
  Vector3D<T> localToGlobal(int deId, const Vector3D<T>& direction) const
  {
    /// Converts direction in local coordinates into global ones
    return mTransformations[deId](direction);
  }
  template <typename T>
  Vector3D<T> globalToLocal(int deId, const Vector3D<T>& direction) const
  {
    /// Converts direction in global coordinates into a local ones
    return mTransformations[deId].ApplyInverse(direction);
  }

 private:
  void init();
  std::array<ROOT::Math::Transform3D, 72> mTransformations; ///< Array of transformation matrices
};
} // namespace mid
} // namespace o2

#endif /* O2_MID_GEOMETRYTRANSFORMER_H */
