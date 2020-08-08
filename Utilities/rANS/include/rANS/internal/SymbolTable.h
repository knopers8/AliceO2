// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   SymbolTable.h
/// @author Michael Lettrich
/// @since  2019-06-21
/// @brief  Container for information needed to encode/decode a symbol of the alphabet

#ifndef RANS_INTERNAL_SYMBOLTABLE_H
#define RANS_INTERNAL_SYMBOLTABLE_H

#include <vector>
#include <cstdint>
#include <cmath>
#include <fairlogger/Logger.h>

#include "SymbolStatistics.h"

namespace o2
{
namespace rans
{
namespace internal
{

template <typename T>
class SymbolTable
{

 public:
  SymbolTable(const SymbolStatistics& symbolStats, size_t probabiltyBits);

  const T& operator[](int64_t index) const;

  const T& getLiteralSymbol() const;

  bool isLiteralSymbol(int64_t index) const;

  size_t getAlphabetRangeBits() const;
  int getMinSymbol() const { return mMin; }
  int getMaxSymbol() const { return mMax; }

 private:
  int mMin;
  int mMax;
  std::vector<T*> mIndex;
  std::vector<T> mSymbols;
};

template <typename T>
SymbolTable<T>::SymbolTable(const SymbolStatistics& symbolStats, size_t probabiltyBits) : mMin(symbolStats.getMinSymbol()), mMax(symbolStats.getMaxSymbol()), mIndex(), mSymbols()
{
  LOG(trace) << "start building symbol table";

  mIndex.reserve(symbolStats.size());
  mSymbols.reserve(symbolStats.size());

  for (const auto& entry : symbolStats) {
    const auto [symFrequency, symCumulated] = entry;
    mSymbols.emplace_back(symCumulated, symFrequency, probabiltyBits);
    mIndex.emplace_back(&mSymbols.back());
  }

// advanced diagnostics for debug builds
#if !defined(NDEBUG)
  LOG(debug2) << "SymbolTableProperties: {"
              << "entries:" << mSymbols.size() << ", "
              << "sizeB: " << mSymbols.size() * sizeof(T) + mIndex.size() * sizeof(T*) << "}";
#endif

  LOG(trace) << "done building symbol table";
}

template <typename T>
inline const T& SymbolTable<T>::operator[](int64_t index) const
{
  const auto idx = index - mMin;
  assert(idx >= 0);
  assert(idx < mSymbols.size());
  return *(mIndex[idx]);
}

template <typename T>
inline bool SymbolTable<T>::isLiteralSymbol(int64_t index) const
{
  const auto idx = index - mMin;
  assert(idx >= 0);
  assert(idx < mSymbols.size());
  return idx == mSymbols.size() - 1;
}

template <typename T>
inline const T& SymbolTable<T>::getLiteralSymbol() const
{
  return mSymbols.back();
}

template <typename T>
inline size_t SymbolTable<T>::getAlphabetRangeBits() const
{
  if (mMax - mMin > 0) {
    return std::ceil(std::log2(mMax - mMin));
  } else if (mMax - mMin == 0) {
    return 1;
  } else {
    return 0;
  }
}

} // namespace internal
} // namespace rans
} // namespace o2

#endif /* RANS_INTERNAL_SYMBOLTABLE_H */
