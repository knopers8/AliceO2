// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataSamplingConditionRandomFast.cxx
/// \brief Implementation of random DataSamplingCondition
///
/// \author Piotr Konopka, piotr.jan.konopka@cern.ch

#include "DataSampling/DataSamplingCondition.h"
#include "DataSampling/DataSamplingConditionFactory.h"
#include "Framework/DataProcessingHeader.h"

#include <boost/property_tree/ptree.hpp>

using namespace o2::framework;

namespace o2::utilities
{

// todo: consider using run number as a seed
using namespace o2::header;

/// \brief A DataSamplingCondition which makes decisions randomly, but with determinism.
class DataSamplingConditionRandomFast : public DataSamplingCondition
{

 public:
  /// \brief Constructor.
  DataSamplingConditionRandomFast() : DataSamplingCondition(),
                                      mSeed(static_cast<uint64_t>(time(nullptr))){};
  /// \brief Default destructor
  ~DataSamplingConditionRandomFast() override = default;

  /// \brief Reads 'fraction' parameter (type double, between 0 and 1) and seed (int).
  void configure(const boost::property_tree::ptree& config) override
  {
    mFraction = config.get<double>("fraction");
    mSeed = config.get<uint64_t>("seed");
  };
  /// \brief Makes pseudo-random, deterministic decision based on TimesliceID.
  /// The reason behind using TimesliceID is to ensure, that data of the same events is sampled even on different FLPs.
  bool decide(const o2::framework::DataRef& dataRef) override
  {
    const auto* dpHeader = get<DataProcessingHeader*>(dataRef.header);
    assert(dpHeader);

    // https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
    // http://xorshift.di.unimi.it/splitmix64.c
    uint64_t x = dpHeader->startTime * mSeed;
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);

    return x < mFraction * static_cast<double>(std::numeric_limits<uint64_t>::max());
  }

  uint64_t rnd(uint64_t i)
  {
    uint64_t x = i * mSeed;
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
  }

 private:
  uint64_t mSeed;
  double mFraction;
};

std::unique_ptr<DataSamplingCondition> DataSamplingConditionFactory::createDataSamplingConditionRandomFast()
{
  return std::make_unique<DataSamplingConditionRandomFast>();
}

} // namespace o2::utilities
