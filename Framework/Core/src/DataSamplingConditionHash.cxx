// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataSamplingConditionHash.cxx
/// \brief Implementation of DataSamplingConditionHash
///
/// \author Piotr Konopka, piotr.jan.konopka@cern.ch

#include <functional>

#include "Framework/DataSamplingCondition.h"
#include "Framework/DataSamplingConditionFactory.h"
#include "Framework/DataProcessingHeader.h"

namespace o2
{
namespace framework
{

using namespace o2::header;

/// \brief A DataSamplingCondition which
class DataSamplingConditionHash : public DataSamplingCondition
{

  public:
  /// \brief Constructor.
  DataSamplingConditionHash() : DataSamplingCondition(), mSeed(static_cast<uint64_t>(time(nullptr))){};
  /// \brief Default destructor
  ~DataSamplingConditionHash() = default;

  /// \brief
  void configure(const boost::property_tree::ptree& config) override
  {
    mFraction = config.get<double>("fraction");
    mSeed = config.get<uint64_t>("seed");
  };
  /// \brief
  bool decide(const o2::framework::DataRef& dataRef) override
  {
    const auto* dpHeader = get<DataProcessingHeader*>(dataRef.header);
    assert(dpHeader);

    // https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
    // http://xorshift.di.unimi.it/splitmix64.c
//    size_t h = std::hash<uint64_t>()(static_cast<uint64_t>(dpHeader->startTime));

    uint64_t x = dpHeader->startTime * mSeed;
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);


    return x < mFraction * std::numeric_limits<uint64_t>::max();
  }

  uint64_t rnd(uint64_t i) override {
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

std::unique_ptr<DataSamplingCondition> DataSamplingConditionFactory::createDataSamplingConditionHash()
{
  return std::make_unique<DataSamplingConditionHash>();
}

} // namespace framework
} // namespace o2