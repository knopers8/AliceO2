// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataSamplingConditionTRandom3.cxx
/// \brief Implementation of random DataSamplingCondition
///
/// \author Piotr Konopka, piotr.jan.konopka@cern.ch

#include "Framework/DataSamplingCondition.h"
#include "Framework/DataSamplingConditionFactory.h"
#include "Framework/DataProcessingHeader.h"

#include <TRandom3.h>

namespace o2
{
namespace framework
{

// todo: choose the best PRNG (TRandom3 is fast, but i am not sure about its statistical soundness and behaviour with
// very small percents) or use completely different decision mechanism (map or formula f(timeslice) -> {0,1})
// todo: consider using run number as a seed

using namespace o2::header;

/// \brief A DataSamplingCondition which makes decisions randomly, but with determinism.
class DataSamplingConditionTRandom3 : public DataSamplingCondition
{

  public:
  /// \brief Constructor.
  DataSamplingConditionTRandom3() : DataSamplingCondition(), mSeed(static_cast<uint64_t>(time(nullptr))), mFraction(0.0), mGenerator(0){};
  /// \brief Default destructor
  ~DataSamplingConditionTRandom3() = default;

  /// \brief Reads 'fraction' parameter (type double, between 0 and 1) and seed (int).
  void configure(const boost::property_tree::ptree& cfg) override
  {
    mFraction = cfg.get<double>("fraction");
    mSeed = cfg.get<int>("seed");
  };
  /// \brief Makes pseudo-random, deterministic decision based on TimesliceID.
  /// The reason behind using TimesliceID is to ensure, that data of the same events is sampled even on different FLPs.
  bool decide(const o2::framework::DataRef& dataRef) override
  {
    const auto* dpHeader = get<DataProcessingHeader*>(dataRef.header);
    assert(dpHeader);

    mGenerator.SetSeed(dpHeader->startTime * mSeed);
    return static_cast<bool>(mGenerator.Binomial(1, mFraction));
  }

  uint64_t rnd(uint64_t i) override {
    mGenerator.SetSeed(i * mSeed);
    return mGenerator.Integer(std::numeric_limits<uint32_t>::max());
  }
  private:
  uint64_t mSeed;
  double mFraction;
  TRandom3 mGenerator;
};

std::unique_ptr<DataSamplingCondition> DataSamplingConditionFactory::createDataSamplingConditionTRandom3()
{
  return std::make_unique<DataSamplingConditionTRandom3>();
}

} // namespace framework
} // namespace o2
