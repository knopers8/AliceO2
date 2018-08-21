// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataSamplingConditionPCG.cxx
/// \brief Implementation of random DataSamplingCondition
///
/// \author Piotr Konopka, piotr.jan.konopka@cern.ch

#include "Framework/DataSamplingCondition.h"
#include "Framework/DataSamplingConditionFactory.h"
#include "Framework/DataProcessingHeader.h"

#include "pcg_random.h"

namespace o2
{
namespace framework
{

using namespace o2::header;

/// \brief A DataSamplingCondition which makes decisions randomly, but with determinism.
class DataSamplingConditionPCG : public DataSamplingCondition
{

 public:
  /// \brief Constructor.
  DataSamplingConditionPCG() : DataSamplingCondition(), mSeed(0), mFraction(0.0), mGenerator(), mCurrentTimesliceID (){};
  /// \brief Default destructor
  ~DataSamplingConditionPCG() = default;

  /// \brief Reads 'fraction' parameter (type double, between 0 and 1) and seed (int).
  void configure(const boost::property_tree::ptree& cfg) override
  {
    mFraction = cfg.get<double>("fraction");
    mSeed = cfg.get<int>("seed");
    mGenerator.seed(mSeed);
  };
  /// \brief Makes pseudo-random, deterministic decision based on TimesliceID.
  /// The reason behind using TimesliceID is to ensure, that data of the same events is sampled even on different FLPs.
  bool decide(const o2::framework::DataRef& dataRef) override
  {
    const auto* dpHeader = get<DataProcessingHeader*>(dataRef.header);
    assert(dpHeader);

    int64_t diff = dpHeader->startTime - mCurrentTimesliceID;
    if (diff > 0) {
      mGenerator.advance(static_cast<uint64_t>(diff));
    } else if (diff < 0) {
      mGenerator.backstep(static_cast<uint64_t>(-diff));
    } else if (diff == -1){
      //return previous result?
    }
    mCurrentTimesliceID = dpHeader->startTime + 1;
//     make sure it works EACH time for 0% and 100% (at least for 100%)

//    auto rnd = mGenerator();
//    LOG(INFO) << "---------";
//    LOG(INFO) << rnd;
    // todo: division by integer might be faster
    return static_cast<uint32_t>(mGenerator()) < mFraction*std::numeric_limits<uint32_t>::max();
//    return static_cast<bool>(mGenerator.Binomial(1, mFraction));
  }

  uint64_t rnd(uint64_t i) override {
    return mGenerator();
  }

  private:
  int mSeed;
  double mFraction;
  pcg32_fast mGenerator;
  DataProcessingHeader::StartTime mCurrentTimesliceID;
};

std::unique_ptr<DataSamplingCondition> DataSamplingConditionFactory::createDataSamplingConditionPCG()
{
  return std::make_unique<DataSamplingConditionPCG>();
}

} // namespace framework
} // namespace o2
