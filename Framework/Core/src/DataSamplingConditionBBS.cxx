// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataSamplingConditionBBS.cxx
/// \brief Implementation of random DataSamplingCondition
///
/// \author Piotr Konopka, piotr.jan.konopka@cern.ch

#include "Framework/DataSamplingCondition.h"
#include "Framework/DataSamplingConditionFactory.h"
#include "Framework/DataProcessingHeader.h"
//#include "gmpbbs.h"

namespace o2
{
namespace framework
{

using namespace o2::header;

/// \brief A DataSamplingCondition which makes decisions randomly, but with determinism.
class DataSamplingConditionBBS : public DataSamplingCondition
{

 public:
  /// \brief Constructor.
  DataSamplingConditionBBS() : DataSamplingCondition() {
//    bbs = rndbbs_new();
  };
  /// \brief Default destructor
  ~DataSamplingConditionBBS() = default;

  /// \brief Reads 'fraction' parameter (type double, between 0 and 1) and seed (int).
  void configure(const boost::property_tree::ptree& cfg) override
  {

  };
  /// \brief Makes pseudo-random, deterministic decision based on TimesliceID.
  /// The reason behind using TimesliceID is to ensure, that data of the same events is sampled even on different FLPs.
  bool decide(const o2::framework::DataRef& dataRef) override
  {
    return true;
  }

 private:
  // rndbbs_t* bbs;
};

std::unique_ptr<DataSamplingCondition> DataSamplingConditionFactory::createDataSamplingConditionBBS()
{
  return std::make_unique<DataSamplingConditionBBS>();
}

} // namespace framework
} // namespace o2
