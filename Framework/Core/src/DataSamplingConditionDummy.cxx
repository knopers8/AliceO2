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
class DataSamplingConditionDummy : public DataSamplingCondition
{

  public:
  /// \brief Constructor.
  DataSamplingConditionDummy() : DataSamplingCondition() {};
  /// \brief Default destructor
  ~DataSamplingConditionDummy() = default;

  /// \brief
  void configure(const boost::property_tree::ptree& config) override
  {
  };
  /// \brief
  bool decide(const o2::framework::DataRef& dataRef) override
  {
    const auto* dpHeader = get<DataProcessingHeader*>(dataRef.header);
    assert(dpHeader);
    return dpHeader->startTime;
  }

  uint64_t rnd(uint64_t i) override {
    return 123;
  }

  private:
};

std::unique_ptr<DataSamplingCondition> DataSamplingConditionFactory::createDataSamplingConditionDummy()
{
  return std::make_unique<DataSamplingConditionDummy>();
}

} // namespace framework
} // namespace o2