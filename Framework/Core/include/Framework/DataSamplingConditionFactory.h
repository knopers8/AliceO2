// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_DATASAMPLINGCONDITIONFACTORY_H
#define ALICEO2_DATASAMPLINGCONDITIONFACTORY_H

/// \file DataSamplingConditionFactory.h
/// \brief A definition of DataSamplingConditionFactory
///
/// \author Piotr Konopka, piotr.jan.konopka@cern.ch

#include "DataSamplingCondition.h"

namespace o2
{
namespace framework
{

/// A factory of DataSamplingConditions children.
class DataSamplingConditionFactory
{
 public:
  /// \brief Creates instance of DataSamplingCondition child, given the name
  static std::unique_ptr<DataSamplingCondition> create(std::string name);

  // a list of getters of specific DataSamplingCondition's, they should be implemented
  // inside particular DataSamplingCondition*.cxx files.
  /// \brief Getter for DataSamplingConditionRandom
  static std::unique_ptr<DataSamplingCondition> createDataSamplingConditionRandom();
  /// \brief Getter for DataSamplingConditionBBS
  static std::unique_ptr<DataSamplingCondition> createDataSamplingConditionBBS();
  /// \brief Getter for DataSamplingConditionPCG
  static std::unique_ptr<DataSamplingCondition> createDataSamplingConditionPCG();
  /// \brief Getter for DataSamplingConditionPayloadSize
  static std::unique_ptr<DataSamplingCondition> createDataSamplingConditionPayloadSize();
  /// \brief Getter for DataSamplingConditionNConsecutive
  static std::unique_ptr<DataSamplingCondition> createDataSamplingConditionNConsecutive();
  /// \brief Getter for DataSamplingConditionHash
  static std::unique_ptr<DataSamplingCondition> createDataSamplingConditionHash();
  /// \brief Getter for DataSamplingConditionHashCombine
  static std::unique_ptr<DataSamplingCondition> createDataSamplingConditionHashCombine();
  /// \brief Getter for DataSamplingConditionTRandom1
  static std::unique_ptr<DataSamplingCondition> createDataSamplingConditionTRandom1();
  /// \brief Getter for DataSamplingConditionTRandom2
  static std::unique_ptr<DataSamplingCondition> createDataSamplingConditionTRandom2();
  /// \brief Getter for DataSamplingConditionTRandom3
  static std::unique_ptr<DataSamplingCondition> createDataSamplingConditionTRandom3();
  /// \brief Getter for DataSamplingConditionDummy
  static std::unique_ptr<DataSamplingCondition> createDataSamplingConditionDummy();


};

} // namespace framework
} // namespace o2

#endif //ALICEO2_DATASAMPLINGCONDITIONFACTORY_H
