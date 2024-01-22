// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file    vectorOfTObjectsExample.cxx
/// \author  Piotr Konopka
///
/// \brief producer-subscriber workflow to see the (de)serialization of std::vector<TObject*>

#include "Framework/RootSerializationSupport.h"

using namespace o2::framework;

#include "Framework/runDataProcessing.h"
#include "Framework/Logger.h"
#include "Framework/ControlService.h"

#include <TH1F.h>
//#include <Mergers/VectorTObject.h>

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  WorkflowSpec specs;

  DataProcessorSpec producer {
    "producer",
    Inputs{},
    Outputs{{{"out"}, "TST", "VECTOR", 0, Lifetime::Sporadic}},
    AlgorithmSpec{[](ProcessingContext& ctx) {
      std::vector<TObject*> vec;
      vec.push_back(new TH1F("histo", "histo", 10, 0, 10));
      static_cast<TH1F*>(vec[0])->Fill(4);

      ctx.outputs().snapshot({"out"}, ROOTSerialized<decltype(vec)>(vec));
      delete vec[0];
      ctx.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    }}
  };
  specs.push_back(producer);

  DataProcessorSpec receiver{
    "receiver",
    Inputs{{ "vector", "TST", "VECTOR", 0, Lifetime::Sporadic }},
    Outputs{},
    AlgorithmSpec{
      [](ProcessingContext& ctx) mutable {
        LOG(info) << "receiver received data";
        auto dataRef = ctx.inputs().get("vector");
        auto vec = DataRefUtils::as<ROOTSerialized<std::vector<TObject*>>>(dataRef);
        LOG(info) << "vector size: " << vec->size();
        if (vec->size() == 1) {
          static_cast<TH1F*>(vec->at(0))->Print();
          delete vec->at(0);
          ctx.services().get<ControlService>().readyToQuit(QuitRequest::All);
        }
      }
    }
  };
  specs.push_back(receiver);

  return specs;
}
