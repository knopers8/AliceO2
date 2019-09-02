// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/runDataProcessing.h"

using namespace o2::framework;

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return
  {
    {
      "A1",
      Inputs{},
      Outputs{{{"A-out"}, "AORG", "ADESCR"}},
      (AlgorithmSpec::ProcessCallback)[](ProcessingContext & ctx){
        sleep(1);
        auto aData = ctx.outputs().make<int>(OutputRef{"A-out"}, 1);
      }
    },
    {
    "A2",
      Inputs{{"A-in", "AORG", "ADESCR"}},
      Outputs{},
      (AlgorithmSpec::ProcessCallback)[](ProcessingContext & ctx) {}
    }
  };
}