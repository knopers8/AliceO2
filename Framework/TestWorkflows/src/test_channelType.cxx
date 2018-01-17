// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.



#include "Framework/DataProcessorSpec.h"
#include "FairMQLogger.h"
#include "Framework/ParallelContext.h"
#include "Framework/ControlService.h"
#include <vector>
#include <TString.h>
#include "Framework/ChannelConfigurationPolicy.h"
#include "Framework/ChannelConfigurationPolicyHelpers.h"
#include "Framework/ChannelSpec.h"

using namespace o2::framework;
using DataHeader = o2::header::DataHeader;

void customize(std::vector<o2::framework::ChannelConfigurationPolicy> &policies){

  ChannelConfigurationPolicy pushPullPolicy;
  pushPullPolicy.match = [](std::string const &producerId, std::string const &consumerId) -> bool {
    return producerId == "producerPush" || consumerId == "consumerPull";
  };
  pushPullPolicy.modifyInput = ChannelConfigurationPolicyHelpers::pullInput;
  pushPullPolicy.modifyOutput = ChannelConfigurationPolicyHelpers::pushOutput;
}

#include "Framework/runDataProcessing.h"

void defineDataProcessing(o2::framework::WorkflowSpec &specs) {

  DataProcessorSpec producerPub{
    "producerPub",
    {},
    {
      OutputSpec{"PS", "STRING", 0, OutputSpec::Timeframe}
    },
    AlgorithmSpec{
      (AlgorithmSpec::ProcessCallback)[](ProcessingContext& ctx){

        sleep(1);

        auto& data = ctx.allocator().make<TObjString>(
          OutputSpec{"PS", "STRING", 0},
          "pub pub pub");
      }
    }
  };

  DataProcessorSpec consumerSub{
    "consumerSub",
    {
      InputSpec{"str", "PS", "STRING", 0, InputSpec::Timeframe}
    },
    {},
    AlgorithmSpec{
      (AlgorithmSpec::ProcessCallback)[](ProcessingContext& ctx){

        auto s = ctx.inputs().get<TObjString>("str");
        LOG(INFO) << "consumerSub string: " << s->GetString().Data();
      }
    }
  };

  specs.push_back(producerPub);
  specs.push_back(consumerSub);


  DataProcessorSpec producerPush{
    "producerPush",
    {},
    {
      OutputSpec{"PP", "STRING", 0, OutputSpec::Timeframe}
    },
    AlgorithmSpec{
      (AlgorithmSpec::ProcessCallback)[](ProcessingContext& ctx){

        sleep(1);

        auto& data = ctx.allocator().make<TObjString>(
          OutputSpec{"PP", "STRING", 0},
          "push push push");
      }
    }
  };

  DataProcessorSpec consumerPull{
    "consumerPull",
    {
      InputSpec{"str", "PP", "STRING", 0, InputSpec::Timeframe}
    },
    {},
    AlgorithmSpec{
      (AlgorithmSpec::ProcessCallback)[](ProcessingContext& ctx){

        auto s = ctx.inputs().get<TObjString>("str");
        LOG(INFO) << "consumerPull string: " << s->GetString().Data();
      }
    }
  };

  specs.push_back(producerPush);
  specs.push_back(consumerPull);

}
