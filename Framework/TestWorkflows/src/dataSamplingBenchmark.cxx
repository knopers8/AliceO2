// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/ConfigParamSpec.h"
#include "Framework/DataSampling.h"
#include <vector>

using namespace o2::framework;
void customize(std::vector<CompletionPolicy>& policies)
{
  DataSampling::CustomizeInfrastructure(policies);
}
void customize(std::vector<ChannelConfigurationPolicy>& policies)
{
  DataSampling::CustomizeInfrastructure(policies);
}

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back( ConfigParamSpec{ "payload-size", VariantType::Int, 10000, { "payload size" } });
  workflowOptions.push_back( ConfigParamSpec{ "producers", VariantType::Int, 1, { "number of producers"} });
  workflowOptions.push_back( ConfigParamSpec{ "dispatchers", VariantType::Int, 1, { "number of dispatchers"} });
  workflowOptions.push_back( ConfigParamSpec{ "usleep", VariantType::Int, 1, { "usleep time of producers"} });
}

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <memory>

#include "Framework/ControlService.h"
#include "Framework/DataSampling.h"
#include "Framework/DataSamplingPolicy.h"
#include "Framework/runDataProcessing.h"

using namespace o2::framework;

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{
  WorkflowSpec specs;

  size_t playloadSize = config.options().get<int>("payload-size");
  size_t producers = config.options().get<int>("producers");
  size_t dispatchers = config.options().get<int>("dispatchers");
  size_t usleeptime = config.options().get<int>("usleep");

  for (size_t i = 0; i < producers; i++) {
    specs.push_back(DataProcessorSpec {
      "dataProducer" + std::to_string(i),
      Inputs{},
      {
        OutputSpec{ "TST", "RAWDATA", i, Lifetime::Timeframe }
      },
      AlgorithmSpec{
        (AlgorithmSpec::ProcessCallback) [usleeptime, playloadSize, i](ProcessingContext& processingContext) {
          usleep(usleeptime);
          auto data = processingContext.outputs().make<char>(Output{ "TST", "RAWDATA", i, Lifetime::Timeframe }, playloadSize);
        }
      }
    });

  }

  std::string configurationSource = std::string("json://") + getenv("BASEDIR") + "/../../O2/Framework/TestWorkflows/src/dataSamplingBenchmark.json";
  DataSampling::GenerateInfrastructure(specs, configurationSource, dispatchers);

  DataProcessorSpec podDataSink{
    "dataSink",
    Inputs{{ "tst", { DataSamplingPolicy::createPolicyDataOrigin(), DataSamplingPolicy::createPolicyDataDescription("benchmark", 0) }}},
    Outputs{},
    AlgorithmSpec{
      (AlgorithmSpec::ProcessCallback) [](ProcessingContext& processingContext){

      }
    }
  };

  specs.push_back(podDataSink);
  return specs;
}
