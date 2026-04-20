// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Mocking.h"
#include <catch_amalgamated.hpp>
#include "../src/KubernetesConfigHelpers.h"
#include "../src/DeviceSpecHelpers.h"
#include "../src/SimpleResourceManager.h"
#include "../src/ComputingResourceHelpers.h"
#include "Framework/DeviceControl.h"
#include "Framework/DeviceSpec.h"
#include "Framework/WorkflowSpec.h"
#include "Framework/DriverConfig.h"
#include "Framework/ConfigContext.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ConfigParamStore.h"

#include <sstream>

using namespace o2::framework;

namespace
{
WorkflowSpec defineDataProcessing()
{
  return {{.name = "A",
           .outputs = Outputs{OutputSpec{"TST", "A1"}},
           .options = {ConfigParamSpec{"channel-config", VariantType::String,
                                       "name=raw_tcp,type=push,method=bind,address=tcp://*:42000,transport=zeromq",
                                       {"Raw TCP channel config"}}}},
          {.name = "B",
           .inputs = Inputs{InputSpec{"x", "TST", "A1"}}}};
}
} // namespace

TEST_CASE("TestKubernetesManifestDump")
{
  auto workflow = defineDataProcessing();
  std::ostringstream ss{""};
  auto configContext = makeEmptyConfigContext();
  auto channelPolicies = makeTrivialChannelPolicies(*configContext);
  std::vector<DeviceSpec> devices;
  std::vector<ComputingResource> resources{ComputingResourceHelpers::getLocalhostResource()};
  SimpleResourceManager rm(resources);
  auto completionPolicies = CompletionPolicy::createDefaultPolicies();
  auto callbacksPolicies = CallbacksPolicy::createDefaultPolicies();
  DeviceSpecHelpers::dataProcessorSpecs2DeviceSpecs(workflow, channelPolicies, completionPolicies, callbacksPolicies, devices, rm, "workflow-id", *configContext, true);
  std::vector<DeviceControl> controls;
  std::vector<DeviceExecution> executions;
  controls.resize(devices.size());
  executions.resize(devices.size());

  std::vector<ConfigParamSpec> workflowOptions = {
    ConfigParamSpec{"jobs", VariantType::Int, 1, {"number of producer jobs"}}};

  std::vector<DataProcessorInfo> dataProcessorInfos = {
    {
      {.name = "A", .executable = "foo", .workflowOptions = workflowOptions},
      {.name = "B", .executable = "foo", .workflowOptions = workflowOptions},
    }};

  DriverConfig driverConfig{
    .batch = false,
  };
  DeviceSpecHelpers::prepareArguments(false, false, false, 8080,
                                      driverConfig,
                                      dataProcessorInfos,
                                      devices, executions, controls, {},
                                      "workflow-id");

  CommandInfo commandInfo{"foo"};
  KubernetesConfigHelpers::dumpPodManifest(ss, "testwf", devices, executions, dataProcessorInfos, commandInfo);
  auto manifest = ss.str();

  REQUIRE(manifest.find("apiVersion: v1") != std::string::npos);
  REQUIRE(manifest.find("kind: Pod") != std::string::npos);
  REQUIRE(manifest.find("name: \"testwf\"") != std::string::npos);
  REQUIRE(manifest.find("mountPath: /dev/shm") != std::string::npos);
  REQUIRE(manifest.find("mountPath: /tmp") != std::string::npos);
  REQUIRE(manifest.find("image: \"aliceo2/o2:latest\"") != std::string::npos);
  REQUIRE(manifest.find("containers:\n  - name: \"a\"") != std::string::npos);
  REQUIRE(manifest.find("containerPort: 42000") != std::string::npos);
}
