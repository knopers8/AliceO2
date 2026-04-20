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
#ifndef O2_FRAMEWORK_KUBERNETESCONFIGHELPERS_H_
#define O2_FRAMEWORK_KUBERNETESCONFIGHELPERS_H_

#include "Framework/CommandInfo.h"
#include "Framework/DataProcessorInfo.h"
#include "Framework/DeviceExecution.h"
#include "Framework/DeviceSpec.h"

#include <iosfwd>
#include <string>
#include <vector>

namespace o2::framework
{

/// Helper to dump a Kubernetes Pod manifest for a DPL workflow.
struct KubernetesConfigHelpers {
  static void dumpDeviceSpec2Kubernetes(std::string const& workflowName,
                                        std::vector<DeviceSpec> const& specs,
                                        std::vector<DeviceExecution> const& executions,
                                        std::vector<DataProcessorInfo> const& dataProcessorInfos,
                                        CommandInfo const& commandInfo);

  static void dumpPodManifest(std::ostream& out,
                              std::string const& workflowName,
                              std::vector<DeviceSpec> const& specs,
                              std::vector<DeviceExecution> const& executions,
                              std::vector<DataProcessorInfo> const& dataProcessorInfos,
                              CommandInfo const& commandInfo);
};

} // namespace o2::framework
#endif // O2_FRAMEWORK_KUBERNETESCONFIGHELPERS_H_
