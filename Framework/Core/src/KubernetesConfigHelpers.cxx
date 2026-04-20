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
#include "KubernetesConfigHelpers.h"
#include "DeviceSpecHelpers.h"

#include "Framework/Logger.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <fmt/core.h>
#include <boost/program_options.hpp>

namespace bfs = std::filesystem;

namespace o2::framework
{

namespace
{
constexpr const char* Indent = "  ";
constexpr std::string_view DefaultImage = "gitlab-registry.cern.ch/aliceo2group/dockerfiles/alma9-flp-node/dpl:latest";
constexpr std::string_view SharedShmVolumeName = "o2-shm";
constexpr std::string_view SharedTmpVolumeName = "o2-tmp";

std::string yamlQuote(std::string_view value)
{
  std::string result;
  result.reserve(value.size() + 2);
  result.push_back('"');
  for (char c : value) {
    switch (c) {
      case '\\':
        result += "\\\\";
        break;
      case '"':
        result += "\\\"";
        break;
      case '\n':
        result += "\\n";
        break;
      case '\r':
        result += "\\r";
        break;
      case '\t':
        result += "\\t";
        break;
      default:
        result.push_back(c);
        break;
    }
  }
  result.push_back('"');
  return result;
}

std::string sanitizeK8sName(std::string_view value)
{
  std::string result;
  result.reserve(value.size());
  for (char c : value) {
    char lowered = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (std::isalnum(static_cast<unsigned char>(lowered)) != 0) {
      result.push_back(lowered);
    } else if (lowered == '-') {
      result.push_back(lowered);
    } else {
      result.push_back('-');
    }
  }
  while (!result.empty() && result.front() == '-') {
    result.erase(result.begin());
  }
  while (!result.empty() && result.back() == '-') {
    result.pop_back();
  }
  if (result.empty()) {
    result = "o2-dpl";
  }
  constexpr size_t maxLen = 63;
  if (result.size() > maxLen) {
    result.resize(maxLen);
    while (!result.empty() && result.back() == '-') {
      result.pop_back();
    }
    if (result.empty()) {
      result = "o2-dpl";
    }
  }
  return result;
}

std::string uniqueName(std::string const& base, std::unordered_map<std::string, size_t>& seen)
{
  auto it = seen.find(base);
  if (it == seen.end()) {
    seen.emplace(base, 1);
    return base;
  }
  auto index = it->second++;
  std::string suffix = "-" + std::to_string(index);
  constexpr size_t maxLen = 63;
  std::string trimmed = base;
  if (trimmed.size() + suffix.size() > maxLen) {
    trimmed.resize(maxLen - suffix.size());
    while (!trimmed.empty() && trimmed.back() == '-') {
      trimmed.pop_back();
    }
    if (trimmed.empty()) {
      trimmed = "o2-dpl";
    }
  }
  return trimmed + suffix;
}

struct RawChannel {
  std::string_view name;
  std::string_view method;
  std::string_view address;
};

std::string_view extractValueFromChannelConfig(std::string_view config, std::string_view token)
{
  size_t tokenStart = config.find(token);
  if (tokenStart == std::string_view::npos) {
    return {};
  }
  size_t valueStart = tokenStart + token.size();
  if (valueStart >= config.size()) {
    return {};
  }
  size_t valueEnd = config.find(',', valueStart);
  return valueEnd == std::string_view::npos ? config.substr(valueStart) : config.substr(valueStart, valueEnd - valueStart);
}

std::vector<RawChannel> extractRawChannels(const DeviceSpec& spec, const DeviceExecution& execution)
{
  std::vector<std::string> dplChannels;
  dplChannels.reserve(spec.inputChannels.size() + spec.outputChannels.size());
  for (const auto& channel : spec.inputChannels) {
    dplChannels.emplace_back(channel.name);
  }
  for (const auto& channel : spec.outputChannels) {
    dplChannels.emplace_back(channel.name);
  }

  std::vector<RawChannel> rawChannels;
  for (size_t i = 0; i < execution.args.size(); i++) {
    if (execution.args[i] != nullptr && strcmp(execution.args[i], "--channel-config") == 0 && i + 1 < execution.args.size()) {
      auto channelConfig = std::string_view{execution.args[i + 1]};
      auto channelName = extractValueFromChannelConfig(channelConfig, "name=");
      if (std::find(dplChannels.begin(), dplChannels.end(), channelName) == dplChannels.end()) {
        rawChannels.push_back({channelName,
                               extractValueFromChannelConfig(channelConfig, "method="),
                               extractValueFromChannelConfig(channelConfig, "address=")});
      }
    }
  }
  return rawChannels;
}

std::optional<int> parseTcpPort(std::string_view address)
{
  constexpr std::string_view tcpPrefix = "tcp://";
  if (!address.starts_with(tcpPrefix)) {
    return std::nullopt;
  }
  auto hostPort = address.substr(tcpPrefix.size());
  auto colon = hostPort.rfind(':');
  if (colon == std::string_view::npos || colon + 1 >= hostPort.size()) {
    return std::nullopt;
  }
  auto portPart = hostPort.substr(colon + 1);
  auto end = portPart.find_first_of(",;");
  if (end != std::string_view::npos) {
    portPart = portPart.substr(0, end);
  }
  if (portPart.empty()) {
    return std::nullopt;
  }
  try {
    return std::stoi(std::string(portPart));
  } catch (...) {
    return std::nullopt;
  }
}

std::set<int> collectTcpBindPorts(const DeviceSpec& spec, const DeviceExecution& execution)
{
  std::set<int> ports;
  for (const auto& channel : spec.outputChannels) {
    if (channel.method == ChannelMethod::Bind && channel.protocol == ChannelProtocol::Network && channel.port != 0) {
      ports.insert(static_cast<int>(channel.port));
    }
  }
  for (const auto& channel : spec.inputChannels) {
    if (channel.method == ChannelMethod::Bind && channel.protocol == ChannelProtocol::Network && channel.port != 0) {
      ports.insert(static_cast<int>(channel.port));
    }
  }
  auto rawChannels = extractRawChannels(spec, execution);
  for (const auto& rawChannel : rawChannels) {
    if (rawChannel.method == "bind") {
      auto port = parseTcpPort(rawChannel.address);
      if (port.has_value()) {
        ports.insert(port.value());
      }
    }
  }
  return ports;
}

std::pair<std::string, std::string> splitEnv(std::string_view env)
{
  auto pos = env.find('=');
  if (pos == std::string_view::npos) {
    return {std::string(env), ""};
  }
  return {std::string(env.substr(0, pos)), std::string(env.substr(pos + 1))};
}

std::unordered_map<std::string, std::string> getDefaultValues(const DeviceSpec& spec,
                                                               const DataProcessorInfo* processorInfo)
{
  std::unordered_map<std::string, std::string> defaults;

  // Add defaults from spec.options (device-specific options)
  for (const auto& option : spec.options) {
    if (option.defaultValue.type() != VariantType::Empty) {
      defaults[option.name] = option.defaultValue.asString();
    }
  }

  // Add defaults from workflowOptions if processorInfo is available
  if (processorInfo != nullptr) {
    for (const auto& option : processorInfo->workflowOptions) {
      if (option.defaultValue.type() != VariantType::Empty) {
        defaults[option.name] = option.defaultValue.asString();
      }
    }
  }

  // Add defaults from forwarded device options
  // We parse an empty command line with the options to get defaults populated
  auto forwardedOptions = DeviceSpecHelpers::getForwardedDeviceOptions();
  boost::program_options::variables_map vm;
  try {
    int emptyArgc = 0;
    char** emptyArgv = nullptr;
    boost::program_options::store(
      boost::program_options::parse_command_line(emptyArgc, emptyArgv, forwardedOptions),
      vm);
    boost::program_options::notify(vm);

    for (const auto& entry : vm) {
      const auto& name = entry.first;
      const auto& value = entry.second;
      if (!value.defaulted()) {
        continue; // Not a default value
      }
      try {
        if (auto strVal = boost::any_cast<std::string>(&value.value())) {
          defaults[name] = *strVal;
        } else if (auto boolVal = boost::any_cast<bool>(&value.value())) {
          defaults[name] = *boolVal ? "true" : "false";
        } else if (auto intVal = boost::any_cast<int>(&value.value())) {
          defaults[name] = std::to_string(*intVal);
        } else if (auto uintVal = boost::any_cast<unsigned int>(&value.value())) {
          defaults[name] = std::to_string(*uintVal);
        } else if (auto ushortVal = boost::any_cast<unsigned short>(&value.value())) {
          defaults[name] = std::to_string(*ushortVal);
        }
      } catch (...) {
        // Skip if we can't convert
      }
    }
  } catch (...) {
    // If parsing fails, skip forwarded options defaults
  }

  return defaults;
}
} // namespace

void KubernetesConfigHelpers::dumpPodManifest(std::ostream& out,
                                              std::string const& workflowName,
                                              std::vector<DeviceSpec> const& specs,
                                              std::vector<DeviceExecution> const& executions,
                                              std::vector<DataProcessorInfo> const& dataProcessorInfos,
                                              CommandInfo const&)
{
  assert(specs.size() == executions.size());
  std::string podName = sanitizeK8sName(workflowName.empty() ? "o2-dpl" : workflowName);

  out << "apiVersion: v1\n";
  out << "kind: Pod\n";
  out << "metadata:\n";
  out << Indent << "name: " << yamlQuote(podName) << "\n";
  out << Indent << "labels:\n";
  out << Indent << Indent << "app: " << yamlQuote(podName) << "\n";
  out << "spec:\n";
  out << Indent << "restartPolicy: Never\n";
  out << Indent << "volumes:\n";
  out << Indent << "- name: " << SharedShmVolumeName << "\n";
  out << Indent << Indent << "emptyDir:\n";
  out << Indent << Indent << Indent << "medium: Memory\n";
  out << Indent << "- name: " << SharedTmpVolumeName << "\n";
  out << Indent << Indent << "emptyDir: {}\n";
  out << Indent << "containers:\n";

  std::unordered_map<std::string, size_t> seenNames;
  for (size_t di = 0; di < specs.size(); ++di) {
    const auto& spec = specs[di];
    const auto& execution = executions[di];
    if (execution.args.empty() || execution.args[0] == nullptr) {
      continue;
    }
    std::string containerName = uniqueName(sanitizeK8sName(spec.id), seenNames);
    out << Indent << "- name: " << yamlQuote(containerName) << "\n";
    out << Indent << Indent << "image: " << yamlQuote(DefaultImage) << "\n";
    out << Indent << Indent << "command:\n";
    out << Indent << Indent << Indent << "- /bin/sh\n";
    out << Indent << Indent << "args:\n";
    out << Indent << Indent << Indent << "- -c\n";

    // Build the full command string with sourcing o2.sh first
    std::string fullCommand = "source /etc/profile.d/o2.sh && ";

    // Extract just the binary name from the path (e.g., "./stage/bin/foo" -> "foo")
    std::string executablePath = execution.args[0];
    auto lastSlash = executablePath.find_last_of("/\\");
    std::string binaryName = (lastSlash != std::string::npos)
                              ? executablePath.substr(lastSlash + 1)
                              : executablePath;
    fullCommand += binaryName;

    // Find the matching DataProcessorInfo for this device
    auto pi = std::find_if(dataProcessorInfos.begin(), dataProcessorInfos.end(),
                           [&spec](auto const& x) { return x.name == spec.name; });
    const DataProcessorInfo* processorInfo = (pi != dataProcessorInfos.end()) ? &(*pi) : nullptr;

    // Get default values for filtering
    auto defaults = getDefaultValues(spec, processorInfo);

    for (size_t ai = 1; ai < execution.args.size(); ++ai) {
      if (execution.args[ai] == nullptr) {
        break;
      }
      const char* option = execution.args[ai];
      const char* value = nullptr;
      // If the subsequent option exists and does not start with -, we assume
      // it is an argument to the previous one. Special case: if it starts with - but is
      // a negative number (e.g., "-1", "-3.14"), treat it as a value, not an option.
      if (ai + 1 < execution.args.size() && execution.args[ai + 1] != nullptr) {
        const char* nextArg = execution.args[ai + 1];
        bool isValue = nextArg[0] != '-'; // doesn't start with dash

        // Check if it's a negative number (starts with - but followed by a digit)
        if (!isValue && nextArg[0] == '-' && nextArg[1] != '\0' &&
            (std::isdigit(nextArg[1]) || nextArg[1] == '.')) {
          isValue = true; // It's a negative number, treat as value
        }

        if (isValue) {
          value = nextArg;
          ai++;
        }
      }

      // Skip options with default values
      std::string optionName = option;
      if (optionName.starts_with("--")) {
        optionName = optionName.substr(2);
      } else if (optionName.starts_with("-")) {
        optionName = optionName.substr(1);
      }

      // Check if this option has a default value and if the current value matches it
      auto defaultIt = defaults.find(optionName);
      if (defaultIt != defaults.end() && value != nullptr) {
        if (defaultIt->second == value) {
          // Skip this argument as it matches the default
          continue;
        }
      }

      fullCommand += " ";
      fullCommand += option;
      if (value) {
        fullCommand += " ";
        fullCommand += fmt::format("'{}'", value);
      }
    }
    out << Indent << Indent << Indent << "- " << yamlQuote(fullCommand) << "\n";
    if (!execution.environ.empty()) {
      out << Indent << Indent << "env:\n";
      for (const auto& env : execution.environ) {
        auto entry = splitEnv(env);
        out << Indent << Indent << Indent << "- name: " << yamlQuote(entry.first) << "\n";
        out << Indent << Indent << Indent << Indent << "value: " << yamlQuote(entry.second) << "\n";
      }
    }
    out << Indent << Indent << "volumeMounts:\n";
    out << Indent << Indent << Indent << "- name: " << SharedShmVolumeName << "\n";
    out << Indent << Indent << Indent << Indent << "mountPath: /dev/shm\n";
    out << Indent << Indent << Indent << "- name: " << SharedTmpVolumeName << "\n";
    out << Indent << Indent << Indent << Indent << "mountPath: /tmp\n";

    auto ports = collectTcpBindPorts(spec, execution);
    if (!ports.empty()) {
      out << Indent << Indent << "ports:\n";
      for (int port : ports) {
        out << Indent << Indent << Indent << "- containerPort: " << port << "\n";
      }
    }
  }
}

void KubernetesConfigHelpers::dumpDeviceSpec2Kubernetes(std::string const& workflowName,
                                                        std::vector<DeviceSpec> const& specs,
                                                        std::vector<DeviceExecution> const& executions,
                                                        std::vector<DataProcessorInfo> const& dataProcessorInfos,
                                                        CommandInfo const& commandInfo)
{
  const char* manifestsDirectory = "kubernetes";

  LOG(info) << "Dumping the workflow configuration for Kubernetes.";
  LOG(info) << "Creating directory '" << manifestsDirectory << "'.";
  bfs::create_directory(manifestsDirectory);
  LOG(info) << "... created.";

  std::string manifestName = workflowName.empty() ? "o2-dpl" : workflowName;
  std::string manifestPath = std::string(manifestsDirectory) + bfs::path::preferred_separator + manifestName + ".yaml";
  LOG(info) << "Creating a Kubernetes manifest dump '" << manifestName << "'.";
  std::ofstream manifestOut(manifestPath);
  dumpPodManifest(manifestOut, manifestName, specs, executions, dataProcessorInfos, commandInfo);
  manifestOut.close();
}

} // namespace o2::framework
