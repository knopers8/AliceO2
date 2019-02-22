// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Dispatcher.cxx
/// \brief Implementation of Dispatcher for O2 Data Sampling
///
/// \author Piotr Konopka, piotr.jan.konopka@cern.ch

#include "Framework/Dispatcher.h"
#include "Framework/RawDeviceService.h"
#include "Framework/DataSamplingPolicy.h"
#include "Framework/DataProcessingHeader.h"
#include "Framework/DataSpecUtils.h"
#include "Framework/CallbackService.h"

#include <Monitoring/Monitoring.h>
#include <Configuration/ConfigurationInterface.h>
#include <Configuration/ConfigurationFactory.h>
#include <fairmq/FairMQDevice.h>
#include <fairmq/FairMQLogger.h>

#include <chrono>
#include <Framework/ControlService.h>

using namespace std::chrono;
using namespace o2::configuration;

namespace o2
{
namespace framework
{

Dispatcher::Dispatcher(std::string name, const std::string reconfigurationSource)
  : mName(name), mReconfigurationSource(reconfigurationSource)
{
}

Dispatcher::~Dispatcher() = default;

void Dispatcher::init(InitContext& ctx)
{
  LOG(DEBUG) << "Reading Data Sampling Policies...";

  std::unique_ptr<ConfigurationInterface> cfg = ConfigurationFactory::getConfiguration(mReconfigurationSource);
  auto policiesTree = cfg->getRecursive("dataSamplingPolicies");
  mPolicies.clear();

  for (auto&& policyConfig : policiesTree) {
    mPolicies.emplace_back(std::make_shared<DataSamplingPolicy>(policyConfig.second));
  }

  ctx.services().get<framework::CallbackService>().set(framework::CallbackService::Id::Stop, [this]() {
    std::ofstream file;
    file.open("data-sampling-benchmark", std::fstream::out | std::fstream::app);
    if (file) {
      if (elapsed_time_ms < 1 || number_of_passed_messages <= 1) {
        std::cerr << "elapsed_time_ms " << elapsed_time_ms << " or number_of_passed_messages " << number_of_passed_messages << " wrong" << std::endl;
      } else {
        LOG(INFO) << "elapsed_time_ms " << elapsed_time_ms << " or number_of_passed_messages " << number_of_passed_messages;
        file << std::setw(20) << number_of_passed_messages * 1000 / elapsed_time_ms;
      }
    } else {
      std::cerr << "could not open file for benchmark results" << std::endl;
    }
    file.close();
  });
}

void Dispatcher::run(ProcessingContext& ctx)
{
  static auto start = steady_clock::now();

  for (const auto& input : ctx.inputs()) {
    if (input.header != nullptr && input.spec != nullptr) {

      for (auto& policy : mPolicies) {
        // todo: consider getting the outputSpec in match to improve performance
        // todo: consider matching (and deciding) in completion policy to save some time
        if (policy->match(*input.spec) && policy->decide(input)) {

          if (!policy->getFairMQOutputChannel().empty()) {
            sendFairMQ(ctx.services().get<RawDeviceService>().device(), input, policy->getFairMQOutputChannelName());
          } else {
            send(ctx.outputs(), input, policy->prepareOutput(*input.spec));
          }
          number_of_passed_messages++;
        }
      }
    }
  }

  auto now = steady_clock::now();
  auto diff = duration_cast<milliseconds>(now - start).count();

  if ( diff > 300*1000) {
    elapsed_time_ms = diff;
    ctx.services().get<ControlService>().readyToQuit(true);
  }
}

void Dispatcher::send(DataAllocator& dataAllocator, const DataRef& inputData, const Output& output) const
{
  const auto* inputHeader = header::get<header::DataHeader*>(inputData.header);
  dataAllocator.snapshot(output, inputData.payload, inputHeader->payloadSize, inputHeader->payloadSerializationMethod);
}

// ideally this should be in a separate proxy device or use Lifetime::External
void Dispatcher::sendFairMQ(FairMQDevice* device, const DataRef& inputData, const std::string& fairMQChannel) const
{
  const auto* dh = header::get<header::DataHeader*>(inputData.header);
  assert(dh);
  const auto* dph = header::get<DataProcessingHeader*>(inputData.header);
  assert(dph);

  header::DataHeader dhout{ dh->dataDescription, dh->dataOrigin, dh->subSpecification, dh->payloadSize };
  dhout.payloadSerializationMethod = dh->payloadSerializationMethod;
  DataProcessingHeader dphout{ dph->startTime, dph->duration };
  o2::header::Stack headerStack{ dhout, dphout };

  auto channelAlloc = o2::pmr::getTransportAllocator(device->Transport());
  FairMQMessagePtr msgHeaderStack = o2::pmr::getMessage(std::move(headerStack), channelAlloc);

  char* payloadCopy = new char[dh->payloadSize];
  memcpy(payloadCopy, inputData.payload, dh->payloadSize);
  auto cleanupFcn = [](void* data, void*) { delete[] reinterpret_cast<char*>(data); };
  FairMQMessagePtr msgPayload(device->NewMessage(payloadCopy, dh->payloadSize, cleanupFcn, payloadCopy));

  FairMQParts message;
  message.AddPart(move(msgHeaderStack));
  message.AddPart(move(msgPayload));

  int64_t bytesSent = device->Send(message, fairMQChannel);
}

void Dispatcher::registerPath(const std::pair<InputSpec, OutputSpec>& path)
{
  //todo: take care of inputs inclusive in others, when subSpec matchers are supported
  auto cmp = [a = path.first](const InputSpec b)
  {
    return a.matcher == b.matcher && a.lifetime == b.lifetime;
  };

  if (std::find_if(inputs.begin(), inputs.end(), cmp) == inputs.end()) {
    inputs.push_back(path.first);
    LOG(DEBUG) << "Registering input " << DataSpecUtils::describe(path.first);
  } else {
    LOG(DEBUG) << "Input " << DataSpecUtils::describe(path.first)
               << " already registered";
  }

  outputs.push_back(path.second);
}

const std::string& Dispatcher::getName()
{
  return mName;
}

Inputs Dispatcher::getInputSpecs()
{
  return inputs;
}

Outputs Dispatcher::getOutputSpecs()
{
  return outputs;
}

} // namespace framework
} // namespace o2
