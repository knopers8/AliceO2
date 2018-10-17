// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// The workflow consists of two data producers, accumulator and printer. Accumulator wants to receive data from two
/// producers at the same time and send some output data each x seconds. Printer shows the result sent by accumulator.


#include "Framework/CompletionPolicy.h"
#include "Framework/CompletionPolicyHelpers.h"
#include "Framework/DeviceSpec.h"
#include "Headers/DataHeader.h"

using namespace o2::framework;
using namespace o2::header;

void customize(std::vector<CompletionPolicy>& policies)
{
  CompletionPolicy timerAndData{
    "timerAndData",
    [](DeviceSpec const& dev) -> bool { return dev.name == "accumulator"; },
    [](gsl::span<PartRef const> const& inputs) -> CompletionPolicy::CompletionOp {

      bool timerPresent = false;
      size_t dataInputsPresent = 0;
      for (auto& input : inputs) {
        if (input.header == nullptr || input.payload == nullptr) {
          continue;
        }

        const DataHeader* inputHeader = get<DataHeader*>(input.header->GetData());
        assert(inputHeader);
        if (!strncmp(inputHeader->dataDescription.str, "TIMER", 5)) {
          timerPresent = true;
        } else {
          dataInputsPresent++;
        }
      }

      if (dataInputsPresent == inputs.size() - 1) {
        // both data inputs present and maybe timer => consume data, execute periodic procedure if timer is present
        return CompletionPolicy::CompletionOp::Consume;
      } else if (timerPresent) {
//         only timer and no data => execute periodic procedure, wait for the rest of data
        if (dataInputsPresent) {
          LOG(WARN) << "DATA INPUTS PRESENT IN THIS TIMESLICE, LOSING DATA PROBABLY";
        }
        return CompletionPolicy::CompletionOp::Process;
      } else {
        // partial data => wait for the rest
        return CompletionPolicy::CompletionOp::Wait;
      }
    }
  };

  policies.push_back(timerAndData);
}

#include "Framework/runDataProcessing.h"
#include <chrono>
#include <Framework/DataProcessingHeader.h>
#include <Framework/CallbackService.h>
#include <Framework/TimesliceIndex.h>

using namespace std::chrono;

WorkflowSpec defineDataProcessing(ConfigContext const&)
{

  DataProcessorSpec producer1{
    "producer1",
    Inputs{},
    Outputs{
      {{ "output1" }, "TST", "DATA1" }
    },
    AlgorithmSpec{
      (AlgorithmSpec::InitCallback) [](InitContext& ic) {
        srand(getpid());
        auto num = std::make_shared<int>(0);
        return (AlgorithmSpec::ProcessCallback) [num](ProcessingContext& ctx) {
          sleep(rand() % 2);
          auto data1 = ctx.outputs().make<int>(OutputRef{ "output1" }, 1);
          data1[0] = (*num)++;
        };
      }
    }
  };

  DataProcessorSpec producer2{
    "producer2",
    Inputs{},
    Outputs{
      {{ "output2" }, "TST", "DATA2" }
    },
    AlgorithmSpec{
      (AlgorithmSpec::InitCallback) [](InitContext& ic) {
        srand(getpid());
        auto num = std::make_shared<int>(0);
        return (AlgorithmSpec::ProcessCallback) [num](ProcessingContext& ctx) {
          sleep(rand() % 2);
          auto data2 = ctx.outputs().make<int>(OutputRef{ "output2" }, 1);
          data2[0] = (*num)++;
        };
      }
    }
  };

  DataProcessorSpec accumulator{
    "accumulator",
    Inputs{
      { "data1", "TST", "DATA1" },
      { "data2", "TST", "DATA2" },
      { "timer", "TST", "TIMER", 0, Lifetime::Timer }
    },
    Outputs{
      {{ "sum" }, "TST", "SUM" }
    },
    AlgorithmSpec{
      (AlgorithmSpec::InitCallback) [](InitContext& ictx) {

        std::shared_ptr<int> sum = std::make_shared<int>(0);

        // note that some data will be lost anyway, when timeslice	is booked on existing partial data. we need some
        // mechanism for preventing that - either separate 'cache' for timers or timesliceIndex API call for overbooking
        // only the obsolete timeslices
        auto bookTimesliceCallback = [&timesliceIndex = ictx.services().get<TimesliceIndex>(),
          time = std::make_shared<steady_clock::time_point>(steady_clock::now())]() {
          auto timeNow = steady_clock::now();
          static TimesliceId i{0};
          if (duration_cast<milliseconds>(timeNow - *time).count() > 5000 ) {
            timesliceIndex.bookTimeslice(i);
            i.value++;
//            LOG(INFO) << "i.value " << i.value << ", timeNow " << duration_cast<milliseconds>(steady_clock::now().time_since_epoch()).count();
            *time = timeNow;
          }
        };

        ictx.services().get<CallbackService>().set(CallbackService::Id::ClockTick, bookTimesliceCallback);
        return (AlgorithmSpec::ProcessCallback) [sum](ProcessingContext& ctx) {

          if (ctx.inputs().isValid("data1") && ctx.inputs().isValid("data2")) {
            LOG(INFO) << "data input: " << ctx.inputs().get<int>("data1") << " " << ctx.inputs().get<int>("data2");
            *sum += ctx.inputs().get<int>("data1");
            *sum += ctx.inputs().get<int>("data2");
          }

          if (ctx.inputs().isValid("timer")) {
            LOG(INFO) << "timer input";
            ctx.outputs().snapshot<int>(Output{ "TST", "SUM" }, *sum);
          }

          LOG(INFO) << "sum: " << *sum;
        };
      }
    },
    Options{
      ConfigParamSpec{ "period-TST-TIMER", VariantType::Int, 5000, { "timer period" }}
    }
  };

  DataProcessorSpec printer{
    "printer",
    Inputs{
      { "sum", "TST", "SUM" }
    },
    Outputs{},
    AlgorithmSpec{
      (AlgorithmSpec::InitCallback) [](InitContext& ictx) {

        auto time = std::make_shared<steady_clock::time_point>(steady_clock::now());
        return (AlgorithmSpec::ProcessCallback) [time](ProcessingContext& ctx) {
          auto timeNow = steady_clock::now();
          LOG(INFO) << "The sum is: " << ctx.inputs().get<int>("sum")
                    << ", time interval of " << duration_cast<milliseconds>(timeNow - *time).count() / 1000.0 << "s";
          *time = timeNow;
        };
      }
    }
  };

  DataProcessorSpec sink{
    "sink",
    Inputs{
      { "data1", "TST", "DATA1" },
      { "data2", "TST", "DATA2" }
    },
    Outputs{},
    AlgorithmSpec{
      (AlgorithmSpec::InitCallback) [](InitContext& ictx) {
        return (AlgorithmSpec::ProcessCallback) [](ProcessingContext& ctx) {

          LOG(INFO) << "=== sink ===";
          if (ctx.inputs().isValid("data1") ){
            auto d = ctx.inputs().get("data1");
            const auto* dpHeader = get<DataProcessingHeader*>(d.header);
            assert(dpHeader);
            LOG(INFO) << "The data1 is: " << ctx.inputs().get<int>("data1") << ", tid: " << dpHeader->startTime;
          }
          if (ctx.inputs().isValid("data2")) {
            auto d = ctx.inputs().get("data2");
            const auto* dpHeader = get<DataProcessingHeader*>(d.header);
            assert(dpHeader);
            LOG(INFO) << "The data2 is: " << ctx.inputs().get<int>("data2") << ", tid: " << dpHeader->startTime;
          }
        };
      }
    }
  };



  WorkflowSpec specs{
    producer1,
    producer2,
    accumulator,
    printer,
    sink
  };

  return specs;
}