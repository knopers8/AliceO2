// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <chrono>

#include "Framework/CompletionPolicy.h"
#include "Framework/CompletionPolicyHelpers.h"
#include "Framework/ConfigParamSpec.h"

using namespace std::chrono;
using namespace o2::framework;

void customize(std::vector<CompletionPolicy>& policies)
{
  policies.push_back(CompletionPolicyHelpers::consumeWhenAny());
}

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{ "lhc-size", VariantType::Int, 5, { "The size of message chain" }});
  workflowOptions.push_back(ConfigParamSpec{ "bunches", VariantType::Int, 1, { "The number of messages to inject" }});
}

#include "Framework/DataProcessorSpec.h"
#include "Framework/runDataProcessing.h"
#include "Monitoring/MonitoringFactory.h"

using namespace o2::monitoring;

WorkflowSpec defineDataProcessing(ConfigContext const& cc)
{
  int lhcSize = cc.options().get<int>("lhc-size");
  int bunchesTotal = cc.options().get<int>("bunches");
  WorkflowSpec specs;

  for (size_t i = 0; i < lhcSize; i++) {
    DataProcessorSpec lhcSector{
      "lhc_sector" + std::to_string(i),
      Inputs{
        InputSpec{ "bunch-in", "LHC", "ENERGY", i }
      },
      {
        OutputSpec{{ "bunch-out" }, "LHC", "ENERGY", (i + 1) % lhcSize },
      },
      AlgorithmSpec{
        (AlgorithmSpec::InitCallback)[i, lhcSize](InitContext& ictx) {

//          auto& monitoring = ictx.services().get<Monitoring>();
//          std::string url = "influxdb-udp://aido2mon-gpn.cern.ch:8087";
//          monitoring.addBackend(MonitoringFactory::GetBackend(url));
//          monitoring.enableProcessMonitoring(1);
//          monitoring.addGlobalTag("taskName", "lhc_sector" + std::to_string(i));

          return (AlgorithmSpec::ProcessCallback) [i, lhcSize](ProcessingContext& ctx) {
            const int& bunchIn = ctx.inputs().get<int>("bunch-in");
            int& bunchOut = ctx.outputs().make<int>(Output{ "LHC", "ENERGY", (i + 1) % lhcSize });
            bunchOut = bunchIn;
          };
        }
      }
    };
    specs.push_back(lhcSector);
  };

  // special lhc sector with injector and statistics log
  specs[0].inputs.push_back(InputSpec{ "bunch-sps", "SPS", "ENERGY" });
  specs[0].algorithm = AlgorithmSpec{
    (AlgorithmSpec::InitCallback)[bunchesTotal](InitContext& ictx) {

      auto& monitoring = ictx.services().get<Monitoring>();
//      std::string url = "influxdb-udp://aido2mon-gpn.cern.ch:8087";
//      monitoring.addBackend(MonitoringFactory::GetBackend(url));
      monitoring.enableProcessMonitoring(1);
      monitoring.addGlobalTag("taskName", "lhc_sector0");

      return (AlgorithmSpec::ProcessCallback) [bunchesTotal](ProcessingContext& ctx) {

        static steady_clock::time_point periodStart = steady_clock::now();
        static long long bunchesPassed = 0;

        auto timeDiff = duration_cast<microseconds>(steady_clock::now() - periodStart).count();
        if (timeDiff > 1000000) {
          double cycles = bunchesPassed * 1000000.0 / timeDiff / bunchesTotal;
          LOG(INFO) << cycles << " cycles per second";
          periodStart = steady_clock::now();
          bunchesPassed = 0;

          auto& monitoring = ctx.services().get<Monitoring>();
//          std::string url = "influxdb-udp://aido2mon-gpn.cern.ch:8087";
//          monitoring.addBackend(MonitoringFactory::GetBackend(url));
//          monitoring.enableProcessMonitoring(1);
//          monitoring.addGlobalTag("taskName", "lhc_sector0");

          monitoring.send({cycles, "cyclesPerSecond"});
        }

        if (ctx.inputs().isValid("bunch-in")) {
          const int& bunchIn = ctx.inputs().get<int>("bunch-in");
          int& bunchOut = ctx.outputs().make<int>(Output{ "LHC", "ENERGY", 1 });
          bunchOut = bunchIn + 1;
        } else { // bunch-sps valid
          ctx.outputs().make<int>(OutputRef{ "bunch-out", 1 }, 1);
        }
        bunchesPassed++;
      };
    }
  };

  DataProcessorSpec sps{
    "sps-injector",
    Inputs{},
    Outputs{
      OutputSpec{{ "bunch-sps" }, "SPS", "ENERGY" }
    },
    AlgorithmSpec{
      (AlgorithmSpec::InitCallback) [bunchesTotal](InitContext&) {

        std::shared_ptr<int> sent = std::make_shared<int>(0);

        return (AlgorithmSpec::ProcessCallback) [sent, bunchesTotal](ProcessingContext& ctx) mutable {
          if (*sent < bunchesTotal) {
            LOG(INFO) << "INJECTION";

            auto bunch = ctx.outputs().make<int>(OutputRef{ "bunch-sps" }, 1);
            (*sent)++;
          } else {
            sleep(1);
          }
        };
      }
    }
  };

  specs.push_back(sps);
  return specs;
}
