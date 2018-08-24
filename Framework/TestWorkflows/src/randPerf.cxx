#include <iostream>
#include <random>

#include <TRandom.h>
#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>


#include <iostream>
#include <vector>
//#include <boost/test/unit_test.hpp>

#include "Framework/DataSamplingConditionFactory.h"
#include "Framework/DataRef.h"
#include "Framework/DataProcessingHeader.h"
#include "Headers/DataHeader.h"

using namespace o2::framework;
using namespace o2::header;
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <boost/program_options.hpp>
#include <TF1.h>

using namespace std::chrono;
namespace bpo = boost::program_options;


int main(int argc, char *argv[]) {

  // Arguments parsing
  bpo::variables_map vm;
  bpo::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Produce help message.")
    ("prng,p",      bpo::value<std::string>()->default_value("pcg"),    R"(prng version, "pcg", "hash", "hashCombine".)")
    ("mode,m",      bpo::value<std::string>()->default_value("produce"),    R"(prng version, "produce", "test".)")
    ("test-size,s",      bpo::value<size_t>()->default_value(10000000));
  bpo::store(parse_command_line(argc, argv, desc), vm);
  bpo::notify(vm);

  std::string prng = vm["prng"].as<std::string>();

  if ( vm["mode"].as<std::string>() == "produce") {
    auto cond = DataSamplingConditionFactory::create(prng);
    uint64_t i = 0;
    size_t bit = 0;
    unsigned char byte = 0;
    while (cond) {
      if (prng == "pcg" || prng.substr(0, 7) == "TRandom") {
        uint32_t rnd = cond->rnd(i);
//      std::cout << rnd << std::endl;
        write(1, &rnd, sizeof(uint32_t));
      } else if (prng == "hash") {
        bool rnd = static_cast<bool>(cond->rnd(i));
        byte |= (static_cast<unsigned char>(rnd) & 0x01) << bit++;
        if (bit >= 8) {
          write(1, &byte, sizeof(unsigned char));
//        std::cout << uint16_t(byte) << std::endl;
          bit = 0;
          byte = 0;
        }
      } else if (prng == "hashCombine") {
        uint64_t rnd = cond->rnd(i);
        write(1, &rnd, sizeof(uint64_t));
      }
      i++;
    }
  } else if (vm["mode"].as<std::string>() == "test") {

    LOG(INFO) << "================================================";
    LOG(INFO) << "Starting non-uniform random bit generation tests";
    LOG(INFO) << "Chosen PRNG: " << prng;

    auto cond = DataSamplingConditionFactory::create(prng);
    size_t testSize = vm["test-size"].as<size_t>();
    // speed test
    LOG(INFO) << "> Speed test starting";

    size_t antiOptimizationDummy = 0;
    steady_clock::time_point timeStart = steady_clock::now();

    for (size_t i = 0; i < testSize; i++) {
      DataProcessingHeader dph{ i, 0 };
      o2::header::Stack headerStack{ dph };
      DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };
      antiOptimizationDummy += cond->decide(dr);
    }

    size_t nsecPerCall = duration_cast<nanoseconds>(steady_clock::now() - timeStart).count() / testSize;

    LOG(INFO) << "Speed test: " << nsecPerCall << " ns/loop";
    LOG(INFO) << "Ignore me " << antiOptimizationDummy;

    std::vector<double> fractions = {0.1, 0.05, 0.01, 0.001, 0.0001};
    std::vector<uint64_t> seeds = {3845983573, 343045, 102130001, 9990234234, 13412404};

    std::vector<double> chiSquareTest;
    std::vector<double> histoChiSquareTest;
    // -- do it for different percentages --
    for (const auto fraction : fractions) {
      LOG(INFO) << "> Testing for fraction " << fraction;

      double chiSquareAll = 0;
      double histoChiSquareAll = 0;

      for (const auto seed : seeds) {

        boost::property_tree::ptree config;
        config.put("fraction", fraction);
        config.put("seed", seed);
        cond->configure(config);

        uint64_t trues = 0;
        uint64_t runLength = 0;
        TH1F runsHisto("runs histo", "runs histo", 100, 0, 1/fraction * 10);

        for (size_t i = 0; i < testSize; i++) {

          DataProcessingHeader dph{ i, 0 };
          o2::header::Stack headerStack{ dph };
          DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };

          bool decision = cond->decide(dr);
          trues += decision;
          // runs test
          if ( decision )  {
            runsHisto.Fill(runLength);
            runLength = 0;
          } else if (!decision){
            runLength++;
          }
        }

        double chiSquare = std::pow(testSize * fraction - trues, 2) / (testSize * fraction);

        TF1 ideal("exp", (std::to_string(fraction) + "*TMath::Exp(-" + std::to_string(fraction) + "*x)").c_str(), 0, 1/fraction * 10);
        runsHisto.Scale(1/runsHisto.Integral("width"));
//        TCanvas* c1 = new TCanvas();
//        ideal.Draw();
//        runsHisto.Draw("SAME");
//        c1->SaveAs("histo.png");
//        delete c1;

        LOG(INFO) << "chiSquare result: " << chiSquare;
        double histoChiSquare = runsHisto.Chisquare(&ideal);
        LOG(INFO) << "histoChiSquare result: " << histoChiSquare;

        chiSquareAll += chiSquare;
        histoChiSquareAll += histoChiSquare;
      }
      chiSquareTest.push_back(chiSquareAll);
      histoChiSquareTest.push_back(histoChiSquareAll);
    }

    std::cout << "===== TEST RESULTS =====" << std::endl;
    std::cout << "Method: " << prng << std::endl;
    std::cout << "Speed: " << nsecPerCall << " ns/loop" << std::endl;
    std::cout << "Tested fracions | ";
    for (auto fr : fractions) std::cout << std::setw(10) << fr << " ";
    std::cout << std::endl;
    std::cout << "------------------------------------------------------------------------" << std::endl;
    std::cout << "Chi square      | ";
    for (auto chi : chiSquareTest) std::cout << std::setw(10) << chi << " ";
    std::cout << std::endl;
    std::cout << "Histo Chi square| ";
    for (auto chi : histoChiSquareTest) std::cout << std::setw(10) << chi << " ";
    std::cout << std::endl;
  }

  return 0;
}