#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <chrono>

#include <boost/program_options.hpp>
#include <TF1.h>
#include <TH1F.h>
//#include <TCanvas.h>

#include "Framework/DataSamplingConditionFactory.h"
#include "Framework/DataRef.h"
#include "Framework/DataProcessingHeader.h"
#include "Headers/DataHeader.h"

using namespace o2::framework;
using namespace o2::header;
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
    LOG(INFO) << "================================================";
    LOG(INFO) << "Chosen PRNG: " << prng;

    auto cond = DataSamplingConditionFactory::create(prng);
    size_t testSize = vm["test-size"].as<size_t>();
    // speed test
    LOG(INFO) << "> Speed test starting";

    size_t antiOptimizationDummy = 0;
    steady_clock::time_point timeStart = steady_clock::now();

    for (size_t i = 0; i < testSize * 10; i++) {
      DataProcessingHeader dph{ i, 0 };
      o2::header::Stack headerStack{ dph };
      DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };
      antiOptimizationDummy += cond->decide(dr);
    }

    size_t nsecPerCall = duration_cast<nanoseconds>(steady_clock::now() - timeStart).count() / (testSize * 10);

    LOG(INFO) << "Speed test: " << nsecPerCall << " ns/loop";
    LOG(INFO) << "Ignore me " << antiOptimizationDummy;

    std::vector<double> fractions = {0.1, 0.05, 0.01, 0.001, 0.0001};
    std::vector<uint64_t> seeds = {
      552289370, 808574509, 214589538, 295179512, 863436240,
      19985442,  402171838, 848812329, 341858571, 191085231,
      819214847, 277661035, 239195551, 509038380, 639543200,
      469553710, 430692463, 972677980, 275793840, 60725136,
      387270561, 792236651, 455508846, 835767893, 741764246,
      978096680, 735658615, 356271575, 559210174, 403252930,
      810372976, 488518630, 252141442, 132143497, 291694362,
      376096303, 711357754, 441666606, 851409926, 830976745,
      93245143,  263517995, 297148760, 731715905, 581928278,
      756284286, 338778058, 840724372, 600526144, 948152916,
      566763840, 811020761, 120018113, 231918662, 417197404,
      21084436,  334675356, 961437177, 241392364, 83520946,
      462322598, 404773570, 525235048, 458286100, 669228482,
      126816757, 584579025, 528267678, 733341089, 545152755,
      994814077, 645667987, 525528174, 353815683, 33435548,
      289834496, 533561188, 182664687, 203100193, 614725360,
      439852892, 98433473,  355483966, 954706423, 354224002,
      449382088, 874226835, 750969153, 173745243, 468339785,
      963621706, 542034923, 824278118, 548449232, 353197413,
      358202811, 992338112, 919798241, 669107447, 125556072
    };

    std::vector<double> chiSquareTest;
    std::vector<double> histoChiSquareTest;
    std::vector<double> runsRmsTest;

    std::vector<double> linearDistribution(seeds.size());
    for(size_t i = 0; i < linearDistribution.size(); i++) {
      linearDistribution[i] = i / (double)(linearDistribution.size()-1);
    }

    // -- do it for different percentages --
    for (const auto fraction : fractions) {
      LOG(INFO) << "> Testing for fraction " << fraction;
      size_t runsHistoXmax = 1/fraction * 10;
      std::vector<double> pValues;
//      double chiSquareAll = 0;
      double histoChiSquareAll = 0;
      double runsRmsSum = 0;

      for (const auto seed : seeds) {

        boost::property_tree::ptree config;
        config.put("fraction", fraction);
        config.put("seed", seed);
        cond->configure(config);

        uint64_t trues = 0;
        uint64_t runLength = 0;
        TH1F runsHisto("runs histo", "runs histo", 100, 0, runsHistoXmax);

        for (size_t i = 0; i < testSize; i++) {

          DataProcessingHeader dph{ i, 0 };
          o2::header::Stack headerStack{ dph };
          DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };

          bool decision = cond->decide(dr);
          trues += decision;
          // runs test
          if ( decision ) {
            runsHisto.Fill(runLength);
            runLength = 0;
          } else if (!decision) {
            runLength++;
          }
        }

        double chiSquare = 2*std::pow(testSize * fraction - trues, 2) / (testSize * fraction);
        pValues.push_back(TMath::Prob(chiSquare, 1));

        TH1F runsIdealHisto("runs ideal histo", "runs ideal histo", 100, 0, runsHistoXmax);
        for (int i = 0 ; i < 100; i++) {
          double x = runsHistoXmax * i / 99;
          int times = fraction*TMath::Exp(-fraction*x)*testSize;
          for ( int j=0; j < times; j++) {
            runsIdealHisto.Fill(x);
          }
        }
        runsIdealHisto.Scale(1/runsIdealHisto.Integral("width"));

        TF1 ideal("exp", (std::to_string(fraction) + "*TMath::Exp(-" + std::to_string(fraction) + "*x)").c_str(), 0, runsHistoXmax);
        runsHisto.Scale(1/runsHisto.Integral("width"));
//        TCanvas* c1 = new TCanvas();
//        ideal.Draw();
//        runsHisto.Draw("SAME");
//        c1->SaveAs("histo.png");
//        delete c1;

//        LOG(INFO) << "chiSquare result: " << chiSquare << " pval: " << TMath::Prob(chiSquare, 1);
        double histoChiSquare = runsHisto.Chisquare(&ideal);
        TH1* asymmetry = runsHisto.GetAsymmetry(&runsIdealHisto);
        asymmetry->SetDefaultSumw2(true);
        double asymmetryRMS = asymmetry->GetRMS();

        std::cout << histoChiSquare << " " << asymmetryRMS << std::endl;
        delete asymmetry;

//        chiSquareAll += chiSquare;
        histoChiSquareAll += histoChiSquare;
        runsRmsSum += asymmetryRMS;
      }

      std::sort(pValues.begin(), pValues.end());

      double kmg = TMath::KolmogorovTest(pValues.size(), pValues.data(), linearDistribution.size(), linearDistribution.data(), "");
      LOG(INFO) << "TMath::KolmogorovTest result: " << kmg;
//      chiSquareTest.push_back(chiSquareAll);
      chiSquareTest.push_back(kmg);
      histoChiSquareTest.push_back(histoChiSquareAll);
      runsRmsTest.push_back(runsRmsSum / seeds.size());
    }

    std::cout << "=================== TEST RESULTS ===================" << std::endl;
    std::cout << "Method: " << prng << std::endl;
    std::cout << "Speed: " << nsecPerCall << " ns/loop" << std::endl;
    std::cout << "Tested fracions | ";
    for (auto fr : fractions) std::cout << std::setw(10) << fr << " ";
    std::cout << std::endl;
    std::cout << "------------------------------------------------------------------------" << std::endl;
    std::cout << "Chi square      | ";
    for (auto chi : chiSquareTest) std::cout << std::setw(10) << chi << " ";
    std::cout << std::endl;
    std::cout << "Runs Histo RMS  | ";
    for (auto rms : runsRmsTest) std::cout << std::setw(10) << rms << " ";
    std::cout << std::endl;
    std::cout << "Runs Histo Chi2 | ";
    for (auto chi : histoChiSquareTest) std::cout << std::setw(10) << chi << " ";
    std::cout << std::endl;
  }

  return 0;
}