#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <chrono>

#include <boost/program_options.hpp>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <PCG/pcg_random.hpp>
//#include <TCanvas.h>

//#include "Framework/DataSamplingConditionFactory.h"
//#include "Framework/DataRef.h"
//#include "Framework/DataProcessingHeader.h"
//#include "Headers/DataHeader.h"

//using namespace o2::framework;
//using namespace o2::header;
using namespace std::chrono;
namespace bpo = boost::program_options;

class PRNG
{
 public:
  PRNG(uint64_t seed, double fraction) : mSeed(seed), mFraction(fraction){};
  virtual ~PRNG() = default;

  virtual bool rand(uint64_t in) = 0;

 protected:
  uint64_t mSeed = 0;
  double mFraction = 0.0;
};

class PRNG_PCG : public PRNG
{
 public:
  PRNG_PCG(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction),
      mGenerator(seed),
      mCurrentTimesliceID(0),
      mThreshold(static_cast<uint32_t>(mFraction * std::numeric_limits<uint32_t>::max())){};
  ~PRNG_PCG() override = default;

  bool rand(uint64_t in)
  {
    int64_t diff = in - mCurrentTimesliceID;
    if (diff == -1) {
      return mLastDecision;
    } else if (diff < -1) {
      mGenerator.backstep(static_cast<uint64_t>(-diff));
    } else if (diff > 0) {
      mGenerator.advance(static_cast<uint64_t>(diff));
    }

    mLastDecision = mGenerator() < mThreshold;
    mCurrentTimesliceID = in + 1;
    return mLastDecision;
  }

 private:
  uint32_t mThreshold = 0;
  pcg32_fast mGenerator;
  uint64_t mCurrentTimesliceID = 0;
  bool mLastDecision = 0;
};

class PRNG_Hash1 : public PRNG
{
 public:
  PRNG_Hash1(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction),
      mThreshold(static_cast<uint64_t>(mFraction * std::numeric_limits<uint64_t>::max())){};
  ~PRNG_Hash1() override = default;

  bool rand(uint64_t in)
  {
    uint64_t x = in * mSeed;
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);

    return x < mThreshold;
  }

 private:
  uint64_t mThreshold = 0;
};

std::unique_ptr<PRNG> createPRNG(std::string classname, uint64_t seed, double fraction)
{
  if (classname == "PRNG_PCG") {
    return std::make_unique<PRNG_PCG>(seed, fraction);
  } else if (classname == "PRNG_Hash1") {
    return std::make_unique<PRNG_Hash1>(seed, fraction);
  } else {
    return nullptr;
  }
}

double testSpeed(std::string prngName, uint64_t N)
{
  // speed test
//    std::cout << "> Speed test starting" << std::endl;

//    size_t antiOptimizationDummy = 0;
//    steady_clock::time_point timeStart = steady_clock::now();
//
//    for (size_t i = 0; i < testSize * 10; i++) {
//      DataProcessingHeader dph{ i, 0 };
//      o2::header::Stack headerStack{ dph };
//      DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };
//      antiOptimizationDummy += cond->decide(dr);
//    }
//
//    size_t nsecPerCall = duration_cast<nanoseconds>(steady_clock::now() - timeStart).count() / (testSize * 10);

//    std::cout << "Speed test: " << nsecPerCall << " ns/loop" << std::endl;
//    std::cout << "Ignore me " << antiOptimizationDummy << std::endl;


  return 0; //todo
}

double testSimpleChiSquare(std::string prngName, std::vector<uint64_t> seeds, double fraction, uint64_t N)
{
  std::vector<double> pValues;

  for (const auto seed : seeds) {

    std::unique_ptr<PRNG> prng = createPRNG(prngName, seed, fraction);
    uint64_t trues = 0;

    for (uint64_t n = 0; n < N; n++) {
      trues += prng->rand(n);
    }

    double chiSquare = 2 * std::pow(N * fraction - trues, 2) / (N * fraction);
    pValues.push_back(TMath::Prob(chiSquare, 1));
  }

  std::sort(pValues.begin(), pValues.end());

  std::vector<double> linearDistribution(seeds.size());
  for (size_t i = 0; i < linearDistribution.size(); i++) {
    linearDistribution[i] = i / (double)(linearDistribution.size() - 1);
  }
  double result = TMath::KolmogorovTest(pValues.size(), pValues.data(), linearDistribution.size(), linearDistribution.data(), "");

  return result;
}

double testRunsChiSquare(std::string prngName, std::vector<uint64_t> seeds, double fraction, uint64_t N)
{
  std::vector<double> pValues;
  const int buckets = 50;
  size_t runsHistoXmax = static_cast<size_t>(std::ceil(0.1 / fraction) * buckets) ;

  // integral from 0 to inf over p (1-p)^x / log(1-p) = -p / log(1 - p)
  // A * -p / log(1-p) = N * p   -----> A = -N * log(1-p)
  // this is to scale the function with relation to the number of samples (A) and buckets size (runsHistoXmax / buckets)
  std::string n = std::to_string((-1) * TMath::Log(1 - fraction) * N * runsHistoXmax / buckets);
  std::string p = std::to_string(fraction);

  // Probability mass function of Negative binomial distribution with 1 - p instead of p, r = 1, k = x
  TF1 ideal("ideal", (n + " * TMath::Power( 1 - " + p + ", x) * " + p).c_str(), 0, runsHistoXmax);

  Float_t xbins[buckets + 1] = {0};
  for (size_t i = 0; i < buckets + 1; i++) {
    xbins[i] = i * runsHistoXmax / buckets;
  }

  for (const auto seed : seeds) {

    std::unique_ptr<PRNG> prng = createPRNG(prngName, seed, fraction);
//    TH1I runsHisto("runs histo", "runs histo", buckets, 0, runsHistoXmax);
    TH1I runsHisto("runs histo", "runs histo", buckets, xbins);

    uint64_t runLength = 0;
    for (uint64_t n = 0; n < N; n++) { //runsHisto.GetSum()
      if (prng->rand(n)) {
        runsHisto.Fill(runLength);
        runLength = 0;
      } else {
        runLength++;
      }
    }

//    TCanvas *c1 = new TCanvas("c1","c1",1800,1000);
//    runsHisto.Scale(1);
//    runsHisto.Draw("same");
//    ideal.Draw("same");
//    c1->Update();
//    c1->Print("runs.png");
//    delete c1;


    double chiSquare = runsHisto.Chisquare(&ideal, "L");
    pValues.push_back(TMath::Prob(chiSquare, buckets - 1));
  }

  std::sort(pValues.begin(), pValues.end());

  std::vector<double> linearDistribution(seeds.size());
  for (size_t i = 0; i < linearDistribution.size(); i++) {
    linearDistribution[i] = i / (double)(linearDistribution.size() - 1);
  }
  double result = TMath::KolmogorovTest(pValues.size(), pValues.data(), linearDistribution.size(), linearDistribution.data(), "");

  return result;
}

double testMatricesRanks(std::string prngName, std::vector<uint64_t> seeds, double fraction, uint64_t N, size_t dimM)
{
  std::vector<double> pValues;
  bool matrix[dimM][dimM] = {0};

  // todo ideal histogram

  for (const auto seed : seeds) {

    std::unique_ptr<PRNG> prng = createPRNG(prngName, seed, fraction);

    for (size_t m = 0; m < N; m++) {
      size_t n = 0;
      for (size_t x = 0; x < dimM; x++) {
        for (size_t y = 0; y < dimM; y++) {
          matrix[x][y] = prng->rand(n++);
        }
      }
      // todo rank
      double rank = dimM;

      // todo put to histogram

    }

    // todo compare histogram with the ideal
    double chiSquare = 0;

    pValues.push_back(TMath::Prob(chiSquare, 1));
  }

  std::sort(pValues.begin(), pValues.end());

  std::vector<double> linearDistribution(seeds.size());
  for (size_t i = 0; i < linearDistribution.size(); i++) {
    linearDistribution[i] = i / (double)(linearDistribution.size() - 1);
  }
  double result = TMath::KolmogorovTest(pValues.size(), pValues.data(), linearDistribution.size(), linearDistribution.data(), "");

  return result;
}

std::vector<double> fractions = { 0.1, 0.05, 0.01, 0.001, 0.0001 };
//std::vector<double> fractions = { 0.2 };
std::vector<uint64_t> seeds = {
  552289370, 808574509, 214589538, 295179512, 863436240,
  19985442, 402171838, 848812329, 341858571, 191085231,
  819214847, 277661035, 239195551, 509038380, 639543200,
  469553710, 430692463, 972677980, 275793840, 60725136,
  387270561, 792236651, 455508846, 835767893, 741764246,
  978096680, 735658615, 356271575, 559210174, 403252930,
  810372976, 488518630, 252141442, 132143497, 291694362,
  376096303, 711357754, 441666606, 851409926, 830976745,
  93245143, 263517995, 297148760, 731715905, 581928278,
  756284286, 338778058, 840724372, 600526144, 948152916,
  566763840, 811020761, 120018113, 231918662, 417197404,
  21084436, 334675356, 961437177, 241392364, 83520946,
  462322598, 404773570, 525235048, 458286100, 669228482,
  126816757, 584579025, 528267678, 733341089, 545152755,
  994814077, 645667987, 525528174, 353815683, 33435548,
  289834496, 533561188, 182664687, 203100193, 614725360,
  439852892, 98433473, 355483966, 954706423, 354224002,
  449382088, 874226835, 750969153, 173745243, 468339785,
  963621706, 542034923, 824278118, 548449232, 353197413,
  358202811, 992338112, 919798241, 669107447, 125556072
};

int main(int argc, char* argv[])
{

  // Arguments parsing
  bpo::variables_map vm;
  bpo::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Produce help message.")
    ("prng,p", bpo::value<std::string>()->default_value("PRNG_PCG"), R"(prng version, "PRNG_PCG", "PRNG_Hash1".)")
    ("mode,m", bpo::value<std::string>()->default_value("test"), R"(prng version, "produce", "test".)")
    ("test-size,s", bpo::value<size_t>()->default_value(10000000));
  bpo::store(parse_command_line(argc, argv, desc), vm);
  bpo::notify(vm);

  std::string prng = vm["prng"].as<std::string>();

  if (vm["mode"].as<std::string>() == "produce") {
    /*
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
    */
  } else if (vm["mode"].as<std::string>() == "test") {

    std::cout << "================================================\n";
    std::cout << "Starting non-uniform random bit generation tests\n";
    std::cout << "================================================\n";
    std::cout << "Chosen PRNG: " << prng << "\n";

//    auto cond = DataSamplingConditionFactory::create(prng);
//    auto prng = createPRNG(prng);
    size_t testSize = vm["test-size"].as<size_t>();

  size_t nsecPerCall = testSpeed(prng, testSize);


    std::vector<double> chiSquareTest;
    std::vector<double> runsChiSquareTest;
    std::vector<double> runsRmsTest;

    // -- do it for different percentages --
    for (const auto fraction : fractions) {
      std::cout << "> Testing for fraction " << fraction << std::endl;

      chiSquareTest.push_back(testSimpleChiSquare(prng, seeds, fraction, testSize));
      runsChiSquareTest.push_back(testRunsChiSquare(prng, seeds, fraction, testSize ));

    }

    std::cout << "============================= TEST RESULTS =============================" << std::endl;
    std::cout << "Method: " << prng << std::endl;
    std::cout << "Speed: " << nsecPerCall << " ns/loop" << std::endl;
    std::cout << "Tested fracions | ";
    for (auto fr : fractions)
      std::cout << std::setw(10) << fr << " ";
    std::cout << std::endl;
    std::cout << "------------------------------------------------------------------------" << std::endl;
    std::cout << "Chi square      | ";
    for (auto chi : chiSquareTest)
      std::cout << std::setw(10) << chi << " ";
    std::cout << std::endl;
    std::cout << "Runs Histo RMS  | ";
    for (auto rms : runsRmsTest)
      std::cout << std::setw(10) << rms << " ";
    std::cout << std::endl;
    std::cout << "Runs Histo Chi2 | ";
    for (auto chi : runsChiSquareTest)
      std::cout << std::setw(10) << chi << " ";
    std::cout << std::endl;
  }

  return 0;
}