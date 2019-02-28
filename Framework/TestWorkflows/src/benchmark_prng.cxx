#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <fstream>

#include <boost/program_options.hpp>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TRandomGen.h>
#include <TVirtualFFT.h>
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
  virtual uint32_t fullRand(uint64_t in) { return 0; };

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
  uint64_t mThreshold;
};

class PRNG_Hash2 : public PRNG
{
 public:
  PRNG_Hash2(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction),
      mThreshold(static_cast<uint64_t>(mFraction * std::numeric_limits<uint64_t>::max()))
  {
    mSeed = hashCombine(0, seed);
  };
  ~PRNG_Hash2() override = default;

  uint64_t hashCombine(uint64_t h, uint64_t k)
  {
    // implementation copied from boost/functional/hash.hpp because it is required that it does not change with future
    // boost versions (to guarantee determinism in many runs and machines)

    const uint64_t m = UINT64_C(0xc6a4a7935bd1e995);
    const int r = 47;

    k *= m;
    k ^= k >> r;
    k *= m;

    h ^= k;
    h *= m;

    // Completely arbitrary number, to prevent 0's
    // from hashing to 0.
    h += 0xe6546b64;

    return h;
  }

  bool rand(uint64_t in)
  {
    return hashCombine(mSeed, in) < mThreshold;
  }

 private:
  uint64_t mThreshold;
};

class PRNG_TRandom : public PRNG
{
 public:
  PRNG_TRandom(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction){};
  ~PRNG_TRandom() override = default;

  bool rand(uint64_t in)
  {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<bool>(mGenerator.Binomial(1, mFraction));
  }

  uint32_t fullRand(uint64_t in) {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<uint32_t>(mGenerator.Rndm() * std::numeric_limits<uint32_t>::max());
  };

 private:
  TRandom mGenerator;
};

class PRNG_TRandom1 : public PRNG
{
 public:
  PRNG_TRandom1(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction){};
  ~PRNG_TRandom1() override = default;

  bool rand(uint64_t in)
  {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<bool>(mGenerator.Binomial(1, mFraction));
  }

  uint32_t fullRand(uint64_t in) {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<uint32_t>(mGenerator.Rndm() * std::numeric_limits<uint32_t>::max());
  };

  private:
  TRandom1 mGenerator;
};

class PRNG_TRandom2 : public PRNG
{
 public:
  PRNG_TRandom2(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction){};
  ~PRNG_TRandom2() override = default;

  bool rand(uint64_t in)
  {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<bool>(mGenerator.Binomial(1, mFraction));
  }

  uint32_t fullRand(uint64_t in) {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<uint32_t>(mGenerator.Rndm() * std::numeric_limits<uint32_t>::max());
  };

 private:
  TRandom2 mGenerator;
};

class PRNG_TRandom3 : public PRNG
{
 public:
  PRNG_TRandom3(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction){};
  ~PRNG_TRandom3() override = default;

  bool rand(uint64_t in)
  {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<bool>(mGenerator.Binomial(1, mFraction));
  }

  uint32_t fullRand(uint64_t in) {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<uint32_t>(mGenerator.Rndm() * std::numeric_limits<uint32_t>::max());
  };

 private:
  TRandom3 mGenerator;
};

class PRNG_TRandomMixMax : public PRNG
{
 public:
  PRNG_TRandomMixMax(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction){};
  ~PRNG_TRandomMixMax() override = default;

  bool rand(uint64_t in)
  {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<bool>(mGenerator.Binomial(1, mFraction));
  }

  uint32_t fullRand(uint64_t in) {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<uint32_t>(mGenerator.Rndm() * std::numeric_limits<uint32_t>::max());
  };

 private:
  TRandomMixMax mGenerator;
};

class PRNG_TRandomMixMax17 : public PRNG
{
 public:
  PRNG_TRandomMixMax17(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction){};
  ~PRNG_TRandomMixMax17() override = default;

  bool rand(uint64_t in)
  {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<bool>(mGenerator.Binomial(1, mFraction));
  }

  uint32_t fullRand(uint64_t in) {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<uint32_t>(mGenerator.Rndm() * std::numeric_limits<uint32_t>::max());
  };

 private:
  TRandomMixMax17 mGenerator;
};

class PRNG_TRandomMT64 : public PRNG
{
 public:
  PRNG_TRandomMT64(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction){};
  ~PRNG_TRandomMT64() override = default;

  bool rand(uint64_t in)
  {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<bool>(mGenerator.Binomial(1, mFraction));
  }

  uint32_t fullRand(uint64_t in) {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<uint32_t>(mGenerator.Rndm() * std::numeric_limits<uint32_t>::max());
  };

 private:
  TRandomMT64 mGenerator;
};

class PRNG_TRandomRanlux48 : public PRNG
{
 public:
  PRNG_TRandomRanlux48(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction){};
  ~PRNG_TRandomRanlux48() override = default;

  bool rand(uint64_t in)
  {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<bool>(mGenerator.Binomial(1, mFraction));
  }

  uint32_t fullRand(uint64_t in) {
    mGenerator.SetSeed(in * mSeed);
    return static_cast<uint32_t>(mGenerator.Rndm() * std::numeric_limits<uint32_t>::max());
  };

 private:
  TRandomRanlux48 mGenerator;
};

class PRNG_Dummy : public PRNG
{
  public:
  PRNG_Dummy(uint64_t seed, double fraction)
    : PRNG::PRNG(seed, fraction){};
  ~PRNG_Dummy() override = default;

  bool rand(uint64_t in)
  {
    return true;
  }
};

std::unique_ptr<PRNG> createPRNG(std::string classname, uint64_t seed, double fraction)
{
  if (classname == "PRNG_PCG") {
    return std::make_unique<PRNG_PCG>(seed, fraction);
  } else if (classname == "PRNG_Hash1") {
    return std::make_unique<PRNG_Hash1>(seed, fraction);
  } else if (classname == "PRNG_Hash2") {
    return std::make_unique<PRNG_Hash2>(seed, fraction);
  } else if (classname == "PRNG_TRandom") {
    return std::make_unique<PRNG_TRandom>(seed, fraction);
  } else if (classname == "PRNG_TRandom1") {
    return std::make_unique<PRNG_TRandom1>(seed, fraction);
  } else if (classname == "PRNG_TRandom2") {
    return std::make_unique<PRNG_TRandom2>(seed, fraction);
  } else if (classname == "PRNG_TRandom3") {
    return std::make_unique<PRNG_TRandom3>(seed, fraction);
  } else if (classname == "PRNG_TRandomMixMax") {
    return std::make_unique<PRNG_TRandomMixMax>(seed, fraction);
  } else if (classname == "PRNG_TRandomMixMax17") {
    return std::make_unique<PRNG_TRandomMixMax17>(seed, fraction);
  } else if (classname == "PRNG_TRandomMT64") {
    return std::make_unique<PRNG_TRandomMT64>(seed, fraction);
  } else if (classname == "PRNG_TRandomRanlux48") {
    return std::make_unique<PRNG_TRandomRanlux48>(seed, fraction);
  } else if (classname == "PRNG_Dummy") {
    return std::make_unique<PRNG_Dummy>(seed, fraction);
  } else {
    throw std::runtime_error("Unknown PRNG: " + classname);
  }
}

double testSpeed(std::string prngName, uint64_t N)
{
  // speed test
  std::cout << "> Speed test starting for " << prngName << std::endl;
  std::unique_ptr<PRNG> prng = createPRNG(prngName, 397583947, 0.01);
  size_t antiOptimizationDummy = 0;
  steady_clock::time_point timeStart = steady_clock::now();

  for (size_t i = 0; i < N; i++) {
    antiOptimizationDummy += prng->rand(i);
  }

  double nsPerCall = duration_cast<nanoseconds>(steady_clock::now() - timeStart).count() / double(N);

  std::cout << "Speed test: " << nsPerCall << " ns/loop" << std::endl;
  std::cout << "Ignore me " << antiOptimizationDummy << std::endl;

  return nsPerCall;
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
  size_t runsHistoXmax = static_cast<size_t>(std::ceil(0.1 / fraction) * buckets);

  // integral from 0 to inf over p (1-p)^x / log(1-p) = -p / log(1 - p)
  // A * -p / log(1-p) = N * p   -----> A = -N * log(1-p)
  // this is to scale the function with relation to the number of samples (A) and buckets size (runsHistoXmax / buckets)
  std::string A = std::to_string((-1) * TMath::Log(1 - fraction) * N * runsHistoXmax / buckets);
  std::string p = std::to_string(fraction);

  // Probability mass function of Negative binomial distribution with 1 - p instead of p, r = 1, k = x
  TF1 ideal("ideal", (A + " * TMath::Power( 1 - " + p + ", x) * " + p).c_str(), 0, runsHistoXmax);

  Float_t xbins[buckets + 1] = { 0 };
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

//    TCanvas* c1 = new TCanvas("c1", "c1", 1800 / 2, 1000 / 2);
    //    runsHisto.Scale(1);
//    runsHisto.Draw("same");
//    ideal.Draw("same");
//    c1->Update();
//    c1->Print("runs.png");
//    delete c1;

    double chiSquare = runsHisto.Chisquare(&ideal);
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

void swap(std::vector<std::vector<int>>& mat, int row1, int row2, int col)
{
  for (int i = 0; i < col; i++) {
    int temp = mat[row1][i];
    mat[row1][i] = mat[row2][i];
    mat[row2][i] = temp;
  }
}

/* function for finding rank of matrix */
// from https://www.geeksforgeeks.org/program-for-rank-of-matrix/
int rankOfMatrix(std::vector<std::vector<int>>& mat, size_t dim)
{
  int rank = dim;

  for (int row = 0; row < rank; row++) {
    if (mat[row][row]) {
      for (int col = 0; col < dim; col++) {
        if (col != row) {
          double mult = (double)mat[col][row] / mat[row][row];
          for (int i = 0; i < rank; i++)
            mat[col][i] -= mult * mat[row][i];
        }
      }
    }
    else {
      bool reduce = true;
      for (int i = row + 1; i < dim; i++) {
        if (mat[i][row]) {
          swap(mat, row, i, rank);
          reduce = false;
          break;
        }
      }
      if (reduce) {
        rank--;
        for (int i = 0; i < dim; i++)
          mat[i][row] = mat[i][rank];
      }
      row--;
    }
  }
  return rank;
}

double testMatricesRanks(std::string prngName, std::vector<uint64_t> seeds, double fraction, uint64_t N, size_t dimM)
{
  std::vector<double> pValues;
  //  int matrix[dimM][dimM] = {0};

  std::vector<std::vector<int>> matrix;
  matrix.resize(dimM);
  for (auto& it : matrix) {
    it.resize(dimM);
  }

  // todo ideal histogram
  auto product = [](int j0, int jn, double t) {
    double prod = 1;
    for (int j = j0; j <= jn; j++) {
      prod *= 1 - 1 / pow(t, j);
    }
    return prod;
  };

  auto ideal = [N, &product, dimM, t=1/(1-fraction)](double *x, double *p) {
    double k = dimM - std::floor(*x);
    return N * product(k + 1, 10000, t) / product(1, k, t) * std::pow(1 / t, k*k);
  };

  TF1 idealf("idealf", ideal, 1, dimM, 1);

  for (const auto seed : seeds) {

    TH1I ranksHisto("ranks histo", "ranks histo", dimM, 0, dimM);

    std::unique_ptr<PRNG> prng = createPRNG(prngName, seed, fraction);

    size_t n = 0;
    for (size_t m = 0; m < N; m++) {
      for (size_t x = 0; x < dimM; x++) {
        for (size_t y = 0; y < dimM; y++) {
          matrix[x][y] = prng->rand(n++);
        }
      }
      double rank = rankOfMatrix(matrix, dimM);
      ranksHisto.Fill(rank);
    }

//    TCanvas* c1 = new TCanvas("c1", "c1", 1800 / 2, 1000 / 2);
////    ranksHisto.Scale(1.0 / ranksHisto.Integral());
//    idealf.Draw();
//    ranksHisto.Draw("same");
//    c1->Update();
//    c1->Print("ranks.png");
//    delete c1;

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

double testRunsAutoCorrelation(std::string prngName, std::vector<uint64_t> seeds, double fraction, uint64_t N)
{
  std::vector<double> stddevs;

  TF1 ideal("ideal", "1", 0, N);
//  const size_t reductionFactor = 1024*16;

  for (const auto seed : seeds) {

    std::unique_ptr<PRNG> prng = createPRNG(prngName, seed, fraction);
    //    TH1I runsHisto("runs histo", "runs histo", buckets, 0, runsHistoXmax);
    TH1I runs("runs", "runs", N, 0, N);

    uint64_t runLength = 0;
    uint64_t j = 0;
    for (uint64_t n = 0; n < N; ) {
      if (prng->rand(j++)) {
        runs.SetBinContent(n+1, runLength);
        runLength = 0;
        n++;
      } else {
        runLength++;
      }
    }

    TCanvas* c1 = new TCanvas("c1", "c1", 1800, 1000 );
    //Compute the transform and look at the magnitude of the output
    TH1 *hm =0;
    TVirtualFFT::SetTransform(0);
    hm = runs.FFT(hm, "MAG");
    if (!hm) {
      throw std::runtime_error("hm null");
    }
    hm->SetBinContent(1, 0); // clear the constant
    hm->Scale(1 / TMath::Sqrt(N));
//    hm->Rebin(reductionFactor);
//    hm->Scale(1 / (TMath::Sqrt(N) * reductionFactor));
//    hm->SetTitle("Magnitude of the 1st transform");
//    hm->Draw();
//    c1->Update();
//    c1->Print("runs.png");

    delete c1;

    double mean = hm->Integral() / (N );
    double stddev = 0;
    for (int i = 0; i < N ; i++) {
      stddev += std::pow(mean - hm->GetBinContent(i + 1), 2) / mean;
    }
    stddev = TMath::Sqrt(stddev);

    delete hm;
    stddevs.push_back(stddev);
  }

  double result = 0;
  result = std::accumulate(stddevs.begin(), stddevs.end(), 0) / stddevs.size();

  return result;
}

std::vector<double> fractions = { 0.1, 0.05, 0.01, 0.001, 0.0001 };
//std::vector<double> fractions = { 0.1 };
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
  desc.add_options()("help,h", "Produce help message.")("prng,p", bpo::value<std::string>()->default_value("PRNG_PCG"), R"(prng version, "PRNG_PCG", "PRNG_Hash1".)")("mode,m", bpo::value<std::string>()->default_value("test"), R"(prng version, "produce", "test".)")("test-size,s", bpo::value<size_t>()->default_value(10000000));
  bpo::store(parse_command_line(argc, argv, desc), vm);
  bpo::notify(vm);



  if (vm["mode"].as<std::string>() == "produce") {

    std::unique_ptr<PRNG> prng = createPRNG(vm["prng"].as<std::string>(), 1928472385, 0.5);

    uint64_t i = 0;
    size_t bit = 0;
    unsigned char byte = 0;
    while (prng) {
      uint32_t rnd = prng->fullRand(i++);
      write(1, &rnd, sizeof(uint32_t));
    }

  } else if (vm["mode"].as<std::string>() == "test") {

    std::string prng = vm["prng"].as<std::string>();

    std::cout << "================================================\n";
    std::cout << "Starting non-uniform random bit generation tests\n";
    std::cout << "================================================\n";
    std::cout << "Chosen PRNG: " << prng << "\n";

    size_t testSize = vm["test-size"].as<size_t>();

    std::vector<double> chiSquareTest;
    std::vector<double> runsChiSquareTest;
    std::vector<double> runsFFTTest;

    for (const auto fraction : fractions) {
      std::cout << "> Testing for fraction " << fraction << std::endl;

      chiSquareTest.push_back(testSimpleChiSquare(prng, seeds, fraction, testSize));
      runsChiSquareTest.push_back(testRunsChiSquare(prng, seeds, fraction, testSize));
//      std::cout << testMatricesRanks(prng, seeds, fraction, 1000, 128) << std::endl;
      runsFFTTest.push_back(testRunsAutoCorrelation(prng, seeds, fraction, testSize * fraction));
    }

    auto publishResults = [&](auto& stream) {

      stream << "PRNG            , " << std::setw(15) << prng << std::endl;
      stream << "Test size       , " << std::setw(15) << testSize << std::endl;
      stream << "Tested fracions , ";
      for (auto fr : fractions)
        stream << std::setw(15) << fr << ", ";
      stream << std::endl;
      stream << "Chi square      , ";
      for (auto chi : chiSquareTest)
        stream << std::setw(15) << chi << ", ";
      stream << std::endl;
      stream << "Runs Histo Chi2 , ";
      for (auto chi : runsChiSquareTest)
        stream << std::setw(15) << chi << ", ";
      stream << std::endl;
      stream << "Runs FFT stddev , ";
      for (auto std : runsFFTTest)
        stream << std::setw(15) << std << ", ";
      stream << std::endl;
    };

    std::cout << "============================= TEST RESULTS =============================" << std::endl;
    publishResults(std::cout);

    std::ofstream file;
    file.open(prng);
    publishResults(file);
    file.close();
  } else if (vm["mode"].as<std::string>() == "speedtest") {

    const std::vector<std::string> prngs = {
      "PRNG_PCG", "PRNG_Hash1", "PRNG_Hash2", "PRNG_TRandom", "PRNG_TRandom1", "PRNG_TRandom2", "PRNG_TRandom3",
      "PRNG_TRandomMixMax", "PRNG_TRandomMixMax17", "PRNG_TRandomMT64", "PRNG_TRandomRanlux48", "PRNG_Dummy"
    };

    std::vector<double> results;
    for( const auto& prng : prngs) {
      size_t testSize = vm["test-size"].as<size_t>();
      if (prng == "PRNG_TRandomMixMax") {
        testSize /= 5000;
      } else if (prng == "PRNG_TRandomMixMax17") {
        testSize /= 100;
      } else if (prng == "PRNG_TRandom3" || prng == "PRNG_TRandomMT64") {
        testSize /= 10;
      } 
      double nsPerCall = testSpeed(prng, testSize);
      results.push_back(nsPerCall);
    }

    auto publishResults = [&](auto& stream) {
      for (size_t i = 0; i < prngs.size(); i++) {
        stream << std::setw(30) << prngs[i] << ", " << std::setw(15) << results[i] << std::endl;
      }
    };

    std::cout << "============================= TEST RESULTS =============================" << std::endl;
    publishResults(std::cout);

    std::ofstream file;
    file.open("speedtest");
    publishResults(file);
    file.close();
  };

  return 0;
}
