#include <iostream>
#include <random>

#include <TRandom.h>
#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>

//todo: test if shuffling i changes results

int rootEachTime(long times){

  TRandom gen;
//  TRandom::Binomial(1, 0.01);
  int res = 0;

  for(long i = 0; i < times; ++i){
    gen.SetSeed(i);
    res += gen.Binomial(1, 0.01);
  }

  return res;
}

int rootEachTime1(long times){

  TRandom1 gen;
//  TRandom::Binomial(1, 0.01);
  int res = 0;

  for(long i = 0; i < times; ++i){
    gen.SetSeed(i);
    res += gen.Binomial(1, 0.01);
  }

  return res;
}

int rootEachTime2(long times){

  TRandom2 gen;
//  TRandom::Binomial(1, 0.01);
  int res = 0;

  for(long i = 0; i < times; ++i){
    gen.SetSeed(i);
    res += gen.Binomial(1, 0.01);
  }

  return res;
}

int rootEachTime3(long times){

  TRandom3 gen;
//  TRandom::Binomial(1, 0.01);
  int res = 0;

  for(long i = 0; i < times; ++i){
    gen.SetSeed(i);
    res += gen.Binomial(1, 0.01);
  }

  return res;
}

int initEachTime(long times){

  int res = 0;
  std::bernoulli_distribution d(0.01);
  for(long i = 0; i < times; ++i){
    std::minstd_rand gen(i);
    res += d(gen);
  }

  return res;
}

int initOnce(long times) {
  int res = 0;
  std::bernoulli_distribution d(0.01);
  std::minstd_rand gen(123);
  for(long i = 0; i < times; ++i) {

//    std::mt19937 gen(i);
    res += d(gen);
  }

  return res;
}

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
#include <chrono>
//#include <boost/program_options/variables_map.hpp>
//#include <boost/program_options/options_description.hpp>
#include <boost/program_options.hpp>


using namespace std::chrono;
namespace bpo = boost::program_options;


int main(int argc, char *argv[]) {

//  long times = 100000000;
//  std::cout << "initEachTime" << std::endl;
//  std::cout << initEachTime(times) << std::endl;
//
//  std::cout << "initOnce" << std::endl;
//  std::cout << initOnce(times) << std::endl;
//
//  std::cout << "rootEachTime" << std::endl;
//  std::cout << rootEachTime(times) << std::endl;
//
//  std::cout << "rootEachTime1" << std::endl;
//  std::cout << rootEachTime1(times) << std::endl;
//
//  std::cout << "rootEachTime2" << std::endl;
//  std::cout << rootEachTime2(times) << std::endl;
//
//  std::cout << "rootEachTime3" << std::endl;
//  std::cout << rootEachTime3(times) << std::endl;

//  auto conditionPCG = DataSamplingConditionFactory::create("pcg");
//
//  boost::property_tree::ptree config;
//  config.put("fraction", 0.5);
//  config.put("seed", 9437537*10);
//  conditionPCG->configure(config);
//
//  uint64_t total = 0;
//  size_t vsize = 10000000;
//
//  std::vector<DataProcessingHeader::StartTime> v(vsize);
//  for(DataProcessingHeader::StartTime id = 0; id < vsize; id++) {
//    v[id] = id;
//  }
////  std::random_shuffle(v.begin(), v.end(), [](int i) { return std::rand()%i;});
//  steady_clock::time_point timeStart = steady_clock::now();
//
//  for (auto id : v) {
////    LOG(INFO) << id;
//    DataProcessingHeader dph{ id, 0 };
//    o2::header::Stack headerStack{ dph };
//    DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };
//
//    total += conditionPCG->decide(dr);
//  }
//
//  std::cout << "total: " << total << ", "
//            << duration_cast<nanoseconds>(steady_clock::now() - timeStart).count() / vsize << "ns/call" << std::endl;
//

  // Arguments parsing
  bpo::variables_map vm;
  bpo::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Produce help message.")
    ("prng,p",      bpo::value<std::string>()->default_value("pcg"),    R"(prng version, "pcg", "hash", "hashCombine".)");
  bpo::store(parse_command_line(argc, argv, desc), vm);
  bpo::notify(vm);

  std::string prng = vm["prng"].as<std::string>();
  auto cond = DataSamplingConditionFactory::create(prng);
  uint64_t i = 0;
  while (cond) {
    if (prng == "pcg") {
      uint32_t rnd = cond->rnd(i);
      write(1, &rnd, sizeof(uint32_t));
    } else if (prng == "hash") {
      uint64_t rnd = cond->rnd(i);
      write(1, &rnd, sizeof(uint64_t));
    } else if (prng == "hashCombine") {
      uint64_t rnd = cond->rnd(i);
      write(1, &rnd, sizeof(uint64_t));
    }
    i++;
  }

  return 0;
}