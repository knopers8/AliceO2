// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#define BOOST_TEST_MODULE Test Framework DataSamplingCondition
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <iostream>
#include <vector>
#include <boost/test/unit_test.hpp>

#include "Framework/DataSamplingConditionFactory.h"
#include "Framework/DataRef.h"
#include "Framework/DataProcessingHeader.h"
#include "Headers/DataHeader.h"

using namespace o2::framework;
using namespace o2::header;
/*
BOOST_AUTO_TEST_CASE(DataSamplingConditionRandom)
{
  auto conditionRandom = DataSamplingConditionFactory::create("random");
  BOOST_REQUIRE(conditionRandom);

  // PRNG should behave the same every time and on every machine.
  // Of course, the test does not cover full range of timesliceIDs.
  std::vector<bool> correctDecision{
    false, true, false, true, true, false, true, false, false, true, false, true, false, false, false, false, false,
    true, false, false, true, true, false, false, true, true, false, false, false, false, true, true, false, false,
    true, true, false, false, false, false, false, true, false, false, false, false, false, true, false
  };
  boost::property_tree::ptree config;
  config.put("fraction", 0.5);
  config.put("seed", 943753948);
  conditionRandom->configure(config);

  for (DataProcessingHeader::StartTime id = 1; id < 50; id++) {
    DataProcessingHeader dph{ id, 0 };
    o2::header::Stack headerStack{ dph };
    DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };
    BOOST_CHECK_EQUAL(correctDecision[id - 1], conditionRandom->decide(dr));
  }
}

BOOST_AUTO_TEST_CASE(DataSamplingConditionPayloadSize)
{
  auto conditionPayloadSize = DataSamplingConditionFactory::create("payloadSize");
  BOOST_REQUIRE(conditionPayloadSize);

  boost::property_tree::ptree config;
  config.put("upperLimit", 500);
  config.put("lowerLimit", 30);
  conditionPayloadSize->configure(config);

  std::vector<std::pair<size_t, bool>> testCases{
    { 0, false },
    { 29, false },
    { 30, true },
    { 200, true },
    { 500, true },
    { 501, false }
  };

  for (const auto& t : testCases) {
    DataHeader dh;
    dh.payloadSize = t.first;
    o2::header::Stack headerStack{ dh };
    DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };
    BOOST_CHECK_EQUAL(conditionPayloadSize->decide(dr), t.second);
  }
}

BOOST_AUTO_TEST_CASE(DataSamplingConditionNConsecutive)
{
  auto conditionNConsecutive = DataSamplingConditionFactory::create("nConsecutive");
  BOOST_REQUIRE(conditionNConsecutive);

  boost::property_tree::ptree config;
  config.put("samplesNumber", 3);
  config.put("cycleSize", 10);
  conditionNConsecutive->configure(config);

  std::vector<std::pair<size_t, bool>> testCases{
    { 0, true },
    { 1, true },
    { 2, true },
    { 3, false },
    { 8, false },
    { 9, false },
    { 9999999999999, false },
    { 10000000000000, true },
    { 10000000000001, true },
    { 10000000000002, true },
    { 10000000000003, false }
  };

  for (const auto& t : testCases) {
    DataProcessingHeader dph{ t.first, 0 };
    o2::header::Stack headerStack{ dph };
    DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };
    BOOST_CHECK_EQUAL(conditionNConsecutive->decide(dr), t.second);
  }
}

 */


#include <algorithm>
#include <chrono>
using namespace std::chrono;

BOOST_AUTO_TEST_CASE(DataSamplingConditionPCG)
{
  auto conditionPCG = DataSamplingConditionFactory::create("pcg");
  BOOST_REQUIRE(conditionPCG);

  boost::property_tree::ptree config;
  config.put("fraction", 0.5);
  config.put("seed", 938475231*10);
  conditionPCG->configure(config);

  uint64_t total = 0;
  size_t vsize = 10000;

  std::vector<DataProcessingHeader::StartTime> v(vsize);
  for(DataProcessingHeader::StartTime id = 0; id < vsize; id++) {
    v[id] = id;
  }
  std::random_shuffle(v.begin(), v.end(), [](int i) { return std::rand()%i;});

  steady_clock::time_point timeStart = steady_clock::now();
  for (auto id : v) {
//    LOG(INFO) << id;
    DataProcessingHeader dph{ id, 0 };
    o2::header::Stack headerStack{ dph };
    DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };

    total += conditionPCG->decide(dr);
  }

  std::cout << "pcg total: " << total << ", "
            << duration_cast<nanoseconds>(steady_clock::now() - timeStart).count() / vsize << "ns/call" << std::endl;

}


#include <algorithm>
#include <chrono>
using namespace std::chrono;

BOOST_AUTO_TEST_CASE(DataSamplingConditionHash)
{
  auto conditionRandom = DataSamplingConditionFactory::create("hash");
  BOOST_REQUIRE(conditionRandom);

  // PRNG should behave the same every time and on every machine.
  // Of course, the test does not cover full range of timesliceIDs.
//  std::vector<bool> correctDecision{
//    false, true, false, true, true, false, true, false, false, true, false, true, false, false, false, false, false,
//      true, false, false, true, true, false, false, true, true, false, false, false, false, true, true, false, false,
//      true, true, false, false, false, false, false, true, false, false, false, false, false, true, false
//  };
  boost::property_tree::ptree config;
  config.put("fraction", 0.001);
  config.put("seed", 938475231);
  conditionRandom->configure(config);

  steady_clock::time_point timeStart = steady_clock::now();
  size_t vsize = 10000000;
  int b = 0;
  for (DataProcessingHeader::StartTime id = 1; id < vsize; id++) {
    DataProcessingHeader dph{ id, 0 };
    o2::header::Stack headerStack{ dph };
    DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };
//    BOOST_CHECK_EQUAL(correctDecision[id - 1], conditionRandom->decide(dr));
//    LOG(INFO) << conditionRandom->decide(dr);
    b += conditionRandom->decide(dr);
  }
  std::cout << "hash total: " << b << ", "
            << duration_cast<nanoseconds>(steady_clock::now() - timeStart).count() / vsize << "ns/call" << std::endl;
}


BOOST_AUTO_TEST_CASE(DataSamplingConditionHashCombine)
{
  auto conditionRandom = DataSamplingConditionFactory::create("hashCombine");
  BOOST_REQUIRE(conditionRandom);

  // PRNG should behave the same every time and on every machine.
  // Of course, the test does not cover full range of timesliceIDs.
//  std::vector<bool> correctDecision{
//    false, true, false, true, true, false, true, false, false, true, false, true, false, false, false, false, false,
//      true, false, false, true, true, false, false, true, true, false, false, false, false, true, true, false, false,
//      true, true, false, false, false, false, false, true, false, false, false, false, false, true, false
//  };
  boost::property_tree::ptree config;
  config.put("fraction", 0.0001);
  config.put("seed", 123123123);
  conditionRandom->configure(config);

  steady_clock::time_point timeStart = steady_clock::now();
  size_t vsize = 10000000;
  int b = 0;
  for (DataProcessingHeader::StartTime id = 1; id < vsize; id++) {
    DataProcessingHeader dph{ id, 0 };
    o2::header::Stack headerStack{ dph };
    DataRef dr{ nullptr, reinterpret_cast<const char*>(headerStack.data()), nullptr };
//    BOOST_CHECK_EQUAL(correctDecision[id - 1], conditionRandom->decide(dr));
//    LOG(INFO) << conditionRandom->decide(dr);
    b += conditionRandom->decide(dr);
  }
  std::cout << "hashCombine total: " << b << ", "
            << duration_cast<nanoseconds>(steady_clock::now() - timeStart).count() / vsize << "ns/call" << std::endl;
}