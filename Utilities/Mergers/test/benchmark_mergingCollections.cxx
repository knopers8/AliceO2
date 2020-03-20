// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include <benchmark/benchmark.h>

#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TTree.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TRandom.h>
#include <TRandomGen.h>

#include <boost/histogram.hpp>

namespace bh = boost::histogram;

#include <ctime>

#define BENCHMARK_RANGE_COLLECTIONS Arg(1)->Arg(1 << 2)->Arg(1 << 4)->Arg(1 << 6)->Arg(1 << 8)->Arg(1 << 10)->Arg(1 << 12)

// A simple test where an input is provided
// and the subsequent InputRecord is immediately requested.
static void BM_RelayMessageCreation(benchmark::State& state)
{
//  Monitoring metrics;

  for (auto _ : state) {
    //state.PauseTiming();

    // do stuff here

    //state.ResumeTiming();
  }
}

//BENCHMARK(BM_RelayMessageCreation);

static void BM_mergingCollectionsTH1F(benchmark::State& state)
{

  size_t collectionSize = state.range(0);
  size_t bins = 62500; // makes 250kB

  TCollection* collection = new TObjArray();
  collection->SetOwner(true);
  TF1* uni = new TF1("uni", "1", 0, 1000000);
  for (size_t i = 0; i < collectionSize; i++) {
    TH1F* h = new TH1F(("test" + std::to_string(i)).c_str(), "test", bins, 0, 1000000);
    h->FillRandom("uni", 50000);
    collection->Add(h);
  }

  TH1F* m = new TH1F("merged", "merged", bins, 0, 1000000);

  for (auto _ : state) {
    m->Merge(collection, "-NOCHECK");
  }

  delete collection;
  delete m;
  delete uni;
}
BENCHMARK(BM_mergingCollectionsTH1F)->BENCHMARK_RANGE_COLLECTIONS;

static void BM_mergingCollectionsTH2F(benchmark::State& state)
{

  size_t collectionSize = state.range(0);
  size_t bins = 250; // 250 bins * 250 bins * 4B makes 250kB

  TCollection* collection = new TObjArray();
  collection->SetOwner(true);
  TF2* uni = new TF2("uni", "1", 0, 1000000, 0, 1000000);
  for (size_t i = 0; i < collectionSize; i++) {
    TH2F* h = new TH2F(("test" + std::to_string(i)).c_str(), "test", bins, 0, 1000000, bins, 0, 1000000);
    h->FillRandom("uni", 50000);
    collection->Add(h);
  }

  TH2F* m = new TH2F("merged", "merged", bins, 0, 1000000, bins, 0, 1000000);

  for (auto _ : state) {
    m->Merge(collection, "-NOCHECK");
  }

  delete collection;
  delete m;
}
//BENCHMARK(BM_mergingCollectionsTH2F)->BENCHMARK_RANGE_COLLECTIONS;

static void BM_mergingCollectionsTH3F(benchmark::State& state)
{
  size_t collectionSize = state.range(0);
  size_t bins = 40; // 250 bins * 250 bins * 250 bins * 4B makes 256kB

  TCollection* collection = new TObjArray();
  collection->SetOwner(true);
  TF3* uni = new TF3("uni", "1", 0, 1000000, 0, 1000000, 0, 1000000);
  for (size_t i = 0; i < collectionSize; i++) {
    TH3F* h = new TH3F(("test" + std::to_string(i)).c_str(), "test",
                       bins, 0, 1000000,
                       bins, 0, 1000000,
                       bins, 0, 1000000);
    h->FillRandom("uni", 50000);
    collection->Add(h);
  }

  TH3F* m = new TH3F("merged", "merged", bins, 0, 1000000, bins, 0, 1000000, bins, 0, 1000000);

  for (auto _ : state) {
    m->Merge(collection, "-NOCHECK");
  }

  delete collection;
  delete m;
}
//BENCHMARK(BM_mergingCollectionsTH3F)->BENCHMARK_RANGE_COLLECTIONS;

static void BM_mergingCollectionsTHNSparse(benchmark::State& state)
{
  size_t collectionSize = state.range(0);

  const Double_t min = 0.0;
  const Double_t max = 1000000.0;
  const size_t dim = 10;
  const Int_t bins = 250; // 250 bins * 250 bins * 250 bins * 4B makes 256kB
  const Int_t binsDims[dim] = {bins, bins, bins, bins, bins, bins, bins, bins, bins, bins};
  const Double_t mins[dim] = { min, min, min, min, min, min, min, min, min, min};
  const Double_t maxs[dim] = { max, max, max, max, max, max, max, max, max, max};

  TRandomMT64 gen;
  gen.SetSeed(std::time(nullptr));
  Double_t randomArray[dim];

  for (auto _ : state) {

    state.PauseTiming();
    TCollection* collection = new TObjArray();
    collection->SetOwner(true);
    for (size_t i = 0; i < collectionSize; i++) {

      auto* h = new THnSparseF(("test" + std::to_string(i)).c_str(), "test", dim, binsDims, mins, maxs);
      for (size_t entry = 0; entry < 50000; entry++) {
        gen.RndmArray(dim, randomArray);
        for (double r : randomArray) {
          r *= max;
        }
        h->Fill(randomArray);
      }
      collection->Add(h);
    }
    auto* m = new THnSparseF("merged", "merged", dim, binsDims, mins, maxs);

    state.ResumeTiming();
    m->Merge(collection);

    state.PauseTiming();

    delete collection;
    delete m;
  }
}
//BENCHMARK(BM_mergingCollectionsTHNSparse)->BENCHMARK_RANGE_COLLECTIONS;

static void BM_mergingPODCollections(benchmark::State& state)
{
  size_t collectionSize = state.range(0);
  size_t bins = 62500; // makes 250kB

  std::vector<std::vector<float>*> collection;
  TF1* uni = new TF1("uni", "1", 0, 1000000);

  const size_t randoms = 50000;
  TRandomMT64 gen;
  gen.SetSeed(std::time(nullptr));
  Double_t randomArray[randoms];

  for (size_t i = 0; i < collectionSize; i++) {
    auto* v = new std::vector<float>(bins, 0);
    gen.RndmArray(randoms, randomArray);
    for (double r : randomArray) {
      size_t idx = r * bins;
      if (idx != bins) {
        (*v)[idx] += 1;
      }
    }
    collection.push_back(v);
  }

  auto* m = new std::vector<float>(bins, 0);

  auto merge = [&](size_t i) {
    auto* v = collection[i];
    for (size_t b = 0; b < bins; b++) {
      (*m)[b] += (*v)[b];
    }
  };

  for (auto _ : state) {
    for (size_t i = 0; i < collectionSize; i++) {
      merge(i);
    }
  }

  for (size_t i = 0; i < collectionSize; i++) {
    delete collection[i];
  }
  delete m;
  delete uni;
}
BENCHMARK(BM_mergingPODCollections)->BENCHMARK_RANGE_COLLECTIONS;

static void BM_mergingBoostCollections(benchmark::State& state)
{
  const double min = 0.0;
  const double max = 1000000.0;
  const size_t collectionSize = state.range(0);
  const size_t bins = 62500; // makes 250kB

  auto merged = bh::make_histogram(bh::axis::regular<>(bins, min, max, "x"));

  std::vector<decltype(merged)> collection;
  TF1* uni = new TF1("uni", "1", 0, 1000000);

  const size_t randoms = 50000;
  TRandomMT64 gen;
  gen.SetSeed(std::time(nullptr));
  Double_t randomArray[randoms];

  for (size_t i = 0; i < collectionSize; i++) {
    //    auto* v = new std::vector<float>(bins, 0);
    collection.emplace_back(std::move(bh::make_histogram(bh::axis::regular<>(bins, min, max, "x"))));

    auto& h = collection.back();
    static_assert(std::is_reference<decltype(h)>::value);

    gen.RndmArray(randoms, randomArray);
    for (double r : randomArray) {
      h(r * max);
    }
  }

  auto merge = [&](size_t i) {
    merged += collection[i];
  };

  for (auto _ : state) {
    for (size_t i = 0; i < collectionSize; i++) {
      merge(i);
    }
  }

  delete uni;
}
BENCHMARK(BM_mergingBoostCollections)->BENCHMARK_RANGE_COLLECTIONS;

BENCHMARK_MAIN();
