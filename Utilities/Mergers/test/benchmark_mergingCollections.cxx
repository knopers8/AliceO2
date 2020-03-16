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
#include <THn.h>
#include <TTree.h>
#include <THnSparse.h>

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

BENCHMARK(BM_RelayMessageCreation);

static void BM_mergingCollectionsTH1F(benchmark::State& state) {

  size_t collectionSize = state.range(0);
  size_t bins = 1000;

  TCollection* collection = new TObjArray();
  collection->SetOwner(true);
  for (size_t i = 0; i < collectionSize; i++) {
    TH1F* h = new TH1F(("test" + std::to_string(i)).c_str(), "test", bins, 0, 1000000);
    h->FillRandom("gaus", 5000); // with that enabled, the benchmark cannot exit.
    collection->Add(h);
  }

  TH1F* m = new TH1F();

  for (auto _ : state) {
    m->Merge(collection);
  }

  delete collection;
  delete m;
}
BENCHMARK(BM_mergingCollectionsTH1F)->Arg(1)->Arg(1 << 2)->Arg(1 << 4)->Arg(1 << 6)->Arg(1 << 8)->Arg(1 << 10)->Arg(1 << 12);

BENCHMARK_MAIN();
