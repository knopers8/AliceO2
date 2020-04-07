// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   ClusterReaderSpec.cxx

#include <vector>

#include "TTree.h"

#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "ITSMFTWorkflow/ClusterReaderSpec.h"

using namespace o2::framework;
using namespace o2::itsmft;

namespace o2
{
namespace itsmft
{

ClusterReader::ClusterReader(o2::detectors::DetID id, bool useMC, bool useClFull, bool useClComp, bool usePatterns)
{
  assert(id == o2::detectors::DetID::ITS || id == o2::detectors::DetID::MFT);
  mDetNameLC = mDetName = id.getName();
  mUseMC = useMC;
  mUseClFull = useClFull;
  mUseClComp = useClComp;
  mUsePatterns = usePatterns;
  std::transform(mDetNameLC.begin(), mDetNameLC.end(), mDetNameLC.begin(), ::tolower);
}

void ClusterReader::init(InitContext& ic)
{
  mInputFileName = ic.options().get<std::string>((mDetNameLC + "-cluster-infile").c_str());
}

void ClusterReader::run(ProcessingContext& pc)
{

  if (mFinished) {
    return;
  }

  read();

  LOG(INFO) << mDetName << "ClusterReader pushes " << mClusROFRec.size() << " ROFRecords,"
            << mClusterArray.size() << " full clusters, " << mClusterCompArray.size()
            << " compact clusters";

  // This is a very ugly way of providing DataDescription, which anyway does not need to contain detector name.
  // To be fixed once the names-definition class is ready
  pc.outputs().snapshot(Output{mOrigin, mOrigin == o2::header::gDataOriginITS ? "ITSClusterROF" : "MFTClusterROF",
                               0, Lifetime::Timeframe},
                        mClusROFRec);
  if (mUseClFull) {
    pc.outputs().snapshot(Output{mOrigin, "CLUSTERS", 0, Lifetime::Timeframe}, mClusterArray);
  }
  if (mUseClComp) {
    pc.outputs().snapshot(Output{mOrigin, "COMPCLUSTERS", 0, Lifetime::Timeframe}, mClusterCompArray);
  }
  if (mUsePatterns) {
    pc.outputs().snapshot(Output{mOrigin, "PATTERNS", 0, Lifetime::Timeframe}, mPatternsArray);
  }
  if (mUseMC) {
    pc.outputs().snapshot(Output{mOrigin, "CLUSTERSMCTR", 0, Lifetime::Timeframe}, mClusterMCTruth);
  }

  mFinished = true;
  pc.services().get<ControlService>().endOfStream();
  pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
}

void ClusterReader::read()
{
  // load data from files
  TFile clFile(mInputFileName.c_str(), "read");
  if (clFile.IsZombie()) {
    LOG(FATAL) << "Failed to open cluster file " << mInputFileName;
  }
  TTree* clTree = (TTree*)clFile.Get(mClusTreeName.c_str());
  if (!clTree) {
    LOG(FATAL) << "Failed to load clusters tree " << mClusTreeName << " from " << mInputFileName;
  }
  if (!clTree->GetBranch((mDetName + mClusROFBranchName).c_str())) {
    LOG(FATAL) << "Failed to load clusters ROFrecords branch "
               << " from " << mInputFileName;
  }

  clTree->SetBranchAddress((mDetName + mClusROFBranchName).c_str(), &mClusROFRecPtr);

  if (mUseClFull) {
    clTree->SetBranchAddress((mDetName + mClusterBranchName).c_str(), &mClusterArrayPtr);
  }
  if (mUseClComp) {
    clTree->SetBranchAddress((mDetName + mClusterCompBranchName).c_str(), &mClusterCompArrayPtr);
  }
  if (mUsePatterns) {
    clTree->SetBranchAddress((mDetName + mClusterPattBranchName).c_str(), &mPatternsArrayPtr);
  }
  if (mUseMC) {
    if (clTree->GetBranch((mDetName + mClustMCTruthBranchName).c_str())) {
      clTree->SetBranchAddress((mDetName + mClustMCTruthBranchName).c_str(), &mClusterMCTruthPtr);
      LOG(INFO) << "Will use MC-truth from " << mDetName + mClustMCTruthBranchName;
    } else {
      LOG(INFO) << "MC-truth is missing";
      mUseMC = false;
    }
  }
  clTree->GetEntry(0);
}

DataProcessorSpec getITSClusterReaderSpec(bool useMC, bool useClFull, bool useClComp, bool usePatterns)
{
  std::vector<OutputSpec> outputSpec;
  outputSpec.emplace_back("ITS", "ITSClusterROF", 0, Lifetime::Timeframe);
  if (useClFull) {
    outputSpec.emplace_back("ITS", "CLUSTERS", 0, Lifetime::Timeframe);
  }
  if (useClComp) {
    outputSpec.emplace_back("ITS", "COMPCLUSTERS", 0, Lifetime::Timeframe);
  }
  if (usePatterns) {
    outputSpec.emplace_back("ITS", "PATTERNS", 0, Lifetime::Timeframe);
  }
  if (useMC) {
    outputSpec.emplace_back("ITS", "CLUSTERSMCTR", 0, Lifetime::Timeframe);
  }

  return DataProcessorSpec{
    "its-cluster-reader",
    Inputs{},
    outputSpec,
    AlgorithmSpec{adaptFromTask<ITSClusterReader>(useMC, useClFull, useClComp)},
    Options{
      {"its-cluster-infile", VariantType::String, "o2clus_its.root", {"Name of the input cluster file"}}}};
}

DataProcessorSpec getMFTClusterReaderSpec(bool useMC, bool useClFull, bool useClComp, bool usePatterns)
{
  std::vector<OutputSpec> outputSpec;
  outputSpec.emplace_back("MFT", "MFTClusterROF", 0, Lifetime::Timeframe);
  if (useClFull) {
    outputSpec.emplace_back("MFT", "CLUSTERS", 0, Lifetime::Timeframe);
  }
  if (useClComp) {
    outputSpec.emplace_back("MFT", "COMPCLUSTERS", 0, Lifetime::Timeframe);
  }
  if (usePatterns) {
    outputSpec.emplace_back("MFT", "PATTERNS", 0, Lifetime::Timeframe);
  }
  if (useMC) {
    outputSpec.emplace_back("MFT", "CLUSTERSMCTR", 0, Lifetime::Timeframe);
  }

  return DataProcessorSpec{
    "mft-cluster-reader",
    Inputs{},
    outputSpec,
    AlgorithmSpec{adaptFromTask<MFTClusterReader>(useMC, useClFull, useClComp)},
    Options{
      {"mft-cluster-infile", VariantType::String, "o2clus_mft.root", {"Name of the input cluster file"}}}};
}

} // namespace itsmft
} // namespace o2
