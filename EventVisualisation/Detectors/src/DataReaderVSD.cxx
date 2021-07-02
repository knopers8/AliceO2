// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataReaderVSD.cxx
/// \brief VSD specific reading from file(s) (Visualisation Summary Data)
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch

#include <EventVisualisationDetectors/DataReaderVSD.h>
#include <TSystem.h>
#include <TEveManager.h>
#include <TFile.h>
#include <TPRegexp.h>
#include <TEveTrackPropagator.h>
#include <TEveEventManager.h>
#include <TKey.h>

namespace o2
{
namespace event_visualisation
{

DataReaderVSD::DataReaderVSD(DataInterpreter* interpreter)
  : DataReader(interpreter),
    mFile(nullptr),
    mMaxEv(-1),
    mCurEv(-1)
{
}

DataReaderVSD::~DataReaderVSD()
{
  if (mEvDirKeys.size() > 0) {
    for (auto obj : mEvDirKeys) {
      delete obj;
    }
    mEvDirKeys.clear();
  }

  if (mFile) {
    mFile->Close();
    delete mFile;
    mFile = nullptr;
  }
}

void DataReaderVSD::open()
{
  TString ESDFileName = "events_0.root";
  Warning("GotoEvent", "OPEN");
  mMaxEv = -1;
  mCurEv = -1;
  mFile = TFile::Open(ESDFileName);
  if (!mFile) {
    Error("VSD_Reader", "Can not open file '%s' ... terminating.",
          ESDFileName.Data());
    gSystem->Exit(1);
  }

  assert(mEvDirKeys.size() == 0);

  TPMERegexp name_re("Event\\d+");
  TObjLink* lnk = mFile->GetListOfKeys()->FirstLink();
  while (lnk) {
    if (name_re.Match(lnk->GetObject()->GetName())) {
      mEvDirKeys.push_back((TKey*)lnk->GetObject());
    }
    lnk = lnk->Next();
  }

  mMaxEv = mEvDirKeys.size();
  if (mMaxEv == 0) {
    Error("VSD_Reader", "No events to show ... terminating.");
    gSystem->Exit(1);
  }
}

TObject* DataReaderVSD::getEventData(int ev)
{
  if (ev < 0 || ev >= this->mMaxEv) {
    Warning("GotoEvent", "Invalid event id %d.", ev);
    return nullptr;
  }
  this->mCurEv = ev;
  return this->mEvDirKeys[this->mCurEv]->ReadObj();
}

} // namespace event_visualisation
} // namespace o2
