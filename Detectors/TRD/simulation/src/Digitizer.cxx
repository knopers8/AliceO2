// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <TGeoManager.h>
#include <TRandom.h>

#include "TRDSimulation/Digitizer.h"
#include "TRDBase/TRDGeometry.h"
#include "TRDBase/TRDPadPlane.h"

#include "FairLogger.h"

using namespace o2::trd;

Digitizer::Digitizer()
{
  // Check if you need more initialization
  mGeom = new TRDGeometry();

  // Get the Ionization energy
  if (TRDCommonParam::Instance()->IsXenon()) {
    mWion = 23.53; // Ionization energy XeCO2 (85/15)
  } else if (TRDCommonParam::Instance()->IsArgon()) {
    mWion = 27.21; // Ionization energy ArCO2 (82/18)
  } else {
    LOG(FATAL) << "Wrong gas mixture";
    // add hard exit here!
  }
}

Digitizer::~Digitizer()
{
}

void Digitizer::process(std::vector<o2::trd::HitType> const& hits, std::vector<o2::trd::Digit>& digits)
{
  // (WIP) Implementation for digitization

  // Check if Geometry and if CCDB are available as they will be requiered
  /*
   mGeom = new TRDGeometry();
   TRDCalibDB * calibration = new TRDCalibDB();
   const int nTimeBins = calibration->GetNumberOfTimeBinsDCS();
  */

  // Loop over all TRD detectors
  // Get the a hit container for all the hits in a given detector then call convertHits for a given detector (0 - 539)
  const int kNdet = 540; // Get this from TRD Geometry
  int totalNumberOfProcessedHits = 0;
  LOG(INFO) << "Start of processing " << hits.size() << " hits";
  for (int det = 0; det < kNdet; ++det) {
    // Loop over all TRD detectors

    // Jump to the next detector if the detector is
    // switched off, not installed, etc
    /*
    if (calibration->IsChamberNoData(det)) {
      continue;
    }
    if (!mGeo->ChamberInGeometry(det)) {
      continue
    }
    */

    mHitContainer.clear();
    // Skip detectors without hits
    if (!getHitContainer(det, hits, mHitContainer)) {
      // move to next det if no hits are found for this det
      continue;
    }
    totalNumberOfProcessedHits += mHitContainer.size();
    int signals = 0; // dummy variable for now
    if (!convertHits(det, mHitContainer, signals)) {
      LOG(INFO) << "TRD converstion of hits failed for detector " << det;
      signals = 0; //
    }

    digits.emplace_back();
  } // end of loop over detectors
  LOG(INFO) << "End of processing " << totalNumberOfProcessedHits << " hits";
}

bool Digitizer::getHitContainer(const int det, const std::vector<o2::trd::HitType>& hits, std::vector<o2::trd::HitType>& hitContainer)
{
  //
  // Fills the hit vector for hits in detector number det
  // Returns false if there are no hits in the dectector
  //
  for (const auto& hit : hits) {
    if (hit.GetDetectorID() == det) {
      hitContainer.push_back(hit);
    }
  }
  if (hitContainer.size() == 0) {
    return false;
  }
  return true;
}

bool Digitizer::convertHits(int det, const std::vector<o2::trd::HitType>& hits, int& arraySignal)
{
  //
  // Convert the detector-wise sorted hits to detector signals
  //
  LOG(INFO) << "Start converting " << hits.size() << " hits for detector " << det;

  // Dummy signal for now
  arraySignal = -1;
  // Number of pads included in the pad response
  constexpr int kNpad = 3;
  // Width of the amplification region
  const float kAmWidth = TRDGeometry::amThick();
  // Width of the drift retion
  const float kDrWidth = TRDGeometry::drThick();
  // Drift + Amplification region
  const float kDrMin = -0.5 * kAmWidth;
  const float kDrMax = kDrWidth + 0.5 * kAmWidth;

  int timeBinTRFend = 0;

  double pos[3];
  double loc[3];
  double padSignal[kNpad];
  double signalOld[kNpad];

  // Number of track dictionary arrays
  // const Int_t kNdict     = AliTRDdigitsManager::kNDict;
  // AliTRDarrayDictionary *dictionary[kNdict];

  TRDSimParam* simParam = TRDSimParam::Instance();
  TRDCommonParam* commonParam = TRDCommonParam::Instance();
  // TRDCalibDB *calibration = TRDCalibDB::Instace();
  if (!simParam) {
    LOG(FATAL) << "TRD Simulation Parameters not available";
    return false;
  }
  if (!commonParam) {
    LOG(FATAL) << "TRD Common Parameters not available";
    return false;
  }
  // if (!calibration) {
  //   LOG(FATAL) << "TRD Calibration database not available";
  //   return false;
  // }

  // Get the detector wise calibration objects
  // TRDCalROC* calVdriftROC = 0;
  float calVdriftDetValue = 0.0;
  // TRDCalDet* calVdriftDet = calibration->GetVdriftDet();
  // TRDCalROC* calT0ROC = 0;
  float calT0DetValue = 0.0;
  // TRDCalDet* calT0Det = calibration->GetT0Det();
  double calExBDetValue = 0.0;
  // TRDCalDet* calExBDet = calibration->GetExBDet();

  if (simParam->TRFOn()) {
    timeBinTRFend = ((int)(simParam->GetTRFhi() * commonParam->GetSamplingFrequency())) - 1;
  }

  const int nTimeTotal = 0; //DigitsManager->GetDigitsParam()->GetNTimeBins(det);
  const float samplingRate = commonParam->GetSamplingFrequency();
  const float elAttachProp = simParam->GetElAttachProp() / 100;

  const TRDPadPlane* padPlane = mGeom->getPadPlane(det); // Check if mGeom is working
  const int layer = mGeom->getLayer(det);
  const float row0 = padPlane->getRow0ROC();
  const int nRowMax = padPlane->getNrows();
  const int nColMax = padPlane->getNcols();

  // Allocate space for signals
  // arraySignal->Allocate(nRowMax,nColMax,nTimeTotal);

  // Create a new array for the dictionary
  // for (int dict = 0; dict < kNdict; dict++) {
  //   dictionary[dict] = (AliTRDarrayDictionary *) fDigitsManager->GetDictionary(det,dict);
  //   dictionary[dict]->Allocate(nRowMax,nColMax,nTimeTotal);
  // }

  // Loop over hits
  for (const auto& hit : hits) {
    pos[0] = hit.GetX();
    pos[1] = hit.GetY();
    pos[2] = hit.GetZ();

    const float eDep = hit.GetEnergyLoss();
    const int qTotal = (int)eDep / mWion;

    gGeoManager->SetCurrentPoint(pos);
    gGeoManager->FindNode();
    // int inDrift = 1;
    // if (strstr(gGeoManager->GetPath(), "/UK")) {
    //   inDrift = 0;
    // }
    // Better use ternary operation and make it const
    const int inDrift = strstr(gGeoManager->GetPath(), "/UK") ? 0 : 1; // If /UK, then not in drift region

    // Get the calibration objects
    // They are all zeros until I have implemented the calibration objects
    // calVdriftROC = 0;      // calibration->GetVdriftROC(det);
    // calVdriftDetValue = 0; // calVdriftDet->GetValue(det);
    // calT0ROC = 0;          // calibration->GetT0ROC(det);
    // calT0DetValue = 0;     // calT0Det->GetValue(det);
    // calExBDetValue = 0;    // calExBDet->GetValue(det);

    // Go to the local coordinate system
    // loc [0] -  col direction in amplification or drift volume
    // loc [1] -  row direction in amplification or drift volume
    // loc [2] -  time direction in amplification or drift volume
    gGeoManager->MasterToLocal(pos, loc);
    if (inDrift) {
      loc[2] = loc[2] - kDrWidth / 2 - kAmWidth / 2;
    }
    // The drift length in cm without diffusion yet!
    const double driftLength = -1 * loc[2];

    /*
      Remember that there is a patch for TR photons
      TR photons are not available yet
    */

    int rowE = padPlane->getPadRowNumberROC(loc[1]);
    if (rowE < 0) {
      continue;
    }

    double rowOffset = padPlane->getPadRowOffsetROC(rowE, loc[1]);
    double offsetTilt = padPlane->getTiltOffset(rowOffset);
    int colE = padPlane->getPadColNumber(loc[0] + offsetTilt);
    if (colE < 0) {
      continue;
    }

    // Normalized drift length
    // Commented out what is still not yet implemented
    double driftVelocity = calVdriftDetValue; // * calVdriftROC->GetValue(colE, rowE);
    double absDriftLength = abs(driftLength);
    if (commonParam->ExBOn()) {
      absDriftLength /= TMath::Sqrt(1 / (1 + calExBDetValue * calExBDetValue));
    }

    // Loop over all created electrons
    const int nElectrons = (int)abs(qTotal);
    for (int el = 0; el < nElectrons; ++el) {
      /* 
      Now the real local coordinate system of the ROC
      column direction: locC
      row direction:    locR 
      time direction:   locT
      locR and locC are identical to the coordinates of the corresponding
      volumina of the drift or amplification region.
      locT is defined relative to the wire plane (i.e. middle of amplification
      region), meaning locT = 0, and is negative for hits coming from the
      drift region. 
      */
      double locC = loc[0];
      double locR = loc[1];
      double locT = loc[2];

      // Electron attachment
      if (TRDSimParam::Instance()->ElAttachOn()) {
        if (gRandom->Rndm() < absDriftLength * elAttachProp) {
          continue;
        }
      }

      // Apply diffusion smearing
      if (simParam->DiffusionOn()) {
        if (!diffusion(driftVelocity, absDriftLength, calExBDetValue, locR, locC, locT)) {
          continue;
        }
      }

      // Apply E x B effects
      if (commonParam->ExBOn()) {
        locC = locC + calExBDetValue * driftLength;
      }

      // The electron position after diffusion and ExB in pad coordinates.
      rowE = padPlane->getPadRowNumberROC(locR);
      if (rowE) {
        continue;
      }
      rowOffset = padPlane->getPadRowOffsetROC(rowE, locR);

      // The
      offsetTilt = padPlane->getTiltOffset(rowOffset);
      colE = padPlane->getPadColNumber(locC + offsetTilt);
      if (colE < 0) {
        continue;
      }
      double colOffset = padPlane->getPadColOffset(colE, locC + offsetTilt);

      // Retrieve drift velocity becuase col and row may have changed
      driftVelocity = calVdriftDetValue; // * calVdriftROC->GetValue(colE, rowE);
      float t0 = calT0DetValue;          // + calT0ROC->getValue(colE, rowE);

      // Convert the position to drift time [mus], using either constant drift velocity or
      // time structure of drift cells (non-isochronity, GARFIELD calculation).
      // Also add absolute time of hits to take pile-up events into account properly
      double driftTime;
      if (simParam->TimeStructOn()) {
        // Get z-position with respect to anode wire
        double zz = row0 - locR + padPlane->getAnodeWireOffset();
        zz -= ((int)(2 * zz)) / 2;
        if (zz > 0.25) {
          zz = 0.5 - zz;
        }
        // Use drift time map (GARFIELD)
        driftTime = commonParam->TimeStruct(driftVelocity, 0.5 * kAmWidth - 1.0 * locT, zz) + hit.GetTime();
      } else {
        // Use constant drift velocity
        driftTime = abs(locT) / driftVelocity + hit.GetTime();
      }

      // Apply the gas gain including fluctuations
      double ggRndm = 0.0;
      do {
        ggRndm = gRandom->Rndm();
      } while (ggRndm <= 0);
      double signal = -(simParam->GetGasGain()) * TMath::Log(ggRndm);

      // Apply the pad response
      if (simParam->PRFOn()) {
        // The distance of the electron to the center of the pad
        // in units of pad width
        double dist = (colOffset - 0.5 * padPlane->getColSize(colE)) / padPlane->getColSize(colE);
        // ********************************************************************************
        // This is a fixed parametrization, i.e. not dependent on calibration values !
        // ********************************************************************************
        // ************ CHECK calibration, added it when you have it
        // if (!(calibration->PadResponse(signal, dist, layer, padSignal))) {
        //   continue;
        // }
      } else {
        padSignal[0] = 0.0;
        padSignal[1] = signal;
        padSignal[2] = 0.0;
      }

      // The time bin (always positive), with t0 distortion
      double timeBinIdeal = driftTime * samplingRate + t0;
      // Protection
      if (abs(timeBinIdeal) > 2 * nTimeTotal) {
        timeBinIdeal = 2 * nTimeTotal;
      }
      int timeBinTruncated = ((int)timeBinIdeal);
      // The distance of the position to the middle of the timebin
      double timeOffset = ((float)timeBinTruncated + 0.5 - timeBinIdeal) / samplingRate;

      // Sample the time response inside the drift region + additional time bins before and after.
      // The sampling is done always in the middle of the time bin
      const int firstTimeBin = TMath::Max(timeBinTruncated, 0);
      const int lastTimeBin = TMath::Min(timeBinTruncated + timeBinTRFend, nTimeTotal);
      for (int iTimeBin = firstTimeBin; iTimeBin < lastTimeBin; ++iTimeBin) {
        // Apply the time response
        double timeResponse = 1;
        double crossTalk = 0;
        double time = (iTimeBin - timeBinTruncated) / samplingRate + timeOffset;
        if (simParam->TRFOn()) {
          timeResponse = simParam->TimeResponse(time);
        }
        if (simParam->CTOn()) {
          crossTalk = simParam->CrossTalk(time);
        }
        signalOld[0] = 0;
        signalOld[1] = 0;
        signalOld[2] = 0;
        for (int iPad = 0; iPad < kNpad; iPad++) {
          int colPos = colE + iPad - 1;
          if (colPos < 0) {
            continue;
          }
          if (colPos >= nColMax) {
            break;
          }
          // Add the signals
          // signalOld[iPad] = arraySignal->GetData(rowE, colPos, iTimeBin); ******* TO BE ADDED LATER REQUIRES ARRAYSIGNALS
          if (colPos != colE) {
            // Cross talk added to non-central pads
            signalOld[iPad] += padSignal[iPad] * (timeResponse + crossTalk);
          } else {
            // Without cross talk at central pad
            signalOld[iPad] += padSignal[iPad] * timeResponse;
          }

          // ******* TO BE ADDED LATER REQUIRES ARRAYSIGNALS
          // arraySignal->SetData(rowE, colPos, iTimeBin, signalOld[iPad]); ******* TO BE ADDED LATER REQUIRES ARRAYSIGNALS
          // ******* TO BE ADDED LATER REQUIRES ARRAYSIGNALS

          // Store the track index in the dictionary
          // Note: We store index+1 in order to allow the array to be compressed
          // Note2: Taking out the +1 in track
          /*  
          if (signalOld[iPad] > 0) {
            for (int dict = 0; dict < kNdict; dict++) {
              int oldTrack = dictionary[dict]->GetData(rowE, colPos, iTimeBin);
              if (oldTrack == track) {
                break;
              }
              if (oldTrack == -1) {
                dictionary[dict]->SetData(rowE, colPos, iTimeBin, track);
                break;
              }
            }
          }
          */
        } // Loop: pads
      }   // Loop: time bins
    }     // end of loop over electrons
  }       // end of loop over hits
  return true;
}

bool Digitizer::diffusion(float vdrift, double absdriftlength, double exbvalue, double& lRow, double& lCol, double& lTime)
{
  //
  // Applies the diffusion smearing to the position of a single electron.
  // Depends on absolute drift length.
  //
  float diffL = 0.0;
  float diffT = 0.0;
  if (TRDCommonParam::Instance()->GetDiffCoeff(diffL, diffT, vdrift)) {
    float driftSqrt = TMath::Sqrt(absdriftlength);
    float sigmaT = driftSqrt * diffT;
    float sigmaL = driftSqrt * diffL;
    lRow = gRandom->Gaus(lRow, sigmaT);
    if (TRDCommonParam::Instance()->ExBOn()) {
      lCol = gRandom->Gaus(lCol, sigmaT * 1.0 / (1.0 + exbvalue * exbvalue));
      lTime = gRandom->Gaus(lTime, sigmaL * 1.0 / (1.0 + exbvalue * exbvalue));
    } else {
      lCol = gRandom->Gaus(lCol, sigmaT);
      lTime = gRandom->Gaus(lTime, sigmaL);
    }
    return true;
  } else {
    return false;
  }
}