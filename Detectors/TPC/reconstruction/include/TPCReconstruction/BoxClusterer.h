// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file BoxClusterer.h
/// \brief Class for TPC cluster finding
#ifndef ALICEO2_TPC_BoxClusterer_H_
#define ALICEO2_TPC_BoxClusterer_H_

#include <vector>
#include <memory>

#include "Rtypes.h"
#include "TPCReconstruction/Clusterer.h"
#include "DataFormatsTPC/Cluster.h"

#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2{

namespace TPC
{

class ClusterContainer;

class BoxClusterer : public Clusterer
{

  using MCLabelContainer = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;

 public:
  /// Constructor
  /// \param output is pointer to vector to be filled with clusters
  /// \param rowsMax Max number of rows to process
  /// \param padsMax Max number of pads to process
  /// \param timeBinsMax Max number of timebins to process
  /// \param minQMax Minimum peak charge for cluster
  /// \param requirePositiveCharge Positive charge is required
  /// \param requireNeighbouringPad Requires at least 2 adjecent pads with charge above threshold
  BoxClusterer(std::vector<o2::TPC::Cluster>* output,
               int rowsMax = 18,
               int padsMax = 138,
               int timeBinsMax = 1024,
               int minQMax = 5,
               bool requirePositiveCharge = true,
               bool requireNeighbouringPad = true);

  /// Destructor
  ~BoxClusterer() override;

  /// Steer conversion of points to digits
  /// \param digits Container with TPC digits
  /// \param mcDigitTruth MC Digit Truth container
  void process(std::vector<o2::TPC::Digit> const& digits, MCLabelContainer const* mcDigitTruth) override;
  void finishProcess(std::vector<o2::TPC::Digit> const& digits, MCLabelContainer const* mcDigitTruth) override{};

  void setRowsMax(int val) { mRowsMax = val; };
  void setPadsMax(int val) { mPadsMax = val; };
  void setTimeBinsMax(int val) { mTimeBinsMax = val; };
  void setMinQMax(float val) { mMinQMax = val; };
  void setRequirePositiveCharge(bool val) { mRequirePositiveCharge = val; };
  void setRequireNeighbouringPad(bool val) { mRequireNeighbouringPad = val; };

  int getRowsMax() const { return mRowsMax; };
  int getPadsMax() const { return mPadsMax; };
  int getTimeBinsMax() const { return mTimeBinsMax; };
  float getMinQMax() const { return mMinQMax; };
  bool hasRequirePositiveCharge() const { return mRequirePositiveCharge; };
  bool hasRequireNeighbouringPad() const { return mRequireNeighbouringPad; };

 private:
  // To be done
  /* BoxClusterer(const BoxClusterer &); */
  /* BoxClusterer &operator=(const BoxClusterer &); */

  void findLocalMaxima(const Int_t iCRU);
  void cleanArrays();
  void getPadAndTimeBin(Int_t bin, Short_t& iPad, Short_t& iTimeBin);
  Int_t update(const Int_t iCRU, const Int_t iRow, const Int_t iPad,
               const Int_t iTimeBin, Float_t signal);
  Float_t getQ(const Float_t* adcArray, const Short_t pad,
               const Short_t time, Short_t& timeMin, Short_t& timeMax,
               Short_t& padMin, Short_t& padMax) const;
  Bool_t updateCluster(Float_t charge, Int_t deltaPad, Int_t deltaTime,
                       Float_t& qTotal, Double_t& meanPad,
                       Double_t& sigmaPad, Double_t& meanTime,
                       Double_t& sigmaTime);

  int mRowsMax;                 ///< Maximum row number
  int mPadsMax;                 ///< Maximum pad number
  int mTimeBinsMax;             ///< Maximum time bin
  float mMinQMax;               ///< Minimun Qmax for cluster
  bool mRequirePositiveCharge;  ///< If true, require charge > 0
  bool mRequireNeighbouringPad; ///< If true, require 2+ pads minimum

  //
  //  Expand buffer
  //
  Float_t** mAllBins;  //!<! Array for digit using random access
  Int_t** mAllSigBins; //!<! Array of pointers to the indexes over threshold
  Int_t* mAllNSigBins; //!<! Array with number of signals in each row

  std::vector<o2::TPC::Cluster>* mClusterArray; ///< Internal cluster storage
};
}
}


#endif
