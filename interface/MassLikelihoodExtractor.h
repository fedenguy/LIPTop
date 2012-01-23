#ifndef masslikelihoodextractor_h
#define masslikelihoodextractor_h

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "TSystem.h"
#include "TObject.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// Includes and forward declarations
#include "TChain.h"
#include "LIP/Top/interface/KinResultsHandler.h"
class TH1F;

#endif

//Pietro
//an exercise you can do:
//
//- create a new mtop histogram called kinideogram (same number of bins
//						  and range as the histograms you just sent)
//
//- loop over all events
// . for each event get the histograms for combination #1 and #2
// . count the entries entries of combination #1 and add to each bin
//the following baseline
//
// (Max-Min)/(10000-N)
//
//   where Max (Min) is the maximum (minimum) top mass in the x range
//      N is the number of entries in the histogram
//     10000 is the total number of attempts for KIN to find solutions
// . repeat the same for histogram #2
// . add both histograms (the result is a likelihood of the event for mtop)
// . add to each bin in the total kinideoagram histogram the logarithm
//of the corresponding bin of the previous histogram
// (the final result should be a likelihood for mtop using all events)
//
//Can you save the output histogram of each event and of the kinideogram
//so we can take a look at it?
// Thanks,
//Pedro

class MassLikelihoodExtractor
{
 public:
  
  /**
     @short CTOR
  */
  MassLikelihoodExtractor(bool);
  ~MassLikelihoodExtractor();
  
  void processEvent(KinResultsHandler&, TChain&, Int_t);
  TH1F getEventLikelihood(/* event parameter */);
  TH1F getFinalLikelihood();

 private:
  // data members
  bool debug_;
  TH1F* totalMassLikelihood_;

  TH1F* hcomb1_;
  TH1F* hcomb2_;
  TH1F* evMassLikelihood_;

  // methods
  void applyBaseline(TH1F&, int);
  Double_t getBaselineFromHisto(TH1F&);
  void fillFinalLikelihood();
  


};

#endif

