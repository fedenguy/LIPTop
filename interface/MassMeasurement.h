
#ifndef _mass_measurement_hh_
#define _mass_measurement_hh_

#if !defined(__CINT__) || defined(__MAKECINT__)

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooLandau.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooPolyVar.h"
#include "Roo1DTable.h"
#include "RooSimPdfBuilder.h"
#include "RooSimultaneous.h"
#include "RooDataHist.h"

#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2DErrors.h"
#include "TMultiGraph.h"
#include "TTree.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TFile.h"
#include "TList.h"
#include "TIterator.h"
#include "TObject.h"
#include "TVirtualFitter.h"
#include "TMatrixD.h"
#include "TObjArray.h"
#include "TList.h"
#include "TKey.h"
#include "TDirectory.h"

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"


#include <set>

#endif

#define MAXEVENTS 10000

// a measurement of the mass for a set of events
struct EnsembleMeasurement_t
{
  Bool_t status;
  Int_t nEvents;
  Float_t mass, err, errLo, errHigh, evMasses[MAXEVENTS], evSample[MAXEVENTS];
};

//a summary of the measurements from the unbinned likelihood fits
struct MassFitResults_t
{
  Bool_t status;
  Float_t npoints,ll;
  Float_t tMass, tMassErrHigh, tMassErrLo, tMassErr;
  Float_t nSig,  nSigErrHigh,  nSigErrLo,  nSigErr;
  Float_t nBckg, nBckgErrHigh, nBckgErrLo, nBckgErr;
  Float_t jes, lltmass, lltmassLo, lltmassHigh;
};


class MassMeasurement
{
 public:

  MassMeasurement(TString parfileURL)
    {
      fitPars_ = ParseParametersFrom(parfileURL);
    }
  
  /**
     @short performs the standard unbinned likelihood fit to a set of mass measurements
  */
  EnsembleMeasurement_t DoMassFit(top::EventSummaryHandler &evHandler, bool debug=false);

  /**
     @short fits the mass to an ensemble
  */
  MassFitResults_t DoMassFit(EnsembleMeasurement_t &em, bool debug=false);

 private:
      
  /**
     @short the Roofit interface to the template fit
  */
  MassFitResults_t CombinedMassFitter(TTree *data, bool debug=false);
  
  /**
     @short parses a file with fit parameters
  */
  std::map<TString, Double_t> fitPars_;
  std::map<TString,Double_t> ParseParametersFrom(TString parfileURL);
};


#endif
