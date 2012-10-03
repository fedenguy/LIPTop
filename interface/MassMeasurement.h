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
#include "RooAddition.h"
#include "RooNLLVar.h"

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
#include "TPaveText.h"

#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuSummaryHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"

#include <string.h>
#include <set>

#endif

#define MAXEVENTS 50000
#define MAXFITCATEGORIES 12

// a measurement of the mass for a set of events
struct EnsembleMeasurement_t
{
  Bool_t status;
  Int_t nEvents;
  Float_t mass, err, evMasses[MAXEVENTS], evCategories[MAXEVENTS];
};

//a summary of the measurements from the unbinned likelihood fits
struct MassFitResults_t
{
  Bool_t status;
  Float_t npoints,ll;
  Float_t tMass, tMassErrHigh, tMassErrLo, tMassErr;
  Float_t iTmass[MAXFITCATEGORIES],iTmassErr[MAXFITCATEGORIES];
  Float_t nSig,  nSigErrHigh,  nSigErrLo,  nSigErr;
  Float_t nBckg, nBckgErrHigh, nBckgErrLo, nBckgErr;
  Float_t jes, lltmass, lltmassLo, lltmassHigh;
};


class MassMeasurement
{
 public:

  TString tag_;  
  MassMeasurement(TString parfileURL,TString tag);

  enum CategoryMode {INCLUSIVE=0,EXCLUSIVE, DILEPTONEXCLUSIVE};
  enum InclusiveCategories {EQ1BTAGS=0,GEQ2BTAGS};
  enum ExclusiveCategories {SF_EQ1BTAGS=0,SF_GEQ2BTAGS,OF_EQ1BTAGS,OF_GEQ2BTAGS};
  enum ExclusiveDileptonCategories {EE_EQ1BTAGS=0,EE_GEQ2BTAGS,MUMU_EQ1BTAGS,MUMU_GEQ2BTAGS,EMU_EQ1BTAGS,EMU_GEQ2BTAGS};
  inline int getCategorizationMode() { return (int) fitPars_["cattype"]; }
  inline int getNumberOfCategories() { return (int) fitPars_["ncategs"]; }
  
  /**
     @short performs the standard unbinned likelihood fit to a set of mass measurements
  */
  EnsembleMeasurement_t DoMassFit(ZZ2l2nuSummaryHandler &evHandler, bool debug=false);
  
  /**
     @short fits the mass to an ensemble
  */
  MassFitResults_t DoMassFit(EnsembleMeasurement_t &em, bool debug=false);

 private:

  /**
     @short start the model
   */
  void InitModel();
  RooDataSet *inclusiveData;
  RooRealVar *recoMass,*category,*topMass;
  std::vector<RooDataSet * > allData;
  std::vector<RooAbsPdf * > allPdfs;
  RooArgSet bckgpdfset,signalpdfset;
  std::vector<TPaveText * > allCaptions;
  std::vector<RooNLLVar *> allLL;
  RooArgSet constrParams,sigYieldParams;


  /**
     @short the Roofit interface to the template fit
  */
  MassFitResults_t CombinedMassFitter(bool debug=false);
  
  /**
     @short parses a file with fit parameters
  */
  std::map<TString, Double_t> fitPars_;
  std::map<int,TString> catTitles_;
  std::map<TString,Double_t> ParseParametersFrom(TString parfileURL);
};


#endif
