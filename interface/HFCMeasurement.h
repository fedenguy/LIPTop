#ifndef hfcmeasurement_hh
#define hfcmeasurement_hh

#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/HeavyFlavorPDF.h"

#include "TFile.h"
#include "TRandom2.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooDataHist.h"
#include "RooAddition.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"

//
#define MAXJETMULT 8
enum JetMultBins{BIN_2=0,BIN_3,BIN_4,BIN_5, BIN_6, BIN_7, BIN_8};
struct CombinedHFCModel_t
{
  RooArgSet pdfSet,constrPDFSet;
  RooRealVar *bmult, *r, *lfacceptance;
  RooRealVar *abseb,*sfeb,*sfeb_mean_constrain,*sfeb_sigma_constrain;
  RooFormulaVar *eb;
  RooGaussian *sfeb_constrain;
  RooRealVar *abseq,*sfeq,*sfeq_mean_constrain,*sfeq_sigma_constrain;
  RooFormulaVar *eq;
  RooGaussian *sfeq_constrain;
  RooRealVar *alpha2,*alpha2_mean_constrain,*alpha2_sigma_constrain;
  RooGaussian *alpha2_constrain;
  RooFormulaVar *alpha1;
  RooRealVar *alpha0,*alpha0_mean_constrain,*alpha0_sigma_constrain;
  RooGaussian *alpha0_constrain;
};

//
class HFCMeasurement
{
 public:

  enum FitTypes { FIT_R, FIT_EB, FIT_R_AND_EB, FIT_R_AND_XSEC, FIT_EB_AND_XSEC, FIT_EB_AND_EQ };
  
  /**
     @short CTOR
   */
  HFCMeasurement(int maxJets=4, int fitType=0) : 
    isInit_(false), fitType_(fitType), maxJets_(maxJets),  smR_(1.0), nMeasurements_(0)  
    {
      bookMonitoringHistograms();
    }

    /**
       @short DTOR
    */
    ~HFCMeasurement() { }
    
    
    /**
       @short steer the fit
    */
    void fitHFCtoEnsemble(EventSummaryHandler &evHandler, TString dilCat);
    
    /**
       @short setters for parameters
    */
    void setStandardModelR(float r=1.0) { smR_=r; }

    void configureBtagAlgo(TString btagAlgo,double cut)
    {
      btagAlgo_ = btagAlgo;
      algoCut_  = cut;
    }
    
    void setBtagEfficiency(double eff, double sfactor, double sfactorUnc,int jetBin=0)
    {
      effb_[jetBin]=eff;
      sfb_[jetBin]=sfactor;
      sfbUnc_[jetBin]=sfactorUnc;
    } 

    void setMistagEfficiency(double eff, double sfactor, double sfactorUnc,int jetBin=0)
    {
      effq_[jetBin]=eff;
      sfq_[jetBin]=sfactor;
      sfqUnc_[jetBin]=sfactorUnc;
    } 

    void setAlpha(double a2, double a2unc, double a0, double a0unc, int jetBin=0)
    {
      alpha2_[jetBin]=a2;
      alpha2Unc_[jetBin]=a2unc;
      alpha0_[jetBin]=a0;
      alpha0Unc_[jetBin]=a0unc;
    }

  
    /**
       @short save results
    */
    void saveMonitoringHistograms(TString tag);
    
    CombinedHFCModel_t model;

 private:

    void initHFCModel();
    void runHFCFit(TString dilCat);

    void bookMonitoringHistograms();
    void resetHistograms();
    void resetModelValues();

    bool isInit_;
    int fitType_;

    int maxJets_;

    //R
    double smR_;
    
    //btag algorithm 
    TString btagAlgo_;    
    double algoCut_;
    std::map<Int_t,Float_t>  effb_, sfb_, sfbUnc_, effq_, sfq_, sfqUnc_;

    //event types
    std::map<Int_t, Float_t> alpha2_,  alpha2Unc_, alpha0_, alpha0Unc_;
    
    SelectionMonitor controlHistos_;    
    int nMeasurements_;
};


#endif
