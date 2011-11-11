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
#define MAXCATEGORIES 8
enum JetMultBins{BIN_2=0,BIN_3, BIN_4};
struct CombinedHFCModel_t
{
  RooArgSet pdfSet,constrPDFSet;
  RooRealVar *bmult, *r, *lfacceptance;
  RooRealVar *abseb,*sfeb,*sfeb_mean_constrain,*sfeb_sigma_constrain,*diff_sfeb[MAXCATEGORIES];
  RooFormulaVar *eb,*diff_eb[MAXCATEGORIES];
  RooGaussian *sfeb_constrain;
  RooRealVar *abseq,*sfeq,*sfeq_mean_constrain,*sfeq_sigma_constrain;
  RooFormulaVar *eq;
  RooGaussian *sfeq_constrain;
  RooRealVar *jetocc[MAXCATEGORIES];
  RooRealVar *fcorrect[MAXCATEGORIES],*fcorrect_mean_constrain[MAXCATEGORIES],*fcorrect_sigma_constrain[MAXCATEGORIES];
  RooGaussian *fcorrect_constrain[MAXCATEGORIES];
  RooRealVar *fttbar[MAXCATEGORIES],*fttbar_mean_constrain[MAXCATEGORIES],*fttbar_sigma_constrain[MAXCATEGORIES];
  RooGaussian *fttbar_constrain[MAXCATEGORIES];
  RooRealVar *fsingletop[MAXCATEGORIES],*fsingletop_mean_constrain[MAXCATEGORIES],*fsingletop_sigma_constrain[MAXCATEGORIES];
  RooGaussian *fsingletop_constrain[MAXCATEGORIES];
};

//
class HFCMeasurement
{
 public:

  enum FitTypes { FIT_R, FIT_EB, FIT_R_AND_EB, FIT_R_AND_XSEC, FIT_EB_AND_XSEC, FIT_EB_AND_EQ };
  
  /**
     @short CTOR
   */
  HFCMeasurement(int fitType=0, int maxJets=3) : 
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
    void fitHFCtoEnsemble(top::EventSummaryHandler &evHandler);
    
    /**
       @short setters for parameters
    */
    void setStandardModelR(float r=1.0) { smR_=r; }

    void configureBtagAlgo(TString btagAlgo,double cut)
    {
      btagAlgo_ = btagAlgo;
      algoCut_  = cut;
    }
    
    void setBtagEfficiency(double eff, double sfactor, double sfactorUnc)
    {
      effb_   = eff;
      sfb_    = sfactor;
      sfbUnc_ = sfactorUnc;
    } 

    void setMistagEfficiency(double eff, double sfactor, double sfactorUnc)
    {
      effq_   = eff;
      sfq_    = sfactor;
      sfqUnc_ = sfactorUnc;
    } 

    void setSelectionFractions(double fcorrect,   double fcorrectunc, 
			       double fttbar,     double fttbarunc,
			       double fsingletop, double fsingletopunc,
			       int jetBin=0,      TString dilChannel="emu")
    {
      dilChannel += jetBin;
      fcorrect_[dilChannel]      = fcorrect;
      fcorrectUnc_[dilChannel]   = fcorrectunc;
      fttbar_[dilChannel]        = fttbar;
      fttbarUnc_[dilChannel]     = fttbarunc;
      fsingletop_[dilChannel]    = fsingletop;
      fsingletopUnc_[dilChannel] = fsingletopunc;
    }

  
    /**
       @short save results
    */
    void saveMonitoringHistograms(TString tag);
    
    CombinedHFCModel_t model;

 private:

    void initHFCModel();
    void runHFCFit();
    void runHFCDiffFit(TString dilCat);

    void bookMonitoringHistograms();
    void resetHistograms();
    void resetModelValues();

    bool isInit_;
    int fitType_;

    int maxJets_;
    std::vector<TString> categoryKeys_;

    //R
    double smR_;
    
    //btag algorithm 
    TString btagAlgo_;    
    double algoCut_;
    Float_t  effb_, sfb_, sfbUnc_, effq_, sfq_, sfqUnc_;

    //parametrizations for purity 
    std::map<TString, Float_t> fcorrect_,  fcorrectUnc_, 
      fttbar_, fttbarUnc_,
      fsingletop_, fsingletopUnc_;
    
    //the histograms
    SelectionMonitor controlHistos_;    

    //a counter for pseudo-experiments
    int nMeasurements_;
};


#endif
