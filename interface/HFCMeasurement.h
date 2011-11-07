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
#define MAXCATEGORIES 5
enum JetMultBins{BIN_2=0,BIN_3,BIN_4,BIN_5, BIN_6, BIN_7, BIN_8};
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
  RooRealVar *jetocc[MAXJETMULT];
  RooRealVar *fcorrect[MAXJETMULT],*fcorrect_mean_constrain[MAXJETMULT],*fcorrect_sigma_constrain[MAXJETMULT];
  RooGaussian *fcorrect_constrain[MAXJETMULT];
  RooRealVar *fttbar[MAXJETMULT],*fttbar_mean_constrain[MAXJETMULT],*fttbar_sigma_constrain[MAXJETMULT];
  RooGaussian *fttbar_constrain[MAXJETMULT];
  RooRealVar *fsingletop[MAXJETMULT],*fsingletop_mean_constrain[MAXJETMULT],*fsingletop_sigma_constrain[MAXJETMULT];
  RooGaussian *fsingletop_constrain[MAXJETMULT];
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
    void fitHFCtoEnsemble(top::EventSummaryHandler &evHandler, TString dilCat);
    
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

    void setSelectionFractions(double fcorrect,   double fcorrectunc, 
			       double fttbar,     double fttbarunc,
			       double fsingletop, double fsingletopunc,
			       int jetBin=0)
    {
      fcorrect_[jetBin]=fcorrect;
      fcorrectUnc_[jetBin]=fcorrectunc;
      fttbar_[jetBin]=fttbar;
      fttbarUnc_[jetBin]=fttbarunc;
      fsingletop_[jetBin]=fsingletop;
      fsingletopUnc_[jetBin]=fsingletopunc;
    }

  
    /**
       @short save results
    */
    void saveMonitoringHistograms(TString tag);
    
    CombinedHFCModel_t model;

 private:

    void initHFCModel();
    void runHFCFit(TString dilCat);
    void runHFCDiffFit(TString dilCat);

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
    std::map<Int_t, Float_t> fcorrect_,  fcorrectUnc_, 
      fttbar_, fttbarUnc_,
      fsingletop_, fsingletopUnc_;
    
    SelectionMonitor controlHistos_;    
    int nMeasurements_;
};


#endif
