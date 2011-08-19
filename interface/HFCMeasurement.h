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

  enum EventCategoriesForMeasurement { ALLDileptons=500, SFDileptons, OFDileptons };
  enum FitTypes { FIT_R, FIT_EB, FIT_R_AND_EB, FIT_R_AND_XSEC, FIT_EB_AND_XSEC, FIT_EB_AND_EQ };
  
  /**
     @short CTOR
   */
  HFCMeasurement(int maxJets=4,TString btagAlgo="TCHEL", int eventCategory=OFDileptons, int fitType=0) : maxJets_(maxJets), btagAlgo_(btagAlgo), eventCategory_(eventCategory), fitType_(fitType), smR_(1.0), nMeasurements_(0)  
    {
      bookMonitoringHistograms();
      
      algoCut["TCHEL"]=1.7; 
      effb["TCHEL"]=0.78;  sfb["TCHEL"]=0.95;  sfbUnc["TCHEL"]=sqrt(pow(0.01,2)+pow(0.1,2)); 
      effq["TCHEL"]=0.1;   sfq["TCHEL"]=1.11;  sfqUnc["TCHEL"]=sqrt(pow(0.01,2)+pow(0.12,2));

      algoCut["TCHEM"]=3.3;   
      effb["TCHEM"]=0.78;  sfb["TCHEM"]=0.94;  sfbUnc["TCHEM"]=sqrt(pow(0.01,2)+pow(0.09,2)); 
      effq["TCHEM"]=0.1;   sfq["TCHEM"]=1.21;  sfqUnc["TCHEM"]=sqrt(pow(0.02,2)+pow(0.17,2));

      algoCut["TCHPT"]=3.41; 
      effb["TCHPT"]=0.78;  sfb["TCHPT"]=0.88;  sfbUnc["TCHPT"]=sqrt(pow(0.02,2)+pow(0.09,2)); 
      effq["TCHPT"]=0.1;   sfq["TCHPT"]=1.21;  sfqUnc["TCHPT"]=sqrt(pow(0.10,2)+pow(0.18,2));

      algoCut["JBPL"]=1.33; 
      effb["JBPL"]=0.78;  sfb["JBPL"]=sfb["TCHEL"]; sfbUnc["JBPL"]=sfbUnc["TCHEL"];
      effq["JBPL"]=0.1;   sfq["JBPL"]=sfq["TCHEL"]; sfqUnc["JBPL"]=sfqUnc["TCHEL"];

      algoCut["JBPM"]=2.55; 
      effb["JBPM"]=0.78;  sfb["JBPM"]=sfb["TCHEM"]; sfbUnc["JBPM"]=sfbUnc["TCHEM"];
      effq["JBPM"]=0.1;   sfq["JBPM"]=sfq["TCHEM"]; sfqUnc["JBPM"]=sfqUnc["TCHEM"];

      algoCut["JBPT"]=3.74; 
      effb["JBPT"]=0.78;  sfb["JBPT"]=sfb["TCHET"]; sfbUnc["JBPT"]=sfbUnc["TCHET"];
      effq["JBPT"]=0.1;   sfq["JBPT"]=sfq["TCHET"]; sfqUnc["JBPT"]=sfqUnc["TCHET"];

      algoCut["SSVHEM"]=1.74;
      effb["SSVHEM"]=0.78;  sfb["SSVHEM"]=0.95;  sfbUnc["SSVHEM"]=sqrt(pow(0.01,2)+pow(0.1,2)); 
      effq["SSVHEM"]=0.1;   sfq["SSVHEM"]=0.91;  sfqUnc["SSVHEM"]=sqrt(pow(0.02,2)+pow(0.10,2));

      alpha2[0] = 0.63;      alpha2Unc[0]=sqrt(pow(0.01,2)+pow((0.03+0.01)*0.5,2));
      alpha2[2] =alpha0[0];  alpha2Unc[2]=alpha2Unc[0];
      alpha2[3] =alpha0[2];  alpha2Unc[3]=alpha2Unc[2];

      alpha0[0]=0.135;       alpha0[0]=sqrt(pow(0.007,2)+pow((0.006+0.003)/2,2));
      alpha0[2] =alpha0[0];  alpha0Unc[2]=alpha0Unc[0];
      alpha0[3] =alpha0[2];  alpha0Unc[3]=alpha0Unc[2];
    }

    /**
       @short DTOR
    */
    ~HFCMeasurement() { }
    
    /**
       @short steer the fit
    */
    void fitHFCtoEnsemble(EventSummaryHandler &evHandler, TString btagAlgo);
    
    /**
       @short setters for parameters
    */
    void setStandardModelR(float r=1.0) { smR_=r; }
    
    /**
       @short save results
    */
    void saveMonitoringHistograms(TString tag);

    CombinedHFCModel_t model;

 private:

    void initHFCModel();
    void bookMonitoringHistograms();
    void runHFCFit();
    void resetHistograms();

    int maxJets_;
    TString btagAlgo_;
    int eventCategory_;
    int fitType_;
    
    SelectionMonitor controlHistos_;
    
    //ugly containers for values
    std::map<TString,Float_t> algoCut, effb, sfb, sfbUnc, effq, sfq, sfqUnc;
    std::map<Int_t, Float_t> alpha2,alpha2Unc, alpha0, alpha0Unc;
    
    double smR_;
    
    int nMeasurements_;
};


#endif
