#ifndef hfcmeasurement_hh
#define hfcmeasurement_hh

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/HeavyFlavorPDF.h"

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

#endif

//
#define MAXJETMULT 5
enum JetMultBins{BIN_2=0,BIN_3,BIN_4,BIN_5};
struct CombinedHFCModel_t
{
  RooArgSet pdfSet,constrPDFSet;
  RooRealVar *bmult, *r, *lfacceptance;
  RooRealVar *eb,*eb_mean_constrain,*eb_sigma_constrain;
  RooGaussian *eb_constrain;
  RooRealVar *eq,*eq_mean_constrain,*eq_sigma_constrain;
  RooGaussian *eq_constrain;
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

  HFCMeasurement(int maxJets=4,TString wp="loose",int fitType=0) : maxJets_(maxJets), wp_(wp), fitType_(fitType)
    {
      bookMonitoringHistograms(maxJets_);
      
      //btagcut["loose"]=1.7;  
      btagcut["loose"]=2.0;  
      btageff["loose"]=0.788;
      btageffunc["loose"]=0.006;
      //      mistagrate["loose"]=0.182;
      //      mistagrateunc["loose"]=0.005;
      mistagrate["loose"]=0.137;
      mistagrateunc["loose"]=0.005;

      btagcut["medium"]=3.3;  
      btageff["medium"]=0.644;
      btageffunc["medium"]=0.006;
      mistagrate["medium"]=0.067;
      mistagrateunc["medium"]=0.003;

      btagcut["tight"]=3.41; 
      btageff["tight"]=0.433;
      btageffunc["tight"]=0.004;
      mistagrate["tight"]=0.031;
      mistagrateunc["tight"]=0.002;

      initHFCModel(maxJets_,wp_,fitType_);
    }


  ~HFCMeasurement()
    {
    }

  void fitHFCtoEnsemble(TTree *t, EventSummary_t *evt, bool debug=true);
  void fitHFCto(std::vector<TH1D *> &bmultH, bool debug=true);

  void bookMonitoringHistograms(int maxJets);
  void showMonitoringHistograms(bool debug=true);

  int maxJets_;
  TString wp_;
  int fitType_;
  TRandom2 rndGen;
  TH1D * biasMonH, *statMonH, *pullMonH;
  CombinedHFCModel_t model;
  std::vector<TH1D *> bmultH;
  std::map<TString,Float_t> btagcut, btageff, btageffunc, mistagrate, mistagrateunc;

 private:

  void initHFCModel(int maxJets=4,TString wp="loose",int fitType=0);

};


#endif
