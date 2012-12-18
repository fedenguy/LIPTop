#ifndef hfcmeasurement_hh
#define hfcmeasurement_hh

#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"
#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuSummaryHandler.h"

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
#include "RooDataSet.h"
#include "RooAddition.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooLognormal.h"
#include "RooUniform.h"

#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ModelConfig.h"



//
#define MAXJETMULT 12
#define MAXCATEGORIES 12
enum JetMultBins{BIN_2=0,BIN_3, BIN_4};
struct CombinedHFCModel_t
{
  RooAbsPdf *pdf;
  RooProdPdf *constrPdf;
  RooArgSet pdfConstrains;
  std::map<TString, RooAbsPdf *> pdfForCategory;
  RooCategory *sample;
  RooRealVar *r;
  RooRealVar *bmult, *bmultObs, *eventYields;
  RooRealVar *acc1,*acc05,*acc0;
  RooFormulaVar *alpha[MAXCATEGORIES], *alpha2[MAXCATEGORIES], *alpha1[MAXCATEGORIES], *alpha0[MAXCATEGORIES];
  RooRealVar *abseb,*sfeb,*sfebExp;
  RooFormulaVar *eb,*ebcorr[MAXCATEGORIES];
  RooGaussian *sfeb_constrain;
  RooRealVar *abseq,*sfeq,*sfeqExp;
  RooFormulaVar *eq,*eqcorr[MAXCATEGORIES];
  RooGaussian *sfeq_constrain;
  RooRealVar *jetocc[MAXCATEGORIES];
  RooRealVar *fcorrect[MAXCATEGORIES],*fcorrectExp[MAXCATEGORIES];
  RooAbsPdf *fcorrect_constrain[MAXCATEGORIES];
  RooRealVar *fttbar[MAXCATEGORIES], *fttbarExp[MAXCATEGORIES];
  RooAbsPdf *fttbar_constrain[MAXCATEGORIES];
  RooRealVar *fsingletop[MAXCATEGORIES], *fsingletopExp[MAXCATEGORIES];
  RooAbsPdf *fsingletop_constrain[MAXCATEGORIES];
  
  RooWorkspace *ws;
  RooStats::ModelConfig *modelConfig;
  
  Double_t rFitLowerLimit,         rFitUpperLimit;
  Double_t rFitResult,           rFitResultAsymmErrHi,           rFitResultAsymmErrLo;
  Double_t rFit[MAXCATEGORIES],  rFitAsymmErrHi[MAXCATEGORIES],  rFitAsymmErrLo[MAXCATEGORIES];
  Double_t rCombFit[MAXCATEGORIES],  rCombFitAsymmErrHi[MAXCATEGORIES],  rCombFitAsymmErrLo[MAXCATEGORIES];
  Double_t ebFit[MAXCATEGORIES], ebFitAsymmErrHi[MAXCATEGORIES], ebFitAsymmErrLo[MAXCATEGORIES];
  Double_t eqFit[MAXCATEGORIES], eqFitAsymmErrHi[MAXCATEGORIES], eqFitAsymmErrLo[MAXCATEGORIES];
  Double_t minLL;
};

//
class HFCMeasurement
{
 public:

  enum FitTypes        {FIT_R, FIT_EB, FIT_R_AND_EB, FIT_R_AND_XSEC, FIT_EB_AND_XSEC, FIT_EB_AND_EQ};
  
  
  /**
     @short CTOR
   */
  HFCMeasurement(int fitType, TString fitConfig, TString wpConfig="");

  /**
     @short parses configuration file and adds/updates values to the model
  */
  void parseFitConfig(TString url);

  /**
     @short init the PDFs for the fit
  */
  void initHFCModel();
    
  /**
     @short DTOR
  */
  ~HFCMeasurement() { }
  
 private:

  RooWorkspace *ws_;
  RooStats::ModelConfig *mc_;
  std::set<string> sampleCats_;
  
 public:
  

  enum EventCategories {EE_2JETS, EE_3JETS, MUMU_2JETS, MUMU_3JETS, EMU_2JETS, EMU_3JETS};
  enum NuisanceTypes   {GAUSSIAN, UNIFORM, LOGNORMAL };
  enum RunMode         {FITINCLUSIVEONLY, FITEXCLUSIVECATEGORIES };
  SelectionMonitor controlHistos_;    
  CombinedHFCModel_t model;

    //check me this point fwd        
    /**
       @short steer the fit
    */
    void fitHFCtoEnsemble(ZZ2l2nuSummaryHandler &evHandler, int runMode, bool debug=false);
    void fitHFCtoMeasurement(std::vector<TH1D *> &btagHistos,  int runMode, bool debug=false);
    
    /**
       @short setters for parameters
    */
    void setStandardModelR(float r=1.0) { smR_=r; }

    /**
       @short configuration of the b-tag algorithm
    */
    void configureBtagAlgo(TString btagAlgo,double cut)
    {
      btagAlgo_ = btagAlgo;
      algoCut_  = cut;
    }
    
    /**
       @short set MC expected efficiency and data-driven measurement of the scale-factor
     */
    void setBtagEfficiency(double eff, double sfactor, double sfactorUnc)
    {
      effb_   = eff;
      sfb_    = sfactor;
      sfbUnc_ = sfactorUnc;
    } 

    /**
       @short set MC expected efficiency and data-driven measurement of the scale-factor
    */
    void setMistagEfficiency(double eff, double sfactor, double sfactorUnc)
    {
      effq_   = eff;
      sfq_    = sfactor;
      sfqUnc_ = sfactorUnc;
    } 

    /**
       @short event modeling per category: used also to instantiate the categories in the fit
     */
    void setParametersForCategory(double fcorrect,   double fcorrectunc, 
				  double fttbar,     double fttbarunc,
				  double fsingletop, double fsingletopunc,
				  double ebcorr,     double eqcorr,
				  int jetBin=0,      TString dilChannel="emu")
    {
      TString tag(dilChannel); tag+=jetBin;      
      fcorrect_[tag]      = fcorrect;
      fcorrectUnc_[tag]   = fcorrectunc;
      fttbar_[tag]        = fttbar;
      fttbarUnc_[tag]     = fttbarunc;
      fsingletop_[tag]    = fsingletop;
      fsingletopUnc_[tag] = fsingletopunc;
      btagEffCorr_[tag]   = ebcorr;
      ltagEffCorr_[tag]   = eqcorr;
      categoryKeys_.insert(tag);
    }

    /**
       @short dump fitter configuration
     */
    void printConfiguration(std::ostream &os);


    /**
       @short resets the current model values to the default ones
     */
    void resetModelValues();

    /**
       @short steers the fit 
     */
    void runHFCFit(int runMode,bool debug);


    //internal parameters
    int fitType_, nuisanceType_;
    TString fitTypeTitle_, fitTypeName_;
    int maxJets_;
    double smR_;
    TString btagAlgo_;    
    double algoCut_;
    Float_t  effb_, sfb_, sfbUnc_, effq_, sfq_, sfqUnc_;
    std::set<TString> categoryKeys_;
    std::map<TString, Float_t> fcorrect_,  fcorrectUnc_, fttbar_, fttbarUnc_,  fsingletop_, fsingletopUnc_, btagEffCorr_, ltagEffCorr_;
    int nMeasurements_;

    RooStats::FeldmanCousins *fc_;

 private :


};


#endif
