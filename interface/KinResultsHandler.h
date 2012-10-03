#ifndef kinresultshandler_hh
#define kinresultshandler_hh

#include <iostream>

#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuSummaryHandler.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"

//
TChain *getKinChainFrom(TString outpath);

class KinResultsHandler
{
 public:
  
  KinResultsHandler();
  void init(TString outpath,bool doWrite,int maxJetMult=2);
  void addResults(int run , int lumi, int event);
  inline TTree *getResultsTree() { return kinTree_; }
  inline TChain *getResultsChain() { return kinChain_; }
  void end();
  std::vector<double> getMPVEstimate(TH1 *);
  inline void getEventInfo(Int_t &run, Int_t &event, Int_t &lumi)
    {
      run=iRun_;
      event=iEvent_;
      lumi=iLumi_;
    }
  void resetHistos();
  TH1F *getHisto(TString var, int nComb);
  inline std::map<std::pair<TString, int>,TH1F *> &getControlHistos() { return kinHistos_; }
  ~KinResultsHandler();
  
 private:
  
  //
  bool doWrite_;

  //aux variables
  TFile *kinFile_;
  TTree *kinTree_;
  TChain *kinChain_;
  TF1 *fitFunc_;
  Int_t iRun_, iEvent_, iLumi_;
  
  std::map<std::pair<TString, int>,TH1F *> kinHistos_;
};


#endif
