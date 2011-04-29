#ifndef kinresultshandler_hh
#define kinresultshandler_hh

#include <iostream>

#include "LIP/Top/interface/EventSummaryHandler.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1D.h"

class KinResultsHandler
{
 public:
  
  KinResultsHandler(TString outpath,bool doWrite);
  void addResults(EventSummary_t &ev);
  inline TTree *getResultsTree() { return kinTree_; }
  void end();
  std::vector<double> getMPVEstimate(TH1 *);
  inline void getEventInfo(Int_t &run, Int_t &event, Int_t &lumi)
    {
      run=iRun_;
      event=iEvent_;
      lumi=iLumi_;
    }
  void bookHistos(int maxJetMult=2);
  void resetHistos();
  TH1D *getHisto(TString var, int nComb);
  inline std::map<std::pair<TString, int>,TH1D *> &getControlHistos() { return kinHistos_; }
  ~KinResultsHandler();
  
 private:
  
  void bookTree();
  void appendTree();

  //
  bool doWrite_;

  //aux variables
  TFile *kinFile_;
  TTree *kinTree_;
  TF1 *fitFunc_;
  Int_t iRun_, iEvent_, iLumi_;
  
  std::map<std::pair<TString, int>,TH1D *> kinHistos_;
};


#endif
