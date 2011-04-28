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
  void end();
  std::vector<double> getMPVEstimate(TH1 *);
  void bookHistos(int maxJetMult=2);
  void resetHistos();
  TH1D *getHisto(TString var, int nComb);
  inline std::map<std::pair<TString, int>,TH1D *> &getControlHistos() { return kinHistos_; }
  ~KinResultsHandler();
  
 private:
  
  void bookTree(EventSummary_t &ev);

  //aux variables
  TFile *kinFile_;
  TTree *kinTree_;
  TF1 *fitFunc_;
  
  std::map<std::pair<TString, int>,TH1D *> kinHistos_;
};


#endif
