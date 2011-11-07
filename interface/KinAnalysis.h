#ifndef kinanalysis_h
#define kinanalysis_h

#include <vector>
#include <algorithm>

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TLorentzVector.h"
#include "TRandom2.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TF1.h"
#include "TH1.h"

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/TopKinSolver.h"
#include "LIP/Top/interface/KinResultsHandler.h"

typedef std::vector<TLorentzVector> TLorentzVectorCollection;
typedef std::pair<TLorentzVector, float> KinCandidate_t;
typedef std::vector<KinCandidate_t> KinCandidateCollection_t;

TLorentzVectorCollection randomlyRotate(TLorentzVectorCollection &leptonsp4, TLorentzVectorCollection &jetsp4); 

class KinAnalysis
{
 public:

  KinAnalysis(TString& scheme, int maxTries=10000, int maxJetMult=2, float mw=80.398, float mb=4.8, TString outpath="KinAnalysis.root", bool doWrite=true);

  static bool sortKinCandidates(KinCandidate_t a, KinCandidate_t b)   {   return (a.second>b.second);  }
  void runOn(top::EventSummary_t &ev, JetResolution *ptResol, JetResolution *etaResol, JetResolution *phiResol, JetCorrectionUncertainty *jecUnc);
  void endAnalysis() { resHandler_.end(); }
  ~KinAnalysis();

 private:
  TString scheme_;
  int maxTries_,maxJetMult_;
  float mw_, mb_;
  TopKinSolver kin_;
  
  TRandom2 rndGen_;
  TF1 *deltaPzFunc_;
  
  KinResultsHandler resHandler_;
};

#endif
