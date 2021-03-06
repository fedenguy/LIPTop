#ifndef kinanalysis_h
#define kinanalysis_h

#include <vector>
#include <algorithm>

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TRandom2.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TF1.h"
#include "TH1.h"

#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuPhysicsEvent.h"
//#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/TopKinSolver.h"
#include "LIP/Top/interface/KinResultsHandler.h"

typedef std::vector<LorentzVector> LorentzVectorCollection;
typedef std::pair<LorentzVector, float> KinCandidate_t;
typedef std::vector<KinCandidate_t> KinCandidateCollection_t;

//
class KinAnalysis
{
 public:

  KinAnalysis(TString& scheme, int maxTries=10000, int maxJetMult=2, float mw=80.398, float mb=4.8, TString outpath="KinAnalysis.root", bool doWrite=true);

  static bool sortKinCandidates(KinCandidate_t a, KinCandidate_t b)   {   return (a.second>b.second);  }
  void runOn(PhysicsEvent_t &ev, JetResolution *ptResol, JetResolution *etaResol, JetResolution *phiResol, JetCorrectionUncertainty *jecUnc, bool isMC);
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
