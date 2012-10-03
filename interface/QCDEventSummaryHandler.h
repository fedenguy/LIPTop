#ifndef qcdeventsummaryhandler_h
#define qcdeventsummaryhandler_h

#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>
#include "Math/LorentzVector.h"
#include "TTree.h"
#include "DataFormats/Math/interface/deltaR.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;

#define MAXQCDJETCANDIDATES 50
#define MAXQCDTRACKSINJETS 200

namespace qcd
{

  struct EventSummary_t
  {
    int Run, nPU, nPUtrue, nPV, nJet, nMuon, nBFromGSplit, nTrkInc;
    int BitTrigger, Jet_flavour[MAXQCDJETCANDIDATES], Jet_nFirstTrkInc[MAXQCDJETCANDIDATES], Jet_nLastTrkInc[MAXQCDJETCANDIDATES], Muon_IdxJet[MAXQCDJETCANDIDATES]; 
    float pthat, PVz, Jet_pt[MAXQCDJETCANDIDATES], Jet_eta[MAXQCDJETCANDIDATES], Jet_phi[MAXQCDJETCANDIDATES], Jet_Ip1P[MAXQCDJETCANDIDATES], Jet_Ip2P[MAXQCDJETCANDIDATES], Jet_Ip3P[MAXQCDJETCANDIDATES];
    float Jet_Svx[MAXQCDJETCANDIDATES], Jet_SvxHP[MAXQCDJETCANDIDATES], Jet_ProbaP[MAXQCDJETCANDIDATES], Jet_ProbaN[MAXQCDJETCANDIDATES], Jet_Proba[MAXQCDJETCANDIDATES];
    float Jet_BprobP[MAXQCDJETCANDIDATES], Jet_BprobN[MAXQCDJETCANDIDATES], Jet_Bprob[MAXQCDJETCANDIDATES], Jet_CombSvxP[MAXQCDJETCANDIDATES], Jet_CombSvxN[MAXQCDJETCANDIDATES], Jet_CombSvx[MAXQCDJETCANDIDATES];
    float Muon_pt[MAXQCDJETCANDIDATES], Muon_eta[MAXQCDJETCANDIDATES], Muon_ptrel[MAXQCDJETCANDIDATES], Muon_IP[MAXQCDJETCANDIDATES], Muon_IPsig[MAXQCDJETCANDIDATES], Muon_Proba[MAXQCDJETCANDIDATES], Muon_chi2[MAXQCDJETCANDIDATES], Muon_chi2Tk[MAXQCDJETCANDIDATES]; 
    int Muon_nMuHit[MAXQCDJETCANDIDATES], Muon_nTkHit[MAXQCDJETCANDIDATES], Muon_nPixHit[MAXQCDJETCANDIDATES], Muon_nOutHit[MAXQCDJETCANDIDATES], Muon_isGlobal[MAXQCDJETCANDIDATES], Muon_nMatched[MAXQCDJETCANDIDATES];
    float bFromGSplit_pT[MAXQCDJETCANDIDATES], bFromGSplit_eta[MAXQCDJETCANDIDATES], bFromGSplit_phi[MAXQCDJETCANDIDATES];
    float TrkInc_pt[MAXQCDTRACKSINJETS], TrkInc_ptrel[MAXQCDTRACKSINJETS], TrkInc_IP[MAXQCDTRACKSINJETS], TrkInc_IPsig[MAXQCDTRACKSINJETS];
  };
  
  class EventSummaryHandler{
  public:
    
    //c/dtor
    EventSummaryHandler();
    ~EventSummaryHandler();
    
    //current event
    EventSummary_t evSummary_;
    EventSummary_t &getEvent() { return evSummary_; }

    //read mode
    bool attachToTree(TTree *t);
    inline void getEntry(int ientry) { if(t_) t_->GetEntry(ientry); }
    inline int getEntries() { return (t_ ? t_->GetEntriesFast() : 0); }
    inline TTree *getTree(){return t_;}
    
    inline void resetTree() 
      { 
	if(!t_) return;
	t_->Delete(); 
	t_=0;
      }
    
  private:
    
    //the tree
    TTree *t_;
    
  };

  //
  class PhysicsObject_Track
  {
  public:
    PhysicsObject_track(float pt_, float ptrel_, float ip_, float ipsig_) :
      pt(pt_),ptrel(ptrel_),ip(ip_),ipsig(ipsig_)
      {
      }
      float pt,ptrel,ip,ipsig;
  }

  //
  class PhysicsObject_SoftMuon 
  {
  public :
    PhysicsObject_Muon(float pt_, float eta_, float ptrel_, float ip_, float ipsig_)
      : pt(pt_), eta(eta_), ptrel(ptrel_), ip(ip_), ipsig(ipsig_) { }
      float pt,eta,ptrel,ip,ipsig;
      int prob,chi2,tk_chi2,n_muHits, n_pxHits, n_tkHits, n_outHits, n_matchedStations;
      bool isGlobal;

  };

  //
  class PhysicsObject_Jet : public LorentzVector
  {
  public :
    PhysicsObject_Jet(LorentzVector vec);
    int flavour;
    float ssv,  ssvhp, jp, jp_p, jp_n, jpb, jbp_p, jbp_n, csv, csv_p, csv_n;
    int fTrack,lTrack, mid;
    float pt_gsplit,eta_gsplit,phi_gsplit;
    std::vector<PhysicsObject_Track> m_tracks;
    std::vector<PhysicsObject_SoftMuon> m_muons;
  };
  
  typedef std::vector<PhysicsObject_Jet>    PhysicsObjectJetCollection;
  
  //
  class QCDEvent_t
  {
  public:
    PhysicsEvent_t() {};
    ~PhysicsEvent_t() {};

    PhysicsObjectJetCollection jets;
    static bool sortJetsByCSV(PhysicsObject_Jet a,PhysicsObject_Jet b)   {   return (a.csv>b.csv);  }  
    static bool sortJetsByJP(PhysicsObject_Jet a,PhysicsObject_Jet b)   {   return (a.jp>b.jp);  }  
    static bool sortJetsByPt(PhysicsObject_Jet a,PhysicsObject_Jet b)     {   return (a.pt()>b.pt());  }  
  };
  
  //
  PhysicsEvent_t getPhysicsEventFrom(EventSummary_t &ev);
}

  
#endif
