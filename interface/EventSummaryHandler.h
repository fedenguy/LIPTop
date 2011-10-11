
#ifndef eventsummaryhandler_h
#define eventsummaryhandler_h

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

#define MAXPARTICLES 25
#define MAXMEASUREMENTS 10
struct EventSummary_t
{
  Int_t run,lumi,event;
  Int_t cat;
  Bool_t isSignal,hasTrigger;
  Int_t nvtx,ngenpu;
  Float_t rho,weight,normWeight;
  Int_t nparticles,nmcparticles;
  Float_t px[MAXPARTICLES], py[MAXPARTICLES], pz[MAXPARTICLES], en[MAXPARTICLES];
  Int_t id[MAXPARTICLES], genid[MAXPARTICLES], genflav[MAXPARTICLES];
  Float_t info1[MAXPARTICLES],info2[MAXPARTICLES],info3[MAXPARTICLES],info4[MAXPARTICLES],info5[MAXPARTICLES],
    info6[MAXPARTICLES],info7[MAXPARTICLES],info8[MAXPARTICLES],info9[MAXPARTICLES];
  Float_t mcpx[MAXPARTICLES], mcpy[MAXPARTICLES], mcpz[MAXPARTICLES], mcen[MAXPARTICLES];
  Int_t mcid[MAXPARTICLES];
  Int_t nmeasurements;
  Float_t evmeasurements[MAXMEASUREMENTS];
};

class EventSummaryHandler{
 public:

  //c/dtor
  EventSummaryHandler();
  ~EventSummaryHandler();
  
  //current event
  EventSummary_t evSummary_;
  EventSummary_t &getEvent() { return evSummary_; }
  
  //write mode
  bool initTree(TTree *t,bool needsToRecreate=true);
  void fillTree();
  void fillTreeWithEvent(const EventSummary_t &ev);
  void fillTreeWithEvent(const EventSummary_t &ev,std::vector<float> &addMeasurements);

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
class PhysicsObject : public LorentzVector
{
public :
  PhysicsObject(LorentzVector vec, Int_t id_):
    LorentzVector(vec), id(id_){ }
    Int_t id;
};


//
class PhysicsObject_Lepton : public LorentzVector
{
public :
  PhysicsObject_Lepton(LorentzVector vec, Int_t id_,Int_t genid_=0, Float_t ptErr_=0, Float_t iso1_=0, Float_t iso2_=0, Float_t iso3_=0):
    LorentzVector(vec), id(id_), genid(genid_), ptErr(ptErr_), iso1(iso1_), iso2(iso2_), iso3(iso3_) { }
    Int_t id,genid;
    Float_t ptErr, iso1, iso2, iso3;
};

class PhysicsObject_Jet : public LorentzVector
{
public :
  PhysicsObject_Jet(LorentzVector vec, Int_t genid_=0, Int_t flavid_=0,
		    Float_t btag1_=0, Float_t btag2_=0, Float_t btag3_=0, Float_t btag4_=0, Float_t btag5_=0, Float_t btag6_=0, Float_t btag7_=0,
		    Float_t neutHadEnFrac_=0, Float_t emEnFrac_=0,
		    Bool_t vtxAssoc_=true):
    LorentzVector(vec), genid(genid_), flavid(flavid_), btag1(btag1_), btag2(btag2_), btag3(btag3_), btag4(btag4_), 
    btag5(btag5_), btag6(btag6_), btag7(btag7_), neutHadEnFrac(neutHadEnFrac_), emEnFrac(emEnFrac_), vtxAssoc(vtxAssoc_) { }
    Int_t genid,flavid;
    Float_t btag1, btag2, btag3, btag4,btag5,btag6,btag7; 
    Float_t neutHadEnFrac, emEnFrac;
    Bool_t vtxAssoc;
};

typedef std::vector<PhysicsObject>        PhysicsObjectCollection;
typedef std::vector<PhysicsObject_Lepton> PhysicsObjectLeptonCollection;
typedef std::vector<PhysicsObject_Jet>    PhysicsObjectJetCollection;

//
struct PhysicsEvent_t
{
  LorentzVector met;
  LorentzVector vtx;
  PhysicsObjectJetCollection jets;
  Int_t nbjets,nljets;
  LorentzVector dil;
  PhysicsObjectLeptonCollection leptons;
  LorentzVector top,antitop;
  PhysicsObjectCollection topdecay,antitopdecay;
};

//
PhysicsEvent_t getPhysicsEventFrom(EventSummary_t &ev);


#endif
