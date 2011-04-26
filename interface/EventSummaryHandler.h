
#ifndef eventsummaryhandler_h
#define eventsummaryhandler_h

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>

#include "TTree.h"

#endif

#define MAXPARTICLES 25
struct EventSummary_t
{
  Int_t run,lumi,event;
  Int_t cat;
  Float_t weight;
  Int_t nparticles;
  Float_t px[MAXPARTICLES], py[MAXPARTICLES], pz[MAXPARTICLES], en[MAXPARTICLES];
  Int_t id[MAXPARTICLES], genid[MAXPARTICLES];
  Float_t info1[MAXPARTICLES],info2[MAXPARTICLES],info3[MAXPARTICLES],info4[MAXPARTICLES],info5[MAXPARTICLES];
};

class EventSummaryHandler{
 public:
  
  EventSummary_t evSummary_;
  
  EventSummaryHandler();
  void initTree(TTree *t);
  void getEntry(int ientry) { if(t_) t_->GetEntry(ientry); }
  int getEntries() { return (t_ ? t_->GetEntriesFast() : 0); }
  EventSummary_t &getEvent() { return evSummary_; }
  void fillTree();
  ~EventSummaryHandler();

 private:
  
  TTree *t_;
};

#endif
