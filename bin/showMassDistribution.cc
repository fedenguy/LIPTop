#include <iostream>
#include <boost/shared_ptr.hpp>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/KinResultsHandler.h"
#include "LIP/Top/interface/KinAnalysis.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCanvas.h"

using namespace std;

//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  //check arguments
  if ( argc < 5 ) {
    std::cout << "Usage : " << argv[0] << " kinDir eventsFile outputDir fixEntries" << std::endl;
    return 0;
  }

  TH1F *hmtop=(TH1F*)formatPlot( new TH1F ("hmtop", "; m_{Top} [GeV/c^{2}]; Events", 100, 0.,500.), 1,1,1,20,0,true,true,1,1,1);
  TH2F *hmtopvsdilmass=(TH2F*)formatPlot( new TH2F ("mtopvsdilmass", "; m_{Top} [GeV/c^{2}]; Mass(l,l') [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,500.), 1,1,1,20,0,true,true,1,1,1);

  //fix entries flag
  TString fixEntriesBuf(argv[4]);
  bool fixEntries=fixEntriesBuf.Atoi();

  //process kin file
  TString url = argv[1];
  KinResultsHandler kinHandler;
  kinHandler.init(url,false);
  TChain *t=kinHandler.getResultsChain();

  //process events file
  TString evurl = argv[2];
  TFile *evfile = TFile::Open(evurl);
  if(evfile==0) return -1;
  if(evfile->IsZombie()) return -1;
  EventSummaryHandler evSummaryHandler;
  if( !evSummaryHandler.attachToTree( (TTree *)evfile->Get("evAnalyzer/data") ) ) 
    {
      evfile->Close();
      return -1;
    }  
  TTree *evTree=evSummaryHandler.getTree();


 /*
 // tree synchronization
  evTree->BuildIndex ("run","event");
  t->BuildIndex ("run","event");
  evTree->AddFriend(t);
*/
 
  //loop over events
  std::map<TString,int> selEvents;
  for (int inum=0; inum < evTree->GetEntriesFast(); ++inum)
  {
    evTree->GetEvent(inum);
    EventSummary_t &ev = evSummaryHandler.getEvent();
    TString key(""); key+= ev.run; key+="-"; key += ev.lumi; key+="-"; key += ev.event;
    selEvents[key]=inum;
  }

  //loop over kin results
  int nresults(0);
  for (int inum=0; inum < t->GetEntries(); ++inum){
    t->GetEvent(inum);
 
    //get original event
    Int_t irun,ievent,ilumi;
    kinHandler.getEventInfo(irun,ievent,ilumi);
    //fix me after reprocessing KIN trees event <-> lumi
    TString key("");  key+= irun; key+="-"; key += ievent; key+="-"; key += ilumi;
    if(selEvents.find(key)==selEvents.end()) continue;
    nresults++;
    evTree->GetEntry( selEvents[key] );
    EventSummary_t &ev = evSummaryHandler.getEvent();

    //get particles from the event
    std::vector<TLorentzVector> leptons, jets, mets;
    for(Int_t ipart=0; ipart<ev.nparticles; ipart++)
      {
	TLorentzVector p4(ev.px[ipart],ev.py[ipart],ev.pz[ipart],ev.en[ipart]);
	switch( ev.id[ipart] )
	  {
	  case 0:
	    mets.push_back(p4);
	    break;
	  case 1:
	    jets.push_back(p4);
	    break;
	  default:
	    leptons.push_back(p4);
	    break;
	  }
      }

    //compute dilepton invariant mass
    TLorentzVector dil = leptons[0]+leptons[1];
    float dilmass = dil.M();


    //fill histos
    TH1F *h1=kinHandler.getHisto("mt",1), *h2=kinHandler.getHisto("mt",2);
    TH1F *hprev=h1->Integral()> h2->Integral()? h1:h2;
    vector<double> res=kinHandler.getMPVEstimate(hprev);
    double mtop=res[1];

    hmtop->Fill(mtop);
    hmtopvsdilmass->Fill(mtop,dilmass);

    if (mtop>350) cout << irun << ":" << ilumi << ":" << ievent << " " << flush;
  }
  kinHandler.end();
  cout << endl;

  //rescale to number of selected events
  double scaleFactor=selEvents.size()/nresults;
  if(fixEntries)
  {
    hmtop->Scale(scaleFactor);
    hmtopvsdilmass->Scale(scaleFactor);
  }

  //save to file
  TString outUrl( argv[3] );
  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  outUrl += gSystem->BaseName(evurl);
  TFile *file=TFile::Open(outUrl, "recreate");
  hmtop->Write();
  hmtopvsdilmass->Write();
  file->Close(); 
}  
