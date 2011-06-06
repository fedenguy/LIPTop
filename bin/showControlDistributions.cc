#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TLorentzVector.h"
#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCanvas.h"
#include "TString.h"
#include "TDirectory.h"

using namespace std;

//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  //check arguments
  if ( argc < 4 ) {
    std::cout << "Usage : " << argv[0] << " eventsFile outputDir isMC" << std::endl;
    return 0;
  }

  //check if this is MC
  TString isMCBuf(argv[3]);
  bool isMC=isMCBuf.Atoi();

  std::map<TString, TH1 *> results;
  TString cats[]={"all","ee","mumu","emu"};
  size_t ncats=sizeof(cats)/sizeof(TString);
  double massAxis[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,400,500,1000,2000};
  int nMassBins=sizeof(massAxis)/sizeof(double)-1;
  for( size_t icat=0; icat<ncats; icat++)
    {
      results[cats[icat]+"_mlj"] = (TH1*) formatPlot(new TH1D(cats[icat]+"_mlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV/c^{2}];Lepton-jet pairs",nMassBins,massAxis) , 1,1,1,20,0,true,true,1,1,1);
    }

  //process kin file
  TString evurl = argv[1];
  gSystem->ExpandPathName(evurl);
  cout << evurl << endl;
  //process events file
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

 
  //loop over events
  for (int inum=0; inum < evTree->GetEntriesFast(); ++inum)
  {
    evTree->GetEvent(inum);
    EventSummary_t &ev = evSummaryHandler.getEvent();

    //classify event
    std::vector<TString> categs;
    categs.push_back("all");
    if(ev.cat==dilepton::MUMU)  categs.push_back("mumu");
    if(ev.cat==dilepton::EE)  categs.push_back("ee");
    if(ev.cat==dilepton::EMU)  categs.push_back("emu");
    
    //get particles from the event
    int njets(0),nbtags(0);
    std::vector<TLorentzVector> leptons, jets, mets;
    for(Int_t ipart=0; ipart<ev.nparticles; ipart++)
      {
	TLorentzVector p4(ev.px[ipart],ev.py[ipart],ev.pz[ipart],ev.en[ipart]);
	if(isnan(p4.Pt()) || isinf(p4.Pt())) continue;
	switch( ev.id[ipart] )
	  {
	  case 0:
	    mets.push_back( p4 );
            break;
          case 1:
            jets.push_back( p4 );
	    njets++;
	    if(ev.info1[ipart]>1.7) nbtags++;
            break;
          default:
            leptons.push_back( p4 );
            break;

	  }
      }

    //fill histos
    float weight = ev.weight;    
    for(std::vector<TString>::iterator cIt = categs.begin(); cIt != categs.end(); cIt++)
      {
	for(int ijet=0; ijet<njets; ijet++)
	  {
	    //get the lepton-jet pairs
	    TH1D *mljH=(TH1D *)results[*cIt+"_mlj"];
	    TLorentzVector lj1=leptons[0]+jets[ijet];
	    float mlj1=lj1.M();
	    int ibin = mljH->FindBin(mlj1);
	    float width = mljH->GetBinWidth(ibin);
	    mljH->Fill(mlj1,weight/width);

	    TLorentzVector lj2=leptons[1]+jets[ijet];
	    float mlj2=lj2.M();
	    ibin = mljH->FindBin(mlj2);
	    width = mljH->GetBinWidth(ibin);
	    mljH->Fill(mlj2,weight/width);
	  }
      }
  }
    
  //save to file
  TString outUrl = argv[2];
  gSystem->ExpandPathName(outUrl);
  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  outUrl += gSystem->BaseName(evurl);
  TFile *file=TFile::Open(outUrl, "recreate");

  float cnorm=1.0;
  if(isMC)
    {
      TString tag=gSystem->BaseName(evurl);
      tag.ReplaceAll(".root","");
      TH1F *cutflowH = (TH1F *) evfile->Get("evAnalyzer/"+tag+"/cutflow");
      if(cutflowH)
	cnorm=cutflowH->GetBinContent(1);
    }

  TDirectory *baseOutDir=file->mkdir("ctrlAnalyzer");
  for( size_t icat=0; icat<ncats; icat++)
    {
      baseOutDir->cd();
      if(icat) baseOutDir->mkdir( cats[icat] )->cd();
      for(std::map<TString,TH1 *>::iterator hIt = results.begin(); hIt != results.end(); hIt++) 
	{
	  if(!hIt->first.BeginsWith(cats[icat])) continue;
	  if(isMC && cnorm>0) hIt->second->Scale(1./cnorm);
	  hIt->second->Write();
	}
    }

  file->Close(); 
}  
