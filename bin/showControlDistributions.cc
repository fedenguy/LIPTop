#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"

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
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
    return 0;
  }

  //configure                                                                                                                                                                                                                                   
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  TString evurl=runProcess.getParameter<std::string>("input");
  TString outdir=runProcess.getParameter<std::string>("outdir");
  bool isMC = runProcess.getParameter<bool>("isMC");
  int evStart=runProcess.getParameter<int>("evStart");
  int evEnd=runProcess.getParameter<int>("evEnd");
  TString dirname = runProcess.getParameter<std::string>("dirName");

  //control histograms
  SelectionMonitor controlHistos;
  double massAxis[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,400,500,1000,2000};
  int nMassBins=sizeof(massAxis)/sizeof(double)-1;
  controlHistos.addHistogram( new TH1D("mlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV/c^{2}];Lepton-jet pairs",nMassBins,massAxis) );
  controlHistos.addHistogram( new TH1D("cutflow",";Cutflow;Events",3,0,3) );
  TString cats[]={"","up","down"};
  for(size_t icat=0;icat<sizeof(cats)/sizeof(TString); icat++)
    controlHistos.addHistogram( new TH1D("cutflow"+cats[icat],";Cutflow;Events",3,0,3) );
   
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
  if(evEnd<0 || evEnd>evSummaryHandler.getEntries() ) evEnd=evSummaryHandler.getEntries();
  if(evStart > evEnd )
    {
      evfile->Close();
      return -1;
    }
  for (int inum=evStart; inum < evEnd; ++inum)
  {
    evTree->GetEvent(inum);
    EventSummary_t &ev = evSummaryHandler.getEvent();

    //classify event
    if(ev.cat!=dilepton::EMU)  continue;
    
    //get particles from the event
    int njets(0);
    int nbtags[]={0.,0.,0.};
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
	    if(ev.info1[ipart]>1.7) nbtags[0]++;
	    if(ev.info1[ipart]>2.2) nbtags[1]++;
	    if(ev.info1[ipart]>1.2) nbtags[2]++;
	    break;
          default:
            leptons.push_back( p4 );
            break;
	  }
      }

    //fill histos
    float weight = ev.weight;    
    for(int ijet=0; ijet<njets; ijet++)
      {
	//get the lepton-jet pairs
	
	TLorentzVector lj1=leptons[0]+jets[ijet];
	float mlj1=lj1.M();
	controlHistos.fillHisto("mlj","all",mlj1,weight,true);

	TLorentzVector lj2=leptons[1]+jets[ijet];
	float mlj2=lj2.M();
	controlHistos.fillHisto("mlj","all",mlj2,weight,true);
      }
   
    //the cutflow
    for(size_t ivar=0;ivar<3; ivar++) 
      {
	if(nbtags[ivar]>2) nbtags[ivar]=2;
	controlHistos.fillHisto("cutflow"+cats[ivar],"all",nbtags[ivar],weight);
      } 
 }

  //overall normalization factor
  float cnorm=1.0;
  if(isMC)
    {
      TString tag=gSystem->BaseName(evurl);
      tag.ReplaceAll(".root","");
      TH1F *cutflowH = (TH1F *) evfile->Get("evAnalyzer/"+tag+"/cutflow");
      if(cutflowH) cnorm=cutflowH->GetBinContent(1);
    }

  //all done with the events file
  evfile->Close();

    
  //save to file
  TString outUrl(outdir);
  gSystem->ExpandPathName(outUrl);
  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  outUrl += gSystem->BaseName(evurl);
  TFile *file=TFile::Open(outUrl, "recreate");
  TDirectory *baseOutDir=file->mkdir("localAnalysis");
  SelectionMonitor::StepMonitor_t &mons=controlHistos.getAllMonitors();
  std::map<TString, TDirectory *> outDirs;
  outDirs["all"]=baseOutDir->mkdir("all");
  outDirs["ee"]=baseOutDir->mkdir("ee");
  outDirs["emu"]=baseOutDir->mkdir("emu");
  outDirs["mumu"]=baseOutDir->mkdir("mumu");
  for(SelectionMonitor::StepMonitor_t::iterator it =mons.begin(); it!= mons.end(); it++)
    {
      TString icat("all");
      if(it->first.BeginsWith("ee")) icat="ee";
      if(it->first.BeginsWith("emu")) icat="emu";
      if(it->first.BeginsWith("mumu")) icat="mumu";
      outDirs[icat]->cd();
      for(SelectionMonitor::Monitor_t::iterator hit=it->second.begin(); hit!= it->second.end(); hit++)
	{
          if(isMC && cnorm>0) hit->second->Scale(1./cnorm);
          if( !((TClass*)hit->second->IsA())->InheritsFrom("TH2")
              && !((TClass*)hit->second->IsA())->InheritsFrom("TGraph") )
            fixExtremities(hit->second,true,true);
	  if(hit->first.BeginsWith("cutflow") && hit->first!="cutflow")
	    {
	      hit->second->Add( controlHistos.getHisto("cutflow","all"), -1);
	    }
	  hit->second->Write();
        }
    }
  
  file->Close(); 
}  
