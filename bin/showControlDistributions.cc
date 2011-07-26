#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "LIP/Top/interface/EventSummaryHandler.h"

#include "CMGTools/HtoZZ2l2nu/interface/BtagUncertaintyComputer.h"
#include "CMGTools/HtoZZ2l2nu/interface/JetEnergyUncertaintyComputer.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"
#include "CMGTools/HtoZZ2l2nu/interface/TransverseMassComputer.h"

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
  int mcTruthMode = runProcess.getParameter<int>("mctruthmode");
  int evStart=runProcess.getParameter<int>("evStart");
  int evEnd=runProcess.getParameter<int>("evEnd");
  TString dirname = runProcess.getParameter<std::string>("dirName");

  //start the uncertainty computers
  btag::UncertaintyComputer bcomp(0.837, 0.95, 0.06,
 				  0.286, 1.11, 0.11);
//   btag::UncertaintyComputer bcomp(0.81, 1.0, 0.06,
// 				  0.25, 1.0, 0.11);

  TString etaFileName = runProcess.getParameter<std::string>("etaResolFileName"); gSystem->ExpandPathName(etaFileName);
  JetResolution stdEtaResol(etaFileName.Data(),false);
  TString phiFileName = runProcess.getParameter<std::string>("phiResolFileName"); gSystem->ExpandPathName(phiFileName);
  JetResolution stdPhiResol(phiFileName.Data(),false);
  TString ptFileName  = runProcess.getParameter<std::string>("ptResolFileName");  gSystem->ExpandPathName(ptFileName);
  JetResolution stdPtResol(ptFileName.Data(),true);
  TString uncFile =  runProcess.getParameter<std::string>("jesUncFileName"); gSystem->ExpandPathName(uncFile);
  JetCorrectionUncertainty jecUnc(uncFile.Data());
  jet::UncertaintyComputer jcomp(&stdPtResol,&stdEtaResol,&stdPhiResol,&jecUnc);
  
  //control histograms
  SelectionMonitor controlHistos;
  double massAxis[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,400,500,1000,2000};
  int nMassBins=sizeof(massAxis)/sizeof(double)-1;
  controlHistos.addHistogram( new TH1D("mlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV/c^{2}];Lepton-jet pairs",nMassBins,massAxis) );
  controlHistos.addHistogram( new TH1D("btagdisc",";TCHE (b-jets);Events",100,-5,40) );
  controlHistos.addHistogram( new TH1D("ltagdisc",";TCHE (udcsg-jets);Events",100,-5,40) );
  controlHistos.addHistogram( new TH1D("dilmass",";M(l,l');Events",100,0,250) );
  controlHistos.addHistogram( new TH1D("mtsum",";M_{T}(l^{(1)},E_{T}^{miss})+M_{T}(l^{(2)},E_{T}^{miss});Events",100,0,1000) );
  controlHistos.addHistogram( new TH1F ("leadjet", "; Leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadjet", "; Sub-leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("leadlepton", "; Leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadlepton", "; Sub-leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("met", "; #slash{E}_{T} [GeV/c]; Events / (10 GeV/c)", 40, 0.,400.) );
  TH1F *bmultH=new TH1F ("btags", ";b-tag multiplicity;Events", 6, 0.,6.) ;
  TH1F *bmultsfH=new TH1F ("btagssf", ";b-tag multiplicity;Events", 6, 0.,6.) ;
  for(int ibin=1; ibin<=bmultH->GetXaxis()->GetNbins(); ibin++)
    {
      TString label(""); label += ibin-1;
      if(ibin==bmultH->GetXaxis()->GetNbins()) label ="#geq" + label;
      bmultH->GetXaxis()->SetBinLabel(ibin,label + " btags");
      bmultsfH->GetXaxis()->SetBinLabel(ibin,label + " btags");
    }
  controlHistos.addHistogram(bmultH);
  controlHistos.addHistogram(bmultsfH);

  //cutflow histograms
  TString labels[]={"start","#geq 2 jets","=0 b-tags", "=1 b-tags","=2 b-tags"};
  int nsteps=sizeof(labels)/sizeof(TString);
  TH1D *cutflowH=new TH1D("evtflow",";Cutflow;Events",nsteps,0,nsteps);
  for(int ibin=0; ibin<nsteps; ibin++) cutflowH->GetXaxis()->SetBinLabel(ibin+1,labels[ibin]);
  controlHistos.addHistogram( cutflowH );
  TString cats[]={"","jer","jesdown","jesup","btagcen","btagup","btagdown"};
  int nvarcats=sizeof(cats)/sizeof(TString);
  for(int icat=0;icat<nvarcats; icat++)  controlHistos.addHistogram( new TH1D("evtflow"+cats[icat],";Cutflow;Events",nsteps,0,nsteps) );

  controlHistos.initMonitorForStep("ee");
  controlHistos.initMonitorForStep("emu");
  controlHistos.initMonitorForStep("mumu");
  controlHistos.initMonitorForStep("etau");
  controlHistos.initMonitorForStep("mutau");
  
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

  //aux
  TransverseMassComputer mtComp;
 
  //loop over events
  if(evEnd<0 || evEnd>evSummaryHandler.getEntries() ) evEnd=evSummaryHandler.getEntries();
  if(evStart > evEnd )
    {
      evfile->Close();
      return -1;
    }

  cout << endl;
  for (int inum=evStart; inum < evEnd; ++inum)
  {
    if(inum%500==0) printf("\r [ %d/100 ] %s",int(100*float(inum-evStart)/float(evEnd)),evurl.Data());

    evTree->GetEvent(inum);
    EventSummary_t &ev = evSummaryHandler.getEvent();
    if(isMC)
      {
        if(mcTruthMode==1 && !ev.isSignal) continue;
        if(mcTruthMode==2 && ev.isSignal) continue;
      }

    float weight = ev.weight;    
   
    TString ch("");
    if(ev.cat==dilepton::MUMU)  ch="mumu";
    if(ev.cat==dilepton::EE)  ch="ee";
    if(ev.cat==dilepton::EMU) ch="emu";
    if(ev.cat==dilepton::ETAU) ch="etau";
    if(ev.cat==dilepton::MUTAU) ch="mutau";
    TString catsToFill[]={"all",ch};

    //get particles from the event
    double btagcut(1.7); 
    int njets(0),nbtags(0),nbjets(0),nljets(0);
    LorentzVectorCollection leptons, jets, mets;
      for(Int_t ipart=0; ipart<ev.nparticles; ipart++)
      {
	LorentzVector p4(ev.px[ipart],ev.py[ipart],ev.pz[ipart],ev.en[ipart]);
	if(isnan(p4.pt()) || isinf(p4.pt())) continue;
	switch( ev.id[ipart] )
	  {
	  case 0:
	    mets.push_back( p4 );
            break;
          case 1:
	    jets.push_back( p4 );
	    njets++;
	    if(p4.pt()>30)
	      {
		nbtags += (ev.info1[ipart]>btagcut);
		nbjets +=(fabs(ev.genid[ipart])==5);
		nljets +=(fabs(ev.genid[ipart])!=5);
		if(isMC)
		  {
		    controlHistos.fillHisto( fabs(ev.genid[ipart])==5 ? "btagdisc" : "ltagdisc",catsToFill[0],ev.info1[ipart],weight);
		    controlHistos.fillHisto( fabs(ev.genid[ipart])==5 ? "btagdisc" : "ltagdisc",catsToFill[1],ev.info1[ipart],weight);
		  }
	      }
	    break;
          default:
            leptons.push_back( p4 );
            break;
	  }
      }
      LorentzVector dilepton = leptons[0]+leptons[1];
      
      //jet variations
      jcomp.compute(jets,mets[0]);
      LorentzVectorCollection jerVariedJets=jcomp.getVariedJets(jet::UncertaintyComputer::JER);
      LorentzVectorCollection jesupVariedJets=jcomp.getVariedJets(jet::UncertaintyComputer::JES_UP);
      LorentzVectorCollection jesdownVariedJets=jcomp.getVariedJets(jet::UncertaintyComputer::JES_DOWN);

      double mtsum=
	mtComp.compute(leptons[0],mets[0]) +
	mtComp.compute(leptons[1],mets[0]);

      
    //btag variations
    bcomp.compute(nbjets,nljets);
    std::vector<btag::Weight_t> wgt = bcomp.getWeights();
    double p0btags = wgt[0].first;  double p0btags_err=wgt[0].second;
    double p1btags = wgt[1].first;  double p1btags_err=wgt[1].second;
    double p2btags = wgt[2].first;  double p2btags_err=wgt[2].second;

    //the cutflow
    for(int ivar=0;ivar<nvarcats; ivar++) 
      {
	for(size_t ictf=0; ictf<2; ictf++)
	  {
	    controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],0,weight);

	    //jet energy variations
	    LorentzVectorCollection jetColl=jets;
	    if(cats[ivar]=="jer") jetColl=jerVariedJets;
	    if(cats[ivar]=="jesdown") jetColl=jesdownVariedJets;
	    if(cats[ivar]=="jesup" ) jetColl=jesupVariedJets;
	    int nseljets(0);
	    for(size_t ijet=0; ijet<jetColl.size(); ijet++) nseljets += (jetColl[ijet].pt()>30);
	    if(nseljets<2) continue;
	    
	    if(cats[ivar]=="")
	      {
		controlHistos.fillHisto("mtsum",catsToFill[ictf],mtsum,weight);
		controlHistos.fillHisto("leadjet",catsToFill[ictf],max(jets[0].pt(),jets[1].pt()),weight);
		controlHistos.fillHisto("subleadjet",catsToFill[ictf],min(jets[0].pt(),jets[1].pt()),weight);
		controlHistos.fillHisto("leadlepton",catsToFill[ictf],max(leptons[0].pt(),leptons[1].pt()),weight);
		controlHistos.fillHisto("subleadlepton",catsToFill[ictf],min(leptons[0].pt(),leptons[1].pt()),weight);
		controlHistos.fillHisto("met",catsToFill[ictf],mets[0].pt(),weight);
		controlHistos.fillHisto("dilmass",catsToFill[ictf],dilepton.mass(),weight);
		controlHistos.fillHisto("btags",catsToFill[ictf],nbtags,weight);
		if(isMC)
		  {
		    controlHistos.fillHisto("btagssf",catsToFill[ictf],0,p0btags*weight);
		    controlHistos.fillHisto("btagssf",catsToFill[ictf],1,p1btags*weight);
		    controlHistos.fillHisto("btagssf",catsToFill[ictf],2,p2btags*weight);
		  }
	      }
	    controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],1,weight);
	    
	    //b-tag variations
	    if(isMC && cats[ivar].Contains("btag"))
	      {
		double p0weight(p0btags), p1weight(p1btags), p2weight(p2btags);
		if(cats[ivar]=="btagup") { 
		  p0weight += p0btags_err; p1weight += p1btags_err; p2weight += p2btags_err;   
		}
		else if(cats[ivar]=="btagdown") { 
		  p0weight = p0weight-p0btags_err; p1weight = p1weight-p1btags_err; p2weight = p2weight-p2btags_err;
		}
		controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],2,weight*p0weight);
		controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],3,weight*p1weight);
		controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],4,weight*p2weight);
	      }
	    
	    //standard b-tag
	    else
	      {
		int btagbin(nbtags);
		if(btagbin>2) btagbin=2;
		controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],2+btagbin,weight);
	      }
	  }
      }
    
    //fill histos
    for(int ijet=0; ijet<njets; ijet++)
      {
	//get the lepton-jet pairs
	LorentzVector lj1=leptons[0]+jets[ijet];
	float mlj1=lj1.mass();

	LorentzVector lj2=leptons[1]+jets[ijet];
	float mlj2=lj2.mass();

	for(int ictf=0; ictf<2; ictf++)
	  {
	    controlHistos.fillHisto("mlj",catsToFill[ictf],mlj1,weight,true);
	    controlHistos.fillHisto("mlj",catsToFill[ictf],mlj2,weight,true);
	  }
      }
  }
  
  //overall normalization factor
  float cnorm=1.0;
  if(isMC)
    {
      TH1F *cutflowH = (TH1F *) evfile->Get("evAnalyzer/top/cutflow");
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
	  //	  if(hit->first.BeginsWith("evtflow") && hit->first!="evtflow")
	  //	    {
	      //	      hit->second->Add( controlHistos.getHisto("evtflow","all"), -1);
	  //	    }
	  hit->second->Write();


 }
    }
  
  file->Close(); 
}  
