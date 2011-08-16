
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
#include "CMGTools/HtoZZ2l2nu/interface/DuplicatesChecker.h"

#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCanvas.h"
#include "TString.h"
#include "TDirectory.h"
#include "TEventList.h"

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

  //
  // configure
  //
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  TString evurl=runProcess.getParameter<std::string>("input");
  TString outdir=runProcess.getParameter<std::string>("outdir");
  bool isMC = runProcess.getParameter<bool>("isMC");
  int mcTruthMode = runProcess.getParameter<int>("mctruthmode");
  int evStart=runProcess.getParameter<int>("evStart");
  int evEnd=runProcess.getParameter<int>("evEnd");
  TString dirname = runProcess.getParameter<std::string>("dirName");
  double xsec = runProcess.getParameter<double>("xsec");
  bool saveSummaryTree = runProcess.getParameter<bool>("saveSummaryTree");
  TString etaFileName = runProcess.getParameter<std::string>("etaResolFileName"); gSystem->ExpandPathName(etaFileName);
  TString phiFileName = runProcess.getParameter<std::string>("phiResolFileName"); gSystem->ExpandPathName(phiFileName);
  TString ptFileName  = runProcess.getParameter<std::string>("ptResolFileName");  gSystem->ExpandPathName(ptFileName);
  TString uncFile =  runProcess.getParameter<std::string>("jesUncFileName"); gSystem->ExpandPathName(uncFile);

  // b-tag working points: https://twiki.cern.ch/twiki/bin/view/CMS/BTagPerformanceOP
  std::map<TString,float> btagCuts;
  btagCuts["TCHEL"]=1.7;
  btagCuts["TCHEM"]=3.3;
  btagCuts["TCHPT"]=3.41;
  btagCuts["JBPL"]=1.33;
  btagCuts["JBPM"]=2.55;
  btagCuts["JBPT"]=3.74;
  btagCuts["SSVHEM"]=1.74;
   
  //
  // start auxiliary computers
  //
  TransverseMassComputer mtComp;
  btag::UncertaintyComputer bcomp(0.837, 0.95, 0.06, 0.286, 1.11, 0.11);
  JetResolution stdEtaResol(etaFileName.Data(),false);
  JetResolution stdPhiResol(phiFileName.Data(),false);
  JetResolution stdPtResol(ptFileName.Data(),true);
  JetCorrectionUncertainty jecUnc(uncFile.Data());
  jet::UncertaintyComputer jcomp(&stdPtResol,&stdEtaResol,&stdPhiResol,&jecUnc);
  
  //
  // control histograms
  //
  SelectionMonitor controlHistos;
  double massAxis[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,400,500,1000,2000};
  int nMassBins=sizeof(massAxis)/sizeof(double)-1;
  controlHistos.addHistogram( new TH1D("mlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV/c^{2}];Lepton-jet pairs",nMassBins,massAxis) );
  controlHistos.addHistogram( new TH1D("dilmass",";M(l,l');Events",100,0,250) );
  controlHistos.addHistogram( new TH1D("dilcharge",";Charge;Events",3,-1.5,1.5) );
  controlHistos.addHistogram( new TH1D("mtsum",";M_{T}(l^{(1)},E_{T}^{miss})+M_{T}(l^{(2)},E_{T}^{miss});Events",100,0,1000) );
  controlHistos.addHistogram( new TH1F ("leadjet", "; Leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadjet", "; Sub-leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("leadlepton", "; Leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadlepton", "; Sub-leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("met", "; #slash{E}_{T} [GeV/c]; Events / (10 GeV/c)", 40, 0.,400.) );

  TString btagalgos[]={"TCHE","TCHP","SSVHE","JBP"};
  int nbtagalgos=sizeof(btagalgos)/sizeof(TString);
  TString btagwps[]={"L","M","T"};
  int nbtagwps=sizeof(btagwps)/sizeof(TString);
  float minba[]={-5,-5,0,0};
  float maxba[]={25,25,6,4};
  for(int iba=0; iba<nbtagalgos; iba++)
    {
      controlHistos.addHistogram( new TH1D(btagalgos[iba],";"+btagalgos[iba]+";Events",100,minba[iba],maxba[iba]) );
      for(int ibwp=0; ibwp<nbtagwps; ibwp++)
	{
	  TString key=btagalgos[iba]+btagwps[ibwp];
	  if(btagCuts.find(key)==btagCuts.end()) continue;
	  TH1F *bmultH=new TH1F (key, ";b-tag multiplicity;Events", 6, 0.,6.) ;
	  for(int ibin=1; ibin<=bmultH->GetXaxis()->GetNbins(); ibin++)
	    {
	      TString label(""); label += ibin-1;
	      if(ibin==bmultH->GetXaxis()->GetNbins()) label ="#geq" + label;
	      bmultH->GetXaxis()->SetBinLabel(ibin,label + " btags");
	    }
	  controlHistos.addHistogram(bmultH);
	}
    }
 
  TString labels[]={"start","#geq 2 jets","=0 b-tags", "=1 b-tags","=2 b-tags"};
  int nsteps=sizeof(labels)/sizeof(TString);
  TH1D *cutflowH=new TH1D("evtflow",";Cutflow;Events",nsteps,0,nsteps);
  for(int ibin=0; ibin<nsteps; ibin++) cutflowH->GetXaxis()->SetBinLabel(ibin+1,labels[ibin]);
  controlHistos.addHistogram( cutflowH );
  TString cats[]={""};//,"jer","jesdown","jesup","btagcen","btagup","btagdown"};
  int nvarcats=sizeof(cats)/sizeof(TString);
  for(int icat=0;icat<nvarcats; icat++)  controlHistos.addHistogram( new TH1D("evtflow"+cats[icat],";Cutflow;Events",nsteps,0,nsteps) );

  controlHistos.initMonitorForStep("ee");
  controlHistos.initMonitorForStep("emu");
  controlHistos.initMonitorForStep("mumu");
  //   controlHistos.initMonitorForStep("etau");
  //   controlHistos.initMonitorForStep("mutau");
  
  ///
  // process events file
  //
  DuplicatesChecker duplicatesChecker;
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

  //total entries to process
  const Int_t totalEntries=evTree->GetEntriesFast();
  if(evEnd<0 || evEnd>totalEntries) evEnd=totalEntries;
  if(evStart > evEnd || totalEntries==0)
    {
      evfile->Close();
      return -1;
    }
  
  //check normalization from first bin
  float cnorm=1.0;
  if(isMC)
    {
      TH1F *cutflowH = (TH1F *) evfile->Get("evAnalyzer/top/cutflow");
      if(cutflowH) cnorm=cutflowH->GetBinContent(1);
    }
  cout << " Xsec x Br=" << xsec << " analyzing " << totalEntries << "/" << cnorm << " original events"<< endl;

  //check PU normalized entries                                                                                                                                                                                                                                                                              
  evTree->Draw(">>elist","normWeight==1");
  TEventList *elist = (TEventList*)gDirectory->Get("elist");
  const Int_t normEntries = elist->GetN();
  if(normEntries==0) cout << "[Warning] Normalized PU entries is 0, check if the PU normalization producer was properly run" << endl;

  //prepare to save summaries
  EventSummaryHandler *spyEvents=0;
  TFile *spyFile=0;
  TDirectory *spyDir=0;
  float summaryWeight(1);
  if(saveSummaryTree && normEntries>0)
    {
      summaryWeight = xsec * float(totalEntries) / (cnorm * float(normEntries) );
      spyEvents = new EventSummaryHandler;
      spyFile = TFile::Open("EventSummaries.root","UPDATE");
      TString evtag=gSystem->BaseName(evurl);
      evtag.ReplaceAll(".root","");
      spyFile->rmdir(evtag);
      spyDir = spyFile->mkdir(evtag);
      TTree *outT = evTree->CloneTree(0);
      outT->SetDirectory(spyDir);
      spyEvents->initTree(outT,false);
    }
  
  //
  // analyze (puf...)
  //
  float selEvents(0);
  int NumberOfDuplicated(0);
  for (int inum=evStart; inum < evEnd; ++inum)
    {
      if(inum%500==0) printf("\r [ %d/100 ] %s",int(100*float(inum-evStart)/float(evEnd)),evurl.Data());

      evTree->GetEvent(inum);
      EventSummary_t &ev = evSummaryHandler.getEvent();
      if(isMC && mcTruthMode>0)
	{
	  if(mcTruthMode==1 && !ev.isSignal) continue;
	  if(mcTruthMode==2 && ev.isSignal)  continue;
	}
      if(duplicatesChecker.isDuplicate(ev.run,ev.lumi, ev.event)){ NumberOfDuplicated++; continue;   }
      float weight = ev.weight/cnorm;    
      float normWeight = ev.normWeight;
      
      bool isSameFlavor(false);
      TString ch("");
      if(ev.cat==dilepton::MUMU)  { isSameFlavor=true; ch="mumu"; }
      if(ev.cat==dilepton::EE)    { isSameFlavor=true; ch="ee"; }
      if(ev.cat==dilepton::EMU)   ch="emu";
      if(ev.cat==dilepton::ETAU)  ch="etau";
      if(ev.cat==dilepton::MUTAU) ch="mutau";
      TString catsToFill[]={"all",ch};
      
      //get particles from event
      PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      
      //reinforce Z-mass window veto for same flavor
      bool isZcand(isSameFlavor && fabs(phys.dil.mass()-91)<15);
      float dilcharge=(phys.leptons[0].id/fabs(phys.leptons[0].id)) *(phys.leptons[1].id/fabs(phys.leptons[1].id));

      //b-tag analysis
      bcomp.compute(phys.nbjets,phys.nljets);
      std::vector<btag::Weight_t> wgt = bcomp.getWeights();
      double p0btags = wgt[0].first;  double p0btags_err=wgt[0].second;
      double p1btags = wgt[1].first;  double p1btags_err=wgt[1].second;
      double p2btags = wgt[2].first;  double p2btags_err=wgt[2].second;
      std::map<TString,int> nbtags;
      for(std::map<TString,float>::iterator it = btagCuts.begin(); it!= btagCuts.end(); it++) nbtags[it->first]=0;
      LorentzVectorCollection jets;
      for(size_t ijet=0; ijet<phys.jets.size(); ijet++)
	{
	  jets.push_back(phys.jets[ijet]);
	  if(!isZcand)
	    {
	      nbtags["TCHEL"]  += (phys.jets[ijet].btag1>btagCuts["TCHEL"]);
	      nbtags["TCHEM"]  += (phys.jets[ijet].btag2>btagCuts["TCHEM"]);
	      nbtags["TCHPT"]  += (phys.jets[ijet].btag2>btagCuts["TCHPT"]);
	      nbtags["SSVHEM"] += (phys.jets[ijet].btag3>btagCuts["SSVHEM"]);		
	      nbtags["JBPL"]   += (phys.jets[ijet].btag4>btagCuts["JBPL"]);
	      nbtags["JBPM"]   += (phys.jets[ijet].btag4>btagCuts["JBPM"]);
	      nbtags["JBPT"]   += (phys.jets[ijet].btag4>btagCuts["JBPT"]);
	      for(size_t ictf=0; ictf<2; ictf++)
		{
		  controlHistos.fillHisto("TCHE",catsToFill[ictf],phys.jets[ijet].btag1,weight);
		  controlHistos.fillHisto("TCHP",catsToFill[ictf],phys.jets[ijet].btag2,weight);
		  controlHistos.fillHisto("SSVHE",catsToFill[ictf],phys.jets[ijet].btag3,weight);
		  controlHistos.fillHisto("JBP",catsToFill[ictf],phys.jets[ijet].btag4,weight);
		}
	    }
	}
      
      //jet variations
      jcomp.compute(jets,phys.met);
      LorentzVectorCollection jerVariedJets     = jcomp.getVariedJets(jet::UncertaintyComputer::JER);
      LorentzVector jerVariedMet                = jcomp.getVariedMet(jet::UncertaintyComputer::JER);
      LorentzVectorCollection jesupVariedJets   = jcomp.getVariedJets(jet::UncertaintyComputer::JES_UP);
      LorentzVector jesUpVariedMet              = jcomp.getVariedMet(jet::UncertaintyComputer::JES_UP);
      LorentzVectorCollection jesdownVariedJets = jcomp.getVariedJets(jet::UncertaintyComputer::JES_DOWN);
      LorentzVector jesDownVariedMet            = jcomp.getVariedMet(jet::UncertaintyComputer::JES_DOWN);    

      //the cutflow
      for(int ivar=0;ivar<(isMC ? nvarcats : 1); ivar++) 
	{
	  for(size_t ictf=0; ictf<2; ictf++)
	    {
	      if(!isZcand) controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],0,weight);

	      //jet energy variations
	      LorentzVectorCollection    jetColl=jets;
	      LorentzVector              theMET=phys.met;     
	      if(cats[ivar]=="jer")      { jetColl=jerVariedJets;      theMET=jerVariedMet;     }
	      if(cats[ivar]=="jesdown")  { jetColl=jesdownVariedJets;  theMET=jesUpVariedMet;   }
	      if(cats[ivar]=="jesup" )   { jetColl=jesupVariedJets;    theMET=jesDownVariedMet; }
	    
	      int nseljets(0);
	      for(size_t ijet=0; ijet<jetColl.size(); ijet++) nseljets += (jetColl[ijet].pt()>30);
	      if(nseljets<2) continue;
	      if(!isZcand) controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],1,weight);

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
		  if(!isZcand)
		    {
		      controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],2,weight*p0weight);
		      controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],3,weight*p1weight);
		      controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],4,weight*p2weight);
		    }
		}
	      //standard b-tag
	      else
		{
		  int btagbin(nbtags["TCHEL"]);
		  if(btagbin>2) btagbin=2;
		  if(!isZcand) controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],2+btagbin,weight);
		}
	      
	      double mtsum=mtComp.compute(phys.leptons[0],theMET) + mtComp.compute(phys.leptons[1],theMET);
	    
	      //standard control histos (no variations)
	      if(ivar==0)
		{
		  if(ictf==0) selEvents+=weight;
		  controlHistos.fillHisto("dilmass",catsToFill[ictf],phys.dil.mass(),weight);
		  if(!isZcand)
		    {
		      controlHistos.fillHisto("mtsum",catsToFill[ictf],mtsum,weight);
		      controlHistos.fillHisto("leadjet",catsToFill[ictf],max(jetColl[0].pt(),jetColl[0].pt()),weight);
		      controlHistos.fillHisto("subleadjet",catsToFill[ictf],min(jetColl[0].pt(),jetColl[1].pt()),weight);
		      controlHistos.fillHisto("leadlepton",catsToFill[ictf],max(phys.leptons[0].pt(),phys.leptons[1].pt()),weight);
		      controlHistos.fillHisto("subleadlepton",catsToFill[ictf],min(phys.leptons[0].pt(),phys.leptons[1].pt()),weight);
		      controlHistos.fillHisto("met",catsToFill[ictf],theMET.pt(),weight);
		      controlHistos.fillHisto("dilmass",catsToFill[ictf],phys.dil.mass(),weight);
		      controlHistos.fillHisto("dilcharge",catsToFill[ictf],dilcharge,weight);
		      for(std::map<TString,int>::iterator it = nbtags.begin(); it!= nbtags.end(); it++) controlHistos.fillHisto(it->first,catsToFill[ictf],it->second,weight);
		      
		      for(size_t ijet=0; ijet<jetColl.size(); ijet++)
			{
			  //get the lepton-jet pairs
			  LorentzVector lj1=phys.leptons[0]+jetColl[ijet];
			  float mlj1=lj1.mass();
			  LorentzVector lj2=phys.leptons[1]+jetColl[ijet];
			  float mlj2=lj2.mass();
			  for(int ictf=0; ictf<2; ictf++)
			    {
			      controlHistos.fillHisto("mlj",catsToFill[ictf],mlj1,weight,true);
			      controlHistos.fillHisto("mlj",catsToFill[ictf],mlj2,weight,true);
			    }
			}
		      		
		      //save summary if required
		      if(spyEvents && normWeight==1)
			{
			  std::vector<float> measurements;
			  ev.weight=summaryWeight;
			  spyEvents->fillTreeWithEvent( ev, measurements );
			}
		    }
		}
	    }
	}
    }
  
  cout << endl << "Selected " << selEvents << " events and found " << NumberOfDuplicated << " duplicates" << endl;

  //
  // close opened files
  // 
  evfile->Close();
  if(spyEvents)
    {
      spyDir->cd();
      spyEvents->getTree()->Write();
      spyFile->Write();
      spyFile->Close();
    }
  
    
  //
  // save histos to local file
  //
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
          if( !((TClass*)hit->second->IsA())->InheritsFrom("TH2") && !((TClass*)hit->second->IsA())->InheritsFrom("TGraph") )
            fixExtremities(hit->second,true,true);
	  hit->second->Write();
	}
    }
  file->Close(); 


  //that's all folks!
}  
