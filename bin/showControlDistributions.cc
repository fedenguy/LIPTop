
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

float getArcCos(LorentzVector &a, LorentzVector &b)
{
  TVector3 mom1(a.px(),a.py(),a.pz());
  TVector3 mom2(b.px(),b.py(),b.pz());
  double cosine = mom1.Dot(mom2)/(mom1.Mag()*mom2.Mag());
  double arcCosine = acos(cosine)-TMath::Pi();
  return arcCosine;
//bool aCosmicMuon(arcCosine>-0.005){ //<----------------------  ACTUAL VALUE
}

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
  btagCuts["SSVHET"]=3.05;
   
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
  controlHistos.addHistogram( new TH1D("drlj","Lepton-jet spectrum;#Delta R(l,j);Lepton-jet pairs",50,0,6) );
  controlHistos.addHistogram( new TH1D("dilmass",";M(l,l');Events",100,0,250) );
  TH1 *lepMult=new TH1D("nleptons",";Leptons;Events",3,0,3);
  lepMult->GetXaxis()->SetBinLabel(1,"=2 leptons");
  lepMult->GetXaxis()->SetBinLabel(2,"=3 leptons");
  lepMult->GetXaxis()->SetBinLabel(3,"#geq 4 leptons");
  controlHistos.addHistogram( lepMult );
  controlHistos.addHistogram( (TH1 *) lepMult->Clone("ssnleptons") );
  controlHistos.addHistogram( new TH1D("dilarccosine",";arcCos(l,l');Events",100,-3.2,3.2) );
  controlHistos.addHistogram( new TH1D("dilcharge",";Charge;Events",3,-1.5,1.5) );
  controlHistos.addHistogram( new TH1D("mtsum",";M_{T}(l^{(1)},E_{T}^{miss})+M_{T}(l^{(2)},E_{T}^{miss});Events",100,0,1000) );
  controlHistos.addHistogram( new TH1D("ptsum",";p_{T}(l^{(1)})+p_{T}(l^{(2)});Events",100,0,250) );
  controlHistos.addHistogram( new TH1D("dphill",";#Delta#phi(l^{(1)},l^{(2)});Events",100,-3.2,3.2) );
  controlHistos.addHistogram( new TH1D("drll",";#Delta R(l^{(1)},l^{(2)});Events",100,0,6) );
  controlHistos.addHistogram( new TH1D("alpha",";#alpha_{i};Events",3,0,3) );
  TH1D * jfH = new TH1D("jetflav",";Jet flavor;Jets",3,0,3);
  jfH->GetXaxis()->SetBinLabel(1,"b");
  jfH->GetXaxis()->SetBinLabel(2,"c");
  jfH->GetXaxis()->SetBinLabel(3,"udsg");
  controlHistos.addHistogram( jfH );
  controlHistos.addHistogram( new TH1F ("leadjet", "; Leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadjet", "; Sub-leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH2F ("leadjetvssubleadjet", "; Leading jet p_{T} [GeV/c]; Sub-leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250., 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("leadlepton", "; Leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadlepton", "; Sub-leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("met", "; #slash{E}_{T} [GeV/c]; Events / (10 GeV/c)", 40, 0.,400.) );
  controlHistos.addHistogram( new TH1F ("pttbar", "; p_{T} (t#bar{t}) [GeV/c]; Events / (10 GeV/c)", 20, 0.,200.) );
  controlHistos.addHistogram( new TH1F ("mjj", "; M(lead jet,sub-lead jet) [GeV/c^{2}]; Events / (25 GeV/c^{2})", 10, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("mt", "; M_{T}(visible,E_{T}^{miss}) [GeV/c^{2}]; Events / (25 GeV/c^{2})", 10, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 25, 0.,25.) );
  TH1F * h = new TH1F ("njets", "; Jet multiplicity; Events", 4, 0.,4.);
  h->GetXaxis()->SetBinLabel(1,"=2 jets");
  h->GetXaxis()->SetBinLabel(2,"=3 jets");
  h->GetXaxis()->SetBinLabel(3,"=4 jets");
  h->GetXaxis()->SetBinLabel(4,"#geq 5 jets");
  controlHistos.addHistogram( h );

  TString btagalgos[]={"TCHE","TCHP","SSVHE","JBP"};
  int nbtagalgos=sizeof(btagalgos)/sizeof(TString);
  TString btagwps[]={"L","M","T"};
  int nbtagwps=sizeof(btagwps)/sizeof(TString);
  float minba[]={-5,-5,0,0};
  float maxba[]={25,25,6,8};
  TString ptbins[]={"p_{T}>30","30<p_{T}<50","50<p_{T}<80","80<p_{T}<120","p_{T}>120"};
  const size_t nptbins=sizeof(ptbins)/sizeof(TString);
  for(int iba=0; iba<nbtagalgos; iba++)
    {
      TH1 *bdisc=new TH1D(btagalgos[iba],";"+btagalgos[iba]+";Events",100,minba[iba],maxba[iba]) ;
      controlHistos.addHistogram( bdisc );
      controlHistos.addHistogram((TH1 *)bdisc->Clone(btagalgos[iba]+"b"));
      controlHistos.addHistogram((TH1 *)bdisc->Clone(btagalgos[iba]+"udscg"));
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
	  controlHistos.addHistogram((TH1 *)bmultH->Clone(key+"lowpu"));
	  controlHistos.addHistogram((TH1 *)bmultH->Clone(key+"highpu"));
	  
	  TH1F *bTagEffH=new TH1F (key+"b", ";Jets with b flavor;Events", 10, 0.,10.) ;
	  for(size_t iptbin=0; iptbin<nptbins; iptbin++)
	    {
	      bTagEffH->GetXaxis()->SetBinLabel(iptbin*2+1,ptbins[iptbin]+" (untagged)");
	      bTagEffH->GetXaxis()->SetBinLabel(iptbin*2+2,ptbins[iptbin]+" (tagged)");
	    }
	  controlHistos.addHistogram(bTagEffH);
	  
	  TH1F *udscgTagEffH=(TH1F *) bTagEffH->Clone(key+"udscg");
	  udscgTagEffH->GetXaxis()->SetTitle("Jets with udscg flavor");
	  controlHistos.addHistogram(udscgTagEffH);
	}
    }
 
  TString labels[]={"E_{T}^{miss}>40,0","#geq 2 jets","op. sign", "=0 b-tags", "=1 b-tags","=2 b-tags"};
  int nsteps=sizeof(labels)/sizeof(TString);
  TH1D *cutflowH=new TH1D("evtflow",";Cutflow;Events",nsteps,0,nsteps);
  for(int ibin=0; ibin<nsteps; ibin++) cutflowH->GetXaxis()->SetBinLabel(ibin+1,labels[ibin]);
  controlHistos.addHistogram( cutflowH );
  TString cats[]={""};//,"jer","jesdown","jesup","btagcen","btagup","btagdown"};
  int nvarcats=sizeof(cats)/sizeof(TString);
  for(int icat=0;icat<nvarcats; icat++)  controlHistos.addHistogram( new TH1D("evtflow"+cats[icat],";Cutflow;Events",nsteps,0,nsteps) );

  TString chargeCats[]={"ss",""};
  controlHistos.initMonitorForStep("ssall");
  for(size_t icc=0; icc<2; icc++)
    {
      controlHistos.initMonitorForStep(chargeCats[icc]+"ee");
      controlHistos.initMonitorForStep(chargeCats[icc]+"emu");
      controlHistos.initMonitorForStep(chargeCats[icc]+"mumu");
      controlHistos.initMonitorForStep(chargeCats[icc]+"eq2jets");
      controlHistos.initMonitorForStep(chargeCats[icc]+"geq3jets");
      //   controlHistos.initMonitorForStep("etau");
      //   controlHistos.initMonitorForStep("mutau");
    }
  
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

      //get particles from event
      PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      
      bool isSameFlavor(false);
      TString ch("");
      if(ev.cat==dilepton::MUMU)  { isSameFlavor=true; ch="mumu"; }
      if(ev.cat==dilepton::EE)    { isSameFlavor=true; ch="ee"; }
      if(ev.cat==dilepton::EMU)   ch="emu";
      if(ev.cat==dilepton::ETAU)  ch="etau";
      if(ev.cat==dilepton::MUTAU) ch="mutau";
      int nseljets(0);
      for(size_t ijet=0; ijet<phys.jets.size(); ijet++) nseljets += (phys.jets[ijet].pt()>30);
      TString jetCat( nseljets==2 ? "eq2jets" : "geq3jets" ) ;
      TString catsToFill[]={"all",ch, jetCat};
      const size_t nCatsToFill=sizeof(catsToFill)/sizeof(TString);            

      //reinforce Z-mass window veto for same flavor
      bool isZcand(isSameFlavor && fabs(phys.dil.mass()-91)<15);
      float dilcharge=(phys.leptons[0].id/fabs(phys.leptons[0].id)) *(phys.leptons[1].id/fabs(phys.leptons[1].id));
      bool isOS(dilcharge<0);
      if(!isOS) { 
	for(size_t ictf=0; ictf<nCatsToFill; ictf++) catsToFill[ictf]="ss"+catsToFill[ictf];
      }

      //b-tag analysis
      bcomp.compute(phys.nbjets,phys.nljets);
      std::vector<btag::Weight_t> wgt = bcomp.getWeights();
      double p0btags = wgt[0].first;  double p0btags_err=wgt[0].second;
      double p1btags = wgt[1].first;  double p1btags_err=wgt[1].second;
      double p2btags = wgt[2].first;  double p2btags_err=wgt[2].second;
      std::map<TString,int> nbtags;
      std::map<TString,int> nFlavBtags[nptbins];
      int nbs[nptbins];
      int nudscg[nptbins];
      for(std::map<TString,float>::iterator it = btagCuts.begin(); it!= btagCuts.end(); it++) 
	{
	  nbtags[ it->first ]=0;
	  for(size_t iptbin=0; iptbin<nptbins; iptbin++) 
	    {
	      nFlavBtags[iptbin][it->first+"b"]=0;
	      nFlavBtags[iptbin][it->first+"udscg"]=0;
	      nbs[iptbin]=0;
	      nudscg[iptbin]=0;
	    }
	}

      LorentzVectorCollection jets;
      int nGoodassigments(0);
      for(size_t iptbin=0; iptbin<nptbins; iptbin++)	{ nbs[iptbin]=0; nudscg[iptbin]=0; }
      for(size_t ijet=0; ijet<phys.jets.size(); ijet++)
	{
	  jets.push_back(phys.jets[ijet]);
	  float pt=phys.jets[ijet].pt();

	  bool hasBflavor(fabs(phys.jets[ijet].flavid)==5);
	  int jetFlavBin(2);
	  if(hasBflavor) jetFlavBin=0;
	  else if(fabs(phys.jets[ijet].flavid)==4) jetFlavBin=1;
	  TString flavCat(( hasBflavor ? "b" : "udscg" ) );
	  int ptbins[2]={0,0};
	  if(pt>30)  ptbins[1]++;
	  if(pt>50)  ptbins[1]++;
	  if(pt>80)  ptbins[1]++;
	  if(pt>120) ptbins[1]++;
	  for(size_t iptbin=0; iptbin<2; iptbin++)
	    {
	      nbs[ptbins[iptbin]]    += hasBflavor;
	      nudscg[ptbins[iptbin]] += !hasBflavor;
	    }

	  int assignCode=(phys.leptons[0].genid*phys.jets[ijet].genid);
	  bool isCorrect( (assignCode<0) && hasBflavor);
	  assignCode=(phys.leptons[1].genid*phys.jets[ijet].genid);
	  isCorrect |= ( (assignCode<0) && hasBflavor);
	  nGoodassigments+=isCorrect;

	  if(!isZcand)
	    {
	      nbtags["TCHEL"]  += (phys.jets[ijet].btag1>btagCuts["TCHEL"]);
	      nbtags["TCHEM"]  += (phys.jets[ijet].btag1>btagCuts["TCHEM"]);
	      nbtags["TCHPT"]  += (phys.jets[ijet].btag2>btagCuts["TCHPT"]);
	      nbtags["SSVHEM"] += (phys.jets[ijet].btag3>btagCuts["SSVHEM"]);		
	      nbtags["SSVHET"] += (phys.jets[ijet].btag3>btagCuts["SSVHET"]);		
	      nbtags["JBPL"]   += (phys.jets[ijet].btag4>btagCuts["JBPL"]);
	      nbtags["JBPM"]   += (phys.jets[ijet].btag4>btagCuts["JBPM"]);
	      nbtags["JBPT"]   += (phys.jets[ijet].btag4>btagCuts["JBPT"]);
	      
	      if(isMC)
		{
		  //b-tag counting per jet flavor and pt-bin
		  for(size_t iptbin=0; iptbin<2; iptbin++)
		    {
		      nFlavBtags[ptbins[iptbin]]["TCHEL"+flavCat]  += (phys.jets[ijet].btag1>btagCuts["TCHEL"]);
		      nFlavBtags[ptbins[iptbin]]["TCHEM"+flavCat]  += (phys.jets[ijet].btag1>btagCuts["TCHEM"]);
		      nFlavBtags[ptbins[iptbin]]["TCHPT"+flavCat]  += (phys.jets[ijet].btag2>btagCuts["TCHPT"]);
		      nFlavBtags[ptbins[iptbin]]["SSVHEM"+flavCat] += (phys.jets[ijet].btag3>btagCuts["SSVHEM"]);		
		      nFlavBtags[ptbins[iptbin]]["SSVHET"+flavCat] += (phys.jets[ijet].btag3>btagCuts["SSVHET"]);		
		      nFlavBtags[ptbins[iptbin]]["JBPL"+flavCat]   += (phys.jets[ijet].btag4>btagCuts["JBPL"]);
		      nFlavBtags[ptbins[iptbin]]["JBPM"+flavCat]   += (phys.jets[ijet].btag4>btagCuts["JBPM"]);
		      nFlavBtags[ptbins[iptbin]]["JBPT"+flavCat]   += (phys.jets[ijet].btag4>btagCuts["JBPT"]);
		    }
		}

	      for(size_t ictf=0; ictf<nCatsToFill; ictf++)
		{
		  controlHistos.fillHisto("jetflav",catsToFill[ictf],jetFlavBin,weight);
		  controlHistos.fillHisto("TCHE",catsToFill[ictf],phys.jets[ijet].btag1,weight);
		  controlHistos.fillHisto("TCHP",catsToFill[ictf],phys.jets[ijet].btag2,weight);
		  controlHistos.fillHisto("SSVHE",catsToFill[ictf],phys.jets[ijet].btag3,weight);
		  controlHistos.fillHisto("JBP",catsToFill[ictf],phys.jets[ijet].btag4,weight);
		  
		  if(!isMC) continue;
		  //per jet flavor
		  controlHistos.fillHisto("TCHE"+flavCat, catsToFill[ictf],phys.jets[ijet].btag1,weight);
		  controlHistos.fillHisto("TCHP"+flavCat,catsToFill[ictf],phys.jets[ijet].btag2,weight);
		  controlHistos.fillHisto("SSVHE"+flavCat,catsToFill[ictf],phys.jets[ijet].btag3,weight);
		  controlHistos.fillHisto("JBP"+flavCat,catsToFill[ictf],phys.jets[ijet].btag4,weight);
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
	  for(size_t ictf=0; ictf<nCatsToFill; ictf++)
	    {
	      //MET cut has been applied
	      if(!isZcand) controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],0,weight);
		
	      //jet energy variations
	      LorentzVectorCollection    jetColl=jets;
	      LorentzVector              theMET=phys.met;     
	      if(cats[ivar]=="jer")      { jetColl=jerVariedJets;      theMET=jerVariedMet;     }
	      if(cats[ivar]=="jesdown")  { jetColl=jesdownVariedJets;  theMET=jesUpVariedMet;   }
	      if(cats[ivar]=="jesup" )   { jetColl=jesupVariedJets;    theMET=jesDownVariedMet; }
	      
	      // check selected jets after variation (if any)
	      int nseljets(0);
	      for(size_t ijet=0; ijet<jetColl.size(); ijet++) nseljets += (jetColl[ijet].pt()>30);
	      if(nseljets<2) continue;
	      double acosine(getArcCos(phys.leptons[0],phys.leptons[1]));
	      controlHistos.fillHisto("dilarccosine",catsToFill[ictf],acosine,weight);
	      if(!isZcand) 
		{
		  controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],1,weight);
		  if(ivar==0) controlHistos.fillHisto("dilcharge",catsToFill[ictf],dilcharge,weight);
		}
		
	      if(!isZcand) controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],2,weight);

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
		      controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],3,weight*p0weight);
		      controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],4,weight*p1weight);
		      controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],5,weight*p2weight);
		    }
		}
	      //standard b-tag
	      else
		{
		  int btagbin(nbtags["TCHEL"]);
		  if(btagbin>2) btagbin=2;
		  if(!isZcand) controlHistos.fillHisto("evtflow"+cats[ivar],catsToFill[ictf],3+btagbin,weight);
		}
	      
	      double mtsum=mtComp.compute(phys.leptons[0],theMET) + mtComp.compute(phys.leptons[1],theMET);
	      double ptsum=phys.leptons[0].pt() + phys.leptons[1].pt();
	      double drll=deltaR(phys.leptons[0],phys.leptons[1]);
	      double dphill=deltaPhi(phys.leptons[0].phi(),phys.leptons[1].phi());
	      LorentzVector dileptonSystem = phys.leptons[0] + phys.leptons[1];
	      LorentzVector jetSystem      = jetColl[0]+jetColl[1];
	      LorentzVector visibleSystem  = jetSystem+dileptonSystem;
	      LorentzVector ttbar_t        = visibleSystem+theMET;
	      double systemMt              = mtComp.compute(ttbar_t,theMET,false);

	      //standard control histos (no variations)
	      if(ivar==0)
		{
		  if(ictf==0) selEvents+=weight;
		  controlHistos.fillHisto("dilmass",catsToFill[ictf],phys.dil.mass(),weight);
		  if(!isZcand)
		    {
		      controlHistos.fillHisto("pttbar",catsToFill[ictf],ttbar_t.pt(),weight);
		      controlHistos.fillHisto("mjj",catsToFill[ictf],jetSystem.mass(),weight);
		      controlHistos.fillHisto("mt", catsToFill[ictf],systemMt,weight);

		      controlHistos.fillHisto("nleptons",catsToFill[ictf],phys.leptons.size()-2,weight);
		      controlHistos.fillHisto("alpha",catsToFill[ictf],min(nGoodassigments,2),weight);
		      controlHistos.fillHisto("mtsum",catsToFill[ictf],mtsum,weight);
		      controlHistos.fillHisto("ptsum",catsToFill[ictf],ptsum,weight);
		      controlHistos.fillHisto("leadjet",catsToFill[ictf],max(jetColl[0].pt(),jetColl[0].pt()),weight);
		      controlHistos.fillHisto("subleadjet",catsToFill[ictf],min(jetColl[0].pt(),jetColl[1].pt()),weight);
		      controlHistos.fill2DHisto("leadjetvssubleadjet",catsToFill[ictf],max(jetColl[0].pt(),jetColl[0].pt()),min(jetColl[0].pt(),jetColl[1].pt()),weight);
		      controlHistos.fillHisto("leadlepton",catsToFill[ictf],max(phys.leptons[0].pt(),phys.leptons[1].pt()),weight);
		      controlHistos.fillHisto("subleadlepton",catsToFill[ictf],min(phys.leptons[0].pt(),phys.leptons[1].pt()),weight);
		      controlHistos.fillHisto("met",catsToFill[ictf],theMET.pt(),weight);
		      controlHistos.fillHisto("dilmass",catsToFill[ictf],phys.dil.mass(),weight);
		      controlHistos.fillHisto("drll",catsToFill[ictf],drll,weight);
		      controlHistos.fillHisto("dphill",catsToFill[ictf],dphill,weight);
		      for(std::map<TString,int>::iterator it = nbtags.begin(); it!= nbtags.end(); it++)
			{
			  int btagCtr=it->second;
			  controlHistos.fillHisto(it->first,catsToFill[ictf],btagCtr,weight);
			  controlHistos.fillHisto(it->first+(ev.nvtx<6 ? "lowpu" :"highpu"),catsToFill[ictf],btagCtr,weight);
			}
		      if(isMC)
			{
			  for(size_t iptbin=0; iptbin<nptbins; iptbin++)	
			    for(std::map<TString,int>::iterator it = nFlavBtags[iptbin].begin(); it!= nFlavBtags[iptbin].end(); it++)
			      {

				int nmatchedflav=nbs[iptbin];
				if(it->first.EndsWith("udscg") )  nmatchedflav=nudscg[iptbin];
				int ntagged=it->second;
				int nnottagged=nmatchedflav-ntagged;
				controlHistos.fillHisto(it->first,catsToFill[ictf],iptbin*2,nnottagged*weight);
				controlHistos.fillHisto(it->first,catsToFill[ictf],iptbin*2+1,ntagged*weight);
			      }
			}
		      
		      for(size_t ijet=0; ijet<jetColl.size(); ijet++)
			{
			  //get the lepton-jet pairs
			  LorentzVector lj1=phys.leptons[0]+jetColl[ijet];
			  float mlj1=lj1.mass();
			  float drlj1=deltaR(phys.leptons[0],jetColl[ijet]);
			  LorentzVector lj2=phys.leptons[1]+jetColl[ijet];
			  float mlj2=lj2.mass();
			  float drlj2=deltaR(phys.leptons[1],jetColl[ijet]);
			  if(!isOS && (mlj1>200 || mlj2>200) )
			    {
			      if(ictf==0) cout << ev.run << ":" << ev.lumi << ":" << ev.event <<  " " << ev.cat << " " << nseljets<< endl;   
			    }
			  for(size_t ictf=0; ictf<nCatsToFill; ictf++)
			    {
			      controlHistos.fillHisto("mlj",catsToFill[ictf],mlj1,weight,true);
			      controlHistos.fillHisto("mlj",catsToFill[ictf],mlj2,weight,true);
			      controlHistos.fillHisto("drlj",catsToFill[ictf],drlj1,weight,true);
			      controlHistos.fillHisto("drlj",catsToFill[ictf],drlj2,weight,true);
			    }
			}
		      
		      controlHistos.fillHisto("njets",catsToFill[ictf],nseljets-2,weight,true);
		      controlHistos.fillHisto("nvertices",catsToFill[ictf],ev.nvtx,weight,true);

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
