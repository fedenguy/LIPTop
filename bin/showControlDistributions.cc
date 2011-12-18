
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "LIP/Top/interface/EventSummaryHandler.h"

#include "CMGTools/HtoZZ2l2nu/interface/BtagUncertaintyComputer.h"
#include "CMGTools/HtoZZ2l2nu/interface/MacroUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/METUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCanvas.h"
#include "TString.h"
#include "TDirectory.h"
#include "TEventList.h"
#include "TRandom.h"

using namespace std;
using namespace top;


bool sortJetsByTCHE(top::PhysicsObject_Jet a,top::PhysicsObject_Jet b) {   return (a.btag1>b.btag1); }
bool sortJetsByCSV(top::PhysicsObject_Jet a,top::PhysicsObject_Jet b) {   return (a.btag7>b.btag7); }
bool sortJetsByJBP(top::PhysicsObject_Jet a,top::PhysicsObject_Jet b) {   return (a.btag4>b.btag4); }
bool sortJetsByPt(top::PhysicsObject_Jet a,top::PhysicsObject_Jet b) {   return (a.pt()>b.pt()); }

std::pair<float,float> getArcCos(LorentzVector &a, LorentzVector &b)
{
  TVector3 mom1(a.px(),a.py(),a.pz());
  TVector3 mom2(b.px(),b.py(),b.pz());
  double cosine = mom1.Dot(mom2)/(mom1.Mag()*mom2.Mag());
  double arcCosine = acos(cosine);
  return std::pair<float,float>(cosine,arcCosine);
}

std::map<TString, int> getJetSortCounters(top::PhysicsObjectJetCollection &jets, LorentzVector met, LorentzVector ttbar_t)
{
  std::map<TString,int> ctrs;

  ctrs["alpha"]          =  0;
  
  //this is the standard
  top::PhysicsObjectJetCollection csv  = jets; sort(csv.begin(),csv.end(),sortJetsByCSV);

  top::PhysicsObjectJetCollection prunedJets;
  for(size_t i=0; i<csv.size(); i++)
    {
      if(csv[i].pt()<30 || fabs(csv[i].eta())>2.5) continue;
      prunedJets.push_back(jets[i]);
      bool isB(fabs(csv[i].genid)==5);
      ctrs["alpha"] += isB;
    }
  if(prunedJets.size()<2) return ctrs;

  ctrs["1_met30"]     = (met.pt()>30      ? ctrs["alpha"] : -1);
  ctrs["2_met40"]     = (met.pt()>40      ? ctrs["alpha"] : -1);
  ctrs["3_met50"]     = (met.pt()>50      ? ctrs["alpha"] : -1);
  ctrs["4_pttbar25"]  = (ttbar_t.pt()>25  ? ctrs["alpha"] : -1);
  ctrs["5_pttbar50"]  = (ttbar_t.pt()>50  ? ctrs["alpha"] : -1);
  ctrs["6_pttbar100"] = (ttbar_t.pt()>100 ? ctrs["alpha"] : -1);
  ctrs["7_jetpt35"]     = (csv[0].pt()>35 && csv[1].pt()>35) ? ctrs["alpha"] : -1;
  ctrs["8_jetpt40"]     = (csv[0].pt()>40 && csv[1].pt()>40) ? ctrs["alpha"] : -1;

  ctrs["9_sortbycsv"] = (ctrs["alpha"]>=2 && fabs(csv[0].genid)==5 && fabs(csv[1].genid)==5) ? 1 : -1;

  top::PhysicsObjectJetCollection tche = prunedJets; sort(tche.begin(),tche.end(),sortJetsByTCHE);
  ctrs["10_sortbytche"] = (ctrs["alpha"]>=2 && fabs(tche[0].genid)==5 && fabs(tche[1].genid)==5) ? 1 : -1;

  top::PhysicsObjectJetCollection jbp  = prunedJets; sort(jbp.begin(),jbp.end(),sortJetsByJBP);
  ctrs["11_sortbyjbp"] = (ctrs["alpha"]>=2 && fabs(jbp[0].genid)==5 && fabs(jbp[1].genid)==5) ? 1 : -1;
  
  top::PhysicsObjectJetCollection pt   = prunedJets; sort(pt.begin(),pt.end(),sortJetsByPt);
  ctrs["12_sortbypt"] = (ctrs["alpha"]>=2 && fabs(pt[0].genid)==5 && fabs(pt[1].genid)==5) ? 1 : -1;

  return ctrs;
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
  bool runSystematics = runProcess.getParameter<bool>("runSystematics");
  TString etaFileName = runProcess.getParameter<std::string>("etaResolFileName"); gSystem->ExpandPathName(etaFileName);
  TString phiFileName = runProcess.getParameter<std::string>("phiResolFileName"); gSystem->ExpandPathName(phiFileName);
  TString ptFileName  = runProcess.getParameter<std::string>("ptResolFileName");  gSystem->ExpandPathName(ptFileName);
  TString uncFile =  runProcess.getParameter<std::string>("jesUncFileName");      gSystem->ExpandPathName(uncFile);

  //pileup reweighter
  TString proctag=gSystem->BaseName(evurl); 
  //proctag=proctag.ReplaceAll(".root","");
  Ssiz_t pos=proctag.Index(".root");
  proctag.Remove(pos,proctag.Length());
  cout << proctag << " @ " << evurl << endl;
 
  edm::LumiReWeighting *LumiWeights=0;
  double maxPuWeight(-1);
  if(isMC)
    {
      LumiWeights = new edm::LumiReWeighting(runProcess.getParameter<std::string>("mcpileup"), 
					     runProcess.getParameter<std::string>("datapileup"), 
					     proctag.Data(),"pileup");
      for(int n0=0; n0<=30; n0++) maxPuWeight = max( LumiWeights->weight(n0) , maxPuWeight);
      cout << "Max. weight: " << maxPuWeight << endl;
    }

  reweight::PoissonMeanShifter PShiftUp(+0.6);
  reweight::PoissonMeanShifter PShiftDown(-0.6);

  //
  // start auxiliary computers
  //
  btag::UncertaintyComputer bcomp(0.837, 0.95, 0.06, 0.286, 1.11, 0.11);
  JetResolution stdEtaResol(etaFileName.Data(),false);
  JetResolution stdPhiResol(phiFileName.Data(),false);
  JetResolution stdPtResol(ptFileName.Data(),true);
  JetCorrectionUncertainty jecUnc(uncFile.Data());
    
  //
  // control histograms
  //
  SelectionMonitor controlHistos;

  //vertex multiplicity
  controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 25, 0.,25.) );
  controlHistos.addHistogram( new TH1F ("nverticesafteros", "; Vertex multiplicity; Events", 25, 0.,25.) );
  
  ///lepton control
  controlHistos.addHistogram( new TH1D("dilmass",";M(l,l') [GeV/c^{2}];Events",100,0,250) );
  TH1D *lepMult=new TH1D("nleptons",";Leptons;Events",3,0,3);
  lepMult->GetXaxis()->SetBinLabel(1,"=2 leptons");
  lepMult->GetXaxis()->SetBinLabel(2,"=3 leptons");
  lepMult->GetXaxis()->SetBinLabel(3,"#geq 4 leptons");
  controlHistos.addHistogram( lepMult );
  TH1 *sslepMult = (TH1 *) lepMult->Clone("ssnleptons");
  controlHistos.addHistogram( sslepMult );
  controlHistos.addHistogram( new TH1F ("leadlepton", "; Leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadlepton", "; Sub-leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1D("dilcharge",";Charge;Events",3,-1.5,1.5) );
  controlHistos.addHistogram( new TH1D("dphill",";#Delta#phi(l^{(1)},l^{(2)});Events",100,-3.2,3.2) );
  controlHistos.addHistogram( new TH1D("drll",";#Delta R(l^{(1)},l^{(2)});Events",100,0,6) );

  //jet control
  controlHistos.addHistogram( new TH1F ("jet", "; Jet p_{T} [GeV/c]; Events / (10 GeV/c)", 50, 0.,500.) );
  controlHistos.addHistogram( new TH1F ("jetafter1btag", "; Jet p_{T} [GeV/c]; Events / (10 GeV/c)", 40, 0.,400.) );
  controlHistos.addHistogram( new TH1F ("leadjet", "; Leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadjet", "; Sub-leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("mjj", "; M(lead jet,sub-lead jet) [GeV/c^{2}]; Events / (25 GeV/c^{2})", 10, 0.,250.) );
  controlHistos.addHistogram( new TH1F("btagdownvar",";b-tag multiplicity;Events", 3, 0.,3.) );
  controlHistos.addHistogram( new TH1F("btagcenvar",";b-tag multiplicity;Events", 3, 0.,3.) );
  controlHistos.addHistogram( new TH1F("btagupvar",";b-tag multiplicity;Events", 3, 0.,3.) );

  controlHistos.addHistogram( new TH1F("nbtags",";b-tag multiplicity (CSVL);Events", 4, 0.,4.) );
  controlHistos.addHistogram( new TH1F("csvb",";b-tag discriminator (CSV);Events", 100, -0.2,1.2) );
  controlHistos.addHistogram( new TH1F("csvlight",";b-tag discriminator (CSV);Events", 100, -0.2,1.2) );

  //final control
  controlHistos.addHistogram( new TH1F ("pttbar", "; p_{T} (t#bar{t}) [GeV/c]; Events / (10 GeV/c)", 20, 0.,200.) );
  double massAxis[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,400,500,1000,2000};
  int nMassBins=sizeof(massAxis)/sizeof(double)-1;
  controlHistos.addHistogram( new TH1D("mlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV/c^{2}];Lepton-jet pairs",nMassBins,massAxis) );
  controlHistos.addHistogram( new TH1D("drlj","Lepton-jet spectrum;#Delta R(l,j);Lepton-jet pairs",50,0,6) );
  controlHistos.addHistogram( new TH1D("mindrlj","Lepton-jet spectrum;min #Delta R(l,j);Events",50,0,6) );

  controlHistos.addHistogram( new TH2D("bjetsort",";N b-jets reconstructed;Event selection",3,0,3,15,0,15) );

  //event selection histogram
  enum SelectionSteps { SEL2LEPTONS, SELDILEPTON, SELJETS, SELMET, SELOS, SEL0BTAGS, SEL1BTAGS, SEL2BTAGS };
  TString labels[]={"2 leptons",
		    "M>12 #wedge |M-M_{Z}|>15",
		    "#geq 2 jets",
		    "E_{T}^{miss}>30/20 (ee,mumu/emu)",
		    "op. sign",
		    "=0 b-tags",
		    "=1 b-tag",
		    "#geq 2 b-tags"
  };
  int nsteps=sizeof(labels)/sizeof(TString);
  TString cats[]={"","jer","jesdown","jesup","pudown","puup"};
  int nvarcats=sizeof(cats)/sizeof(TString); 
  if(!runSystematics) { nvarcats=1; }
  else                { cout << "Running sytematics: this will take a while..." << endl; }

  //control histograms per variation
  for(int ivar=0;ivar<nvarcats; ivar++) 
    {
      TH1D *cutflowH=new TH1D("evtflow"+cats[ivar],";Cutflow;Events",nsteps,0,nsteps);
      for(int ibin=0; ibin<nsteps; ibin++) cutflowH->GetXaxis()->SetBinLabel(ibin+1,labels[ibin]);
      controlHistos.addHistogram( cutflowH );
      controlHistos.addHistogram( new TH1D("dilmassctr"+cats[ivar],";Region;Events",2,0,2) );
      controlHistos.addHistogram( new TH1D("mtsum"+cats[ivar],";M_{T}(l^{(1)},E_{T}^{miss})+M_{T}(l^{(2)},E_{T}^{miss});Events",100,0,1000) );
      controlHistos.addHistogram( new TH1D("ptsum"+cats[ivar],";p_{T}(l^{(1)})+p_{T}(l^{(2)});Events",100,0,500) );

      TH1F * h = new TH1F ("njets"+cats[ivar], "; Jet multiplicity; Events", 4, 0.,4.);
      h->GetXaxis()->SetBinLabel(1,"=0 jets");
      h->GetXaxis()->SetBinLabel(2,"=1 jets");
      h->GetXaxis()->SetBinLabel(3,"=2 jets");
      h->GetXaxis()->SetBinLabel(4,"#geq 3 jets");
      controlHistos.addHistogram( h );

      controlHistos.addHistogram( new TH1F ("met"+cats[ivar], "; #slash{E}_{T} [GeV/c]; Events / (10 GeV/c)", 50, 0.,500.) );
      
      TString thetallcats[]={"","eq0jets","eq1jets"};
      for(size_t k=0; k<sizeof(thetallcats)/sizeof(TString); k++)
	{
	  controlHistos.addHistogram( new TH1D(thetallcats[k]+"dilarccosine"+cats[ivar],";arcCos(l,l');Events",50,0,3.2) );
	  controlHistos.addHistogram( new TH1D(thetallcats[k]+"lowmetdilarccosine"+cats[ivar],";arcCos(l,l');Events",50,0,3.2) );
	}
    }
  
  controlHistos.initMonitorForStep("ee");
  controlHistos.initMonitorForStep("emu");
  controlHistos.initMonitorForStep("mumu");
  
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
  

  //get the base cut flow histograms
  std::map<TString, TH1F *>  cutflowhistos;
  cutflowhistos["all"]  = (TH1F *) evfile->Get("evAnalyzer/top/cutflow");
  cutflowhistos["ee"]   = (TH1F *) evfile->Get("evAnalyzer/top/ee/ee_cutflow");
  cutflowhistos["mumu"] = (TH1F *) evfile->Get("evAnalyzer/top/mumu/mumu_cutflow");
  cutflowhistos["emu"]  = (TH1F *) evfile->Get("evAnalyzer/top/emu/emu_cutflow");
  //normalization from first bin of the inclusive cut flow
  float cnorm=1.0;
  if(isMC) cnorm=cutflowhistos["all"]->GetBinContent(1);

  for(std::map<TString, TH1F *>::iterator hit = cutflowhistos.begin(); hit != cutflowhistos.end(); hit++)
    {
      for(int ivar=0;ivar<nvarcats; ivar++) 
	{
	  TH1 *h=controlHistos.getHisto("evtflow"+cats[ivar],hit->first);
	  h->SetBinContent(1,hit->second->GetBinContent(2)/cnorm);
	  h->SetBinError(1,hit->second->GetBinError(2)/cnorm);
	  h->SetBinContent(2,hit->second->GetBinContent(3)/cnorm);
	  h->SetBinError(2,hit->second->GetBinError(3)/cnorm);
	}
    }
  cout << " Xsec x Br=" << xsec << " analyzing " << totalEntries << "/" << cnorm << " original events"<< endl;

  //prepare to save summaries
  EventSummaryHandler *spyEvents=0;
  TFile *spyFile=0;
  TDirectory *spyDir=0;
  float summaryWeight(1);
  if(saveSummaryTree)
    {
      gSystem->Exec("mkdir -p " + outdir);
      TString summaryName = outdir + "/" + proctag + "_summary.root";
      gSystem->ExpandPathName(summaryName);
      summaryWeight = xsec / cnorm;
      spyEvents = new EventSummaryHandler;
      spyFile = TFile::Open(summaryName,"RECREATE");
      TString evtag=proctag;
      spyFile->rmdir(evtag);
      spyDir = spyFile->mkdir(evtag);
      TTree *outT = evTree->CloneTree(0);
      outT->SetDirectory(spyDir);
      spyEvents->initTree(outT,false);
      cout << "Creating event summary file:" << summaryName << endl;
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
      top::EventSummary_t &ev = evSummaryHandler.getEvent();
      if(isMC && mcTruthMode>0)
	{
	  if(mcTruthMode==1 && !ev.isSignal) continue;
	  if(mcTruthMode==2 && ev.isSignal)  continue;
	}
      
      if( duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) )
	{ 
	  NumberOfDuplicated++; continue; 
	}
      
      float puweight=1;
      if(LumiWeights) puweight = LumiWeights->weight( ev.ngenpu );

      //get particles from event
      top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      
      bool isSameFlavor(false);
      TString ch("");
      if(ev.cat==MUMU)  { 
	isSameFlavor=true; ch="mumu"; 
      }
      else if(ev.cat==EE)    { 
	if(!phys.leptons[0].hasElectronId(top::PhysicsObject_Lepton::VBTF90) ||
	   !phys.leptons[1].hasElectronId(top::PhysicsObject_Lepton::VBTF90) ) continue;
	isSameFlavor=true; ch="ee"; 
      }
      else if(ev.cat==EMU)  {
	if(!phys.leptons[0].hasElectronId(top::PhysicsObject_Lepton::VBTF90) &&
	   !phys.leptons[1].hasElectronId(top::PhysicsObject_Lepton::VBTF90) ) continue;
	isSameFlavor=false; ch="emu";
      }

      //order jets to apply variations
      top::PhysicsObjectJetCollection orderedJetColl = phys.jets;
      sort(orderedJetColl.begin(),orderedJetColl.end(),sortJetsByCSV);
      LorentzVectorCollection jets;
      for(size_t ijet=0; ijet<orderedJetColl.size(); ijet++) jets.push_back(orderedJetColl[ijet]);
    
      
      //variations
      std::vector<LorentzVectorCollection> jetsVar;
      LorentzVectorCollection metsVar;
      double p0btags(0),p0btags_err(0),p1btags(0),p1btags_err(0),p2btags(0),p2btags_err(0);
      if(runSystematics) 
	{
	  jet::computeVariation(jets,phys.met,jetsVar,metsVar,&stdPtResol,&stdEtaResol,&stdPhiResol,&jecUnc);
	  
	  bcomp.compute(phys.nbjets,phys.nljets);
	  std::vector<btag::Weight_t> wgt = bcomp.getWeights();
	  p0btags = wgt[0].first;  p0btags_err=wgt[0].second;
	  p1btags = wgt[1].first;  p1btags_err=wgt[1].second;
	  p2btags = wgt[2].first;  p2btags_err=wgt[2].second;
	}
      

      //the cutflow with variations
      for(int ivar=0;ivar<(isMC ? nvarcats : 1); ivar++) 
	{
	  float weight=puweight;
	  if(isMC) weight /= cnorm;

	  //objects to use
	  LorentzVectorCollection jetColl = jets;
	  LorentzVector theMET       = phys.met;     
	  if(cats[ivar]=="jer")      { jetColl=jetsVar[jet::JER];      theMET=metsVar[jet::JER]; }
	  if(cats[ivar]=="jesdown")  { jetColl=jetsVar[jet::JES_DOWN]; theMET=metsVar[jet::JES_DOWN]; }
	  if(cats[ivar]=="jesup" )   { jetColl=jetsVar[jet::JES_UP];   theMET=metsVar[jet::JES_UP]; }
	  if(cats[ivar]=="puup" )    { double TotalWeight_plus  = PShiftUp.ShiftWeight( ev.ngenpu );   weight *= TotalWeight_plus;  }
	  if(cats[ivar]=="pudown" )  { double TotalWeight_minus = PShiftDown.ShiftWeight( ev.ngenpu ); weight *= TotalWeight_minus; }

	  LorentzVector l1             = (phys.leptons[0].pt() > phys.leptons[1].pt() ? phys.leptons[0] : phys.leptons[1]);
	  LorentzVector l2             = (phys.leptons[0].pt() < phys.leptons[1].pt() ? phys.leptons[0] : phys.leptons[1]);	  
	  float dilcharge              = (phys.leptons[0].id/fabs(phys.leptons[0].id)) *(phys.leptons[1].id/fabs(phys.leptons[1].id));
	  LorentzVector dileptonSystem = l1+l2;
	  LorentzVector jetSystem      = (jetColl.size() > 1 ? jetColl[0]+jetColl[1] : LorentzVector(0,0,0,0));
	  LorentzVector visible_t      = jetSystem+dileptonSystem;
	  LorentzVector ttbar_t        = jetSystem+dileptonSystem+theMET;

	  std::map<TString, int> jetSortCtr= getJetSortCounters(phys.jets,phys.met,ttbar_t);

	  std::vector<TString> catsToFill;
	  catsToFill.push_back("all");
	  catsToFill.push_back(ch);
	  const size_t nCatsToFill=catsToFill.size();
	  
	  //kinematics with propagated variations
	  int nseljetsLoose(0),nbtags(0);
	  LorentzVectorCollection prunedJetColl;
	  for(size_t ijet=0; ijet<jetColl.size(); ijet++) 
	    {

	      float pt=jetColl[ijet].pt();
	      if(pt<=30 || fabs(jetColl[ijet].eta())>=2.5) continue;
	      //	      nbtags += (orderedJetColl[ijet].btag1>1.7);
	      nbtags += (orderedJetColl[ijet].btag7>0.244);
	      nseljetsLoose ++;
	      prunedJetColl.push_back(jetColl[ijet]);
	    }
	  jetColl=prunedJetColl;
	  int btagbin(nbtags);
	  if(btagbin>2) btagbin=2;
	  
	  double acosine      = getArcCos(l1,l2).second;
	  double ptsum        = l1.pt()+l2.pt();
	  double drll         = deltaR(l1,l2);
	  double dphill       = deltaPhi(l1.phi(),l2.phi());
	  double leadlepmt    = METUtils::transverseMass(l1,theMET,false);
	  double subleadlepmt = METUtils::transverseMass(l2,theMET,false);
	  double mtsum        = leadlepmt+subleadlepmt;
	  double mt           = METUtils::transverseMass(visible_t,theMET,false);
	  
	  //lepton-jet pairs (inclusive)
	  float mindrlj(99999.);
	  std::vector<float> mljs, drljs;
	  for(size_t ijet=0; ijet<jetColl.size(); ijet++)
	    {
	      LorentzVector lj1=phys.leptons[0]+jetColl[ijet];
	      float mlj1=lj1.mass();
	      float drlj1=deltaR(phys.leptons[0],jetColl[ijet]);
	      mljs.push_back(mlj1);
	      drljs.push_back(drlj1);
	      
	      LorentzVector lj2=phys.leptons[1]+jetColl[ijet];
	      float mlj2=lj2.mass();
	      float drlj2=deltaR(phys.leptons[1],jetColl[ijet]);
	      mljs.push_back(mlj2);
	      drljs.push_back(drlj2);
	      
	      mindrlj=min(mindrlj,min(drlj1,drlj2));
	    }
	  
	  
	  bool isInQuarkoniaRegion( dileptonSystem.mass()<12 );
	  bool isZcand(isSameFlavor && fabs(dileptonSystem.mass()-91)<15);
	  bool isEmuInZRegion(!isSameFlavor  && fabs(dileptonSystem.mass()-91)<15);
	  bool passLooseJets(nseljetsLoose>1);
	  bool passMet( theMET.pt()>30 );
	  bool isOS(dilcharge<0);
	  
	  //fill selection histograms
	  for(size_t ictf=0; ictf<nCatsToFill; ictf++)
	    {
	      TString ctf=catsToFill[ictf];

	      if(ivar==0)
		controlHistos.fillHisto("nvertices",ctf,ev.nvtx,weight,true);

	      //reinforce dilepton selection
	      if(isInQuarkoniaRegion) continue;
	      
	      //control for theta_ll and jet multiplicity
	      if(!isZcand && isOS)   
		{
		  TString jetCat="";
		  if(nseljetsLoose==0) jetCat="eq0jets";
		  if(nseljetsLoose==1) jetCat="eq1jets"; 
		  TString metCat( passMet ? "" : "lowmet");
		  controlHistos.fillHisto(jetCat+metCat+"dilarccosine"+cats[ivar],ctf,acosine,weight);
			      
		  if(passMet)
		    {
		      controlHistos.fillHisto("njets"+cats[ivar],ctf,nseljetsLoose,weight);
		    }
		}
	      
	      //>= 2 jets
	      if(!passLooseJets) continue;
	      if(!isZcand)
		{
		  controlHistos.fillHisto("evtflow"+cats[ivar],ctf,SELJETS,weight);
		  controlHistos.fillHisto("met"+cats[ivar],ctf,theMET.pt(),weight);
		}
	      
	      //MET
	      if(!passMet) continue;
	      if(!isZcand)
		{
		  controlHistos.fillHisto("evtflow"+cats[ivar],ctf,SELMET,weight);
		  if(ivar==0) controlHistos.fillHisto("dilcharge",ctf,dilcharge,weight);
		}
 
	      // OS dilepton
	      if(!isOS) continue;
	      if(ivar==0)  controlHistos.fillHisto("dilmass",ctf,dileptonSystem.mass(),weight);
	      controlHistos.fillHisto("dilmassctr"+cats[ivar],ctf,isZcand||isEmuInZRegion ,weight);
	      if(!isZcand)
		{
		  if(ictf==0 && ivar==0) selEvents+=weight;
		  controlHistos.fillHisto("evtflow"+cats[ivar],ctf,SELOS,weight);
		  controlHistos.fillHisto("evtflow"+cats[ivar],ctf,SEL0BTAGS+btagbin,weight);
		  controlHistos.fillHisto("mtsum"+cats[ivar],ctf,mtsum,weight);
		  controlHistos.fillHisto("ptsum"+cats[ivar],ctf,ptsum,weight);
		  if(ivar==0)
		    {
		      if(isMC)
			{
			  int binCtr(0);
			  for(std::map<TString, int>::iterator it = jetSortCtr.begin(); it != jetSortCtr.end(); it++,binCtr++)
			    {
			      if(it->second<0) continue;
			      controlHistos.fill2DHisto("bjetsort",ctf,it->second,binCtr,weight);
			      controlHistos.getHisto("bjetsort",ctf)->GetYaxis()->SetBinLabel(binCtr+1,it->first);
			    }

			  //uncertainty on b-tag selection
			  controlHistos.fillHisto("btagupvar",ctf,0,weight*(p0btags+p0btags_err));
			  controlHistos.fillHisto("btagupvar",ctf,1,weight*(p1btags+p1btags_err));
			  controlHistos.fillHisto("btagupvar",ctf,2,weight*(p2btags+p2btags_err));
			  controlHistos.fillHisto("btagdownvar",ctf,0,weight*(p0btags-p0btags_err));
			  controlHistos.fillHisto("btagdownvar",ctf,1,weight*(p1btags-p1btags_err));
			  controlHistos.fillHisto("btagdownvar",ctf,2,weight*(p2btags-p2btags_err));
			  controlHistos.fillHisto("btagcenvar",ctf,0,weight*p0btags);
			  controlHistos.fillHisto("btagcenvar",ctf,1,weight*p1btags);
			  controlHistos.fillHisto("btagcenvar",ctf,2,weight*p2btags);

			  for(size_t ijet=0; ijet<jetColl.size(); ijet++) 
			    {
			      float pt=jetColl[ijet].pt();
			      if(pt<=30 || fabs(jetColl[ijet].eta())>=2.5) continue;
			      float btagdisc=orderedJetColl[ijet].btag7;
			      bool isMatchedToB( fabs(orderedJetColl[ijet].flavid)==5 );
			      controlHistos.fillHisto( (isMatchedToB ? "csvb" : "csvlight") ,ctf,btagdisc,weight);
			    }	 
			}
			  
  		      controlHistos.fillHisto("nbtags",ctf,nbtags,weight);
		      controlHistos.fillHisto("nverticesafteros",ctf,ev.nvtx,weight,true);
		      controlHistos.fillHisto("pttbar",ctf,ttbar_t.pt(),weight);
		      controlHistos.fillHisto("mjj",ctf,jetSystem.mass(),weight);
		      controlHistos.fillHisto("mt", ctf,mt,weight);
		      controlHistos.fillHisto("mindrlj",ctf,mindrlj,weight,true);
		      controlHistos.fillHisto("nleptons",ctf,phys.leptons.size()-2,weight);
		      for(size_t ijet=0; ijet<jetColl.size(); ijet++) 
			{
			  float pt=jetColl[ijet].pt();
			  if(pt<=30 || fabs(jetColl[ijet].eta())>=2.5) continue;
			  controlHistos.fillHisto("jet",ctf,pt,weight);
			  if(nbtags>0) controlHistos.fillHisto("jetafter1btag",ctf,pt,weight);
			}
		      controlHistos.fillHisto("leadjet",ctf,max(jetColl[0].pt(),jetColl[0].pt()),weight);
		      controlHistos.fillHisto("subleadjet",ctf,min(jetColl[0].pt(),jetColl[1].pt()),weight);
		      controlHistos.fill2DHisto("leadjetvssubleadjet",ctf,max(jetColl[0].pt(),jetColl[0].pt()),min(jetColl[0].pt(),jetColl[1].pt()),weight);
		      controlHistos.fillHisto("leadlepton",ctf,max(phys.leptons[0].pt(),phys.leptons[1].pt()),weight);
		      controlHistos.fillHisto("subleadlepton",ctf,min(phys.leptons[0].pt(),phys.leptons[1].pt()),weight);
		      controlHistos.fillHisto("drll",ctf,drll,weight);
		      controlHistos.fillHisto("dphill",ctf,dphill,weight);
		      for(size_t ilj=0; ilj<mljs.size(); ilj++)
			{
			  controlHistos.fillHisto("mlj",ctf,mljs[ilj],weight,true);
			  controlHistos.fillHisto("drlj",ctf,drljs[ilj],weight,true);
			}
		    }
		}
	      
	      //save summary if required
	      if(isZcand) continue;
	      if(ictf==0 && ivar==0 && spyEvents)
		{
		  //sample the MC according to the PU in data (generate unweighted sample for pseudo-exp.)
		  float rnd=gRandom->Uniform();
		  if(isMC && maxPuWeight>0 && rnd > puweight/maxPuWeight) { ev.normWeight=0; ev.xsecWeight=summaryWeight;}
		  else if (isMC)                                          { ev.normWeight=1; ev.xsecWeight=summaryWeight;}
		  std::vector<float> measurements;
		  spyEvents->fillTreeWithEvent( ev, measurements );
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
      cout << "Finishing summary tree with " << spyEvents->getTree()->GetEntriesFast() << " events" << endl; 
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
  outUrl += proctag + ".root";
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
