#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>
#include <algorithm>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/KinResultsHandler.h"
#include "LIP/Top/interface/HistogramAnalyzer.h"
#include "LIP/Top/interface/MassLikelihoodExtractor.h"
#include "LIP/Top/interface/KinAnalysis.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"
#include "CMGTools/HtoZZ2l2nu/interface/MacroUtils.h"

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

using namespace std;
using namespace top;

//
struct bTagSorter{
  bool operator() (PhysicsObject_Jet a, PhysicsObject_Jet b) {   return (a.btag7>b.btag7); }
} sortByBtag;

//
struct ptSorter{
  bool operator() (PhysicsObject_Lepton a, PhysicsObject_Lepton b) {   return (a.pt()>b.pt()); }
} sortByPt;


//
int getPreferredCombination(TH1F *h1,TH1F *h2, int minCounts=100)
{
  int prefComb(-1);
  if(h1==0 || h2==0) return prefComb;

  //assign the preferred combination 
  // 1 - combinations with larger number of solutions are preferred
  // 2 - in case of ambiguity (number of solutions differ by <10%) take the combination with the more pronounced peak
  // 3 - reject combination if number of solutions < minCounts
  double prefMax[2]    = {h1->GetMaximum(), h2->GetMaximum() };
  double prefCounts[2] = {h1->Integral(),   h2->Integral() };
  if(prefCounts[0]>prefCounts[1]) 
    {
      prefComb=1;
      if(prefCounts[0]<minCounts) prefComb=-1;
      else if(prefCounts[1]/prefCounts[0]>0.90 && prefMax[1]>prefMax[0] && prefMax[1]>minCounts)  prefComb=2;
    }           
  else if(prefCounts[1]>minCounts)
    {
      prefComb=2;
      if(prefCounts[0]/prefCounts[1]>0.90 && prefMax[0]>prefMax[1] && prefMax[0]>minCounts) prefComb=1;
    }
  
  //all done
  return prefComb;
}



//
int main(int argc, char* argv[])
{
  SelectionMonitor controlHistos; //plot storage
  HistogramAnalyzer histoAnalyzer; 

  DuplicatesChecker duplicatesChecker;

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
  TString evurl               = runProcess.getParameter<std::string>("input");
  TString outUrl              = runProcess.getParameter<std::string>("outdir");
  TString kindir              = runProcess.getParameter<std::string>("kindir");
  bool isMC                   = runProcess.getParameter<bool>("isMC");
  int mcTruthMode             = runProcess.getParameter<int>("mctruthmode");
  TString dirname             = runProcess.getParameter<std::string>("dirName");
  bool saveSummaryTree        = runProcess.getParameter<bool>("saveSummaryTree");
  double xsec                 = runProcess.getParameter<double>("xsec");

  //pileup reweighter                                                        
  TString proctag=gSystem->BaseName(evurl); proctag=proctag.ReplaceAll(".root","");
  edm::LumiReWeighting *LumiWeights=0;
  double maxPuWeight(-1);
  if(isMC) 
    {
      LumiWeights = new edm::LumiReWeighting(runProcess.getParameter<std::string>("mcpileup"),
					     runProcess.getParameter<std::string>("datapileup"),
					     //      proctag.Data(),"pileup");
					     "pileup","pileup");
      for(int n0=0; n0<=30; n0++) maxPuWeight = max( LumiWeights->weight(n0) , maxPuWeight);
      cout << "Max. weight: " << maxPuWeight << endl;
    }

  //book histos
  controlHistos.addHistogram( new TH1F ("njets", ";Jets;Events", 6, 0.,6.) );
  controlHistos.addHistogram( new TH1F ("btags", ";b-tag multiplicity (CSVL);Events", 6, 0.,6.) );
  controlHistos.addHistogram( new TH1F ("btagsraw", ";b-tag multiplicity (CSVL);Events", 6, 0.,6.) );
  for(int ibin=1; ibin<=controlHistos.getHisto("njets","all")->GetXaxis()->GetNbins(); ibin++)
    {
      TString label(""); label += ibin-1; 
      if(controlHistos.getHisto("njets","all")->GetXaxis()->GetNbins()==ibin) label ="#geq" + label;
      controlHistos.getHisto("njets","all")->GetXaxis()->SetBinLabel(ibin,label + " jets");
      controlHistos.getHisto("btags","all")->GetXaxis()->SetBinLabel(ibin,label + " b-tags");
      controlHistos.getHisto("btagsraw","all")->GetXaxis()->SetBinLabel(ibin,label + " b-tags");
    }

  controlHistos.addHistogram( new TH1F("totalMassLikelihood","M_{t} likelihood from KINb solutions;M_{t} [GeV/c^{2}];Likelihood / (2.5 GeV/c^{2})",800,0,2000) );
  TH1F * h=new TH1F ("evtflow", "; Cutflow; Events", 5, 0.,5.);
  h->GetXaxis()->SetBinLabel(1,"E_{T}^{miss}>30");
  h->GetXaxis()->SetBinLabel(2,"KIN");
  h->GetXaxis()->SetBinLabel(3,"=0 b-tags");
  h->GetXaxis()->SetBinLabel(4,"=1 b-tags");
  h->GetXaxis()->SetBinLabel(5,"#geq 2 b-tags");
  controlHistos.addHistogram( h );

  controlHistos.addHistogram( new TH1F ("taupt", "; p_{T} (#tau); Events / (20 GeV/c)", 20, 0.,400.) );
  controlHistos.addHistogram( new TH1F ("taueta", "; #eta (#tau}; Events", 20, 0,2.5) );
  controlHistos.addHistogram( new TH1F ("leadjet", "; Leading jet p_{T} [GeV/c]; Events / (20 GeV/c)", 20, 0.,400.) );
  controlHistos.addHistogram( new TH1F ("subleadjet", "; Sub-leading jet p_{T} [GeV/c]; Events / (20 GeV/c)", 20, 0.,400.) );
  controlHistos.addHistogram( new TH1F ("leadlepton", "; Leading lepton p_{T} [GeV/c]; Events / (20 GeV/c)", 20, 0.,400.) );
  controlHistos.addHistogram( new TH1F ("subleadlepton", "; Sub-leading lepton p_{T} [GeV/c]; Events / (20 GeV/c)", 20, 0.,400.) );
  controlHistos.addHistogram( new TH1F ("met", "; #slash{E}_{T} [GeV/c]; Events / (50 GeV/c)", 10, 0.,500.) );
  controlHistos.addHistogram( new TH1F ("ht", "; #sum_{jets} [GeV/c]; Events / (40 GeV/c)",25, 0.,1000.) );
  controlHistos.addHistogram( new TH1F ("st", "; #sum_{leptons,E_{T}^{miss}} p_{T} [GeV/c]; Events / (40 GeV/c)",25, 0.,1000.) );
  controlHistos.addHistogram( new TH1F ("sumpt", "; #sum_{leptons} p_{T} [GeV/c]; Events / (20 GeV/c)",25, 0.,500.) );
  controlHistos.addHistogram( new TH1F ("htlep", "; #sum_{jets,leptons,E_{T}^{miss}} [GeV/c]; Events / (20 GeV/c)",70, 0.,1400.) );
  controlHistos.addHistogram( new TH1F ("mtop", "; m_{Top} [GeV/c^{2}]; Events / (20 GeV/c^{2})", 50, 0.,1000.) );
  controlHistos.addHistogram( new TH1F ("ptttbar", "; p_{t#bar{t}} [GeV/c]; Events / (40 GeV/c)", 25, 0.,1000.) );
  controlHistos.addHistogram( new TH1F ("mttbar", "; Mass(t,#bar{t}) [GeV/c^{2}]; Events / (100 GeV/c^{2})", 30, 0.,3000.) );
  controlHistos.addHistogram( new TH2F ("mtopvsdilmass", "; m_{Top} [GeV/c^{2}]; Mass(l,l') [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,500.) );
  controlHistos.addHistogram( new TH2F ("mtopvsmlj", "; m_{Top} [GeV/c^{2}]; Mass(l,j) [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,500.) );
  controlHistos.addHistogram( new TH2F ("mtopvsmttbar", "; m_{Top} [GeV/c^{2}]; Mass(t,#bar{t}) [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,2000.) );
  controlHistos.addHistogram( new TH2F ("mtopvsnvtx", "; m_{Top} [GeV/c^{2}]; Vertices; Events", 100, 0.,500.,30,0,30) );
  controlHistos.addHistogram( new TH2F ("mtopvsmet", "; m_{Top} [GeV/c^{2}]; E_{T}^{miss} [GeV/c]; Events", 100, 0.,500.,10,0.,500.) );
  controlHistos.addHistogram( new TH2F ("mtopvsptjet", "; m_{Top} [GeV/c^{2}]; p_{T}^{jet}; Events", 100, 0.,500.,4,30.,50.) );
  
  TString cats[]={"ee","mumu","emu"};//"etau","mutau"};

  size_t ncats=sizeof(cats)/sizeof(TString);
  TString subcats[]={"","eq0btags","eq1btags","geq2btags","zcands","ss"};
  size_t nsubcats=sizeof(subcats)/sizeof(TString);
  for(size_t icat=0; icat<ncats; icat++)
    for(size_t jcat=0; jcat<nsubcats; jcat++)
      controlHistos.initMonitorForStep(cats[icat]+subcats[jcat]);

  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  //  evurl.ReplaceAll(".root","_summary.root");
  outUrl += gSystem->BaseName(evurl);
  TFile *file=TFile::Open(outUrl, "recreate");
  
  //fix entries flag
  ofstream *outf=0;
  if(!isMC) outf=new ofstream("highmassevents.txt",ios::app);
  
  //process events file
  gSystem->ExpandPathName(TString(gSystem->BaseName(evurl)));

  TString baseTag = gSystem->BaseName(evurl);
  //  baseTag.ReplaceAll("_summary.root","");
  baseTag.ReplaceAll(".root","");

  TFile *evfile = TFile::Open(evurl);
  //  cout << "evurl: " << evurl << ", baseTag: " << baseTag << endl;
  if(evfile==0) return -1;
  if(evfile->IsZombie()) return -1;
  EventSummaryHandler evSummaryHandler;
  if( !evSummaryHandler.attachToTree( (TTree *)evfile->Get(TString("evAnalyzer/data")) ) ) 
  //if( !evSummaryHandler.attachToTree( (TTree *)evfile->Get(baseTag+TString("/")+dirname) ) ) 
    {
      evfile->Close();
      return -1;
    }  
  TTree *evTree=evSummaryHandler.getTree();

  float cnorm=1;
  if(isMC)
    {
      TH1F *cutflowH = (TH1F *) evfile->Get("evAnalyzer/top/evtflow");
      if(cutflowH==0) cutflowH=(TH1F *) evfile->Get("evAnalyzer/top/cutflow");
      if(cutflowH)    cnorm=cutflowH->GetBinContent(1);
    }



  //init event spy
  EventSummaryHandler *spyEvents=0;
  TFile *spyFile=0;
  TDirectory *spyDir=0;  
  float summaryWeight(1);
  if(saveSummaryTree)
    {
      spyEvents = new EventSummaryHandler;
      spyFile = TFile::Open("EventSummaries.root","UPDATE");
      TString evtag=gSystem->BaseName(evurl);
      //      evtag.ReplaceAll("_summary.root","");
      evtag.ReplaceAll(".root","");
      if(kindir!="std") evtag += "_" + kindir;
      spyFile->rmdir(evtag);
      spyDir = spyFile->mkdir(evtag);
      TTree *outT = new TTree("data","Event summary");
      outT->SetDirectory(spyDir);
      spyEvents->initTree(outT);
      summaryWeight = xsec / cnorm;
    }

  //process kin file
  TString kinUrl(evurl);
  //  kinUrl.ReplaceAll("_summary.root","/"+kindir);
  kinUrl.ReplaceAll(".root","/"+kindir);
  gSystem->ExpandPathName(kinUrl);
  cout << "Kin results from " << kinUrl << " to be processed with summary from " << evurl << endl;
 
  KinResultsHandler kinHandler;
  kinHandler.init(kinUrl,false);
  TChain *t=kinHandler.getResultsChain();

  std::map<TString,int> selEvents;
  for (int inum=0; inum < evTree->GetEntriesFast(); ++inum)
    {
      evTree->GetEvent(inum);
      EventSummary_t &ev = evSummaryHandler.getEvent();
      if(isMC)
	{
	  if(mcTruthMode==1 && !ev.isSignal) continue;
	  if(mcTruthMode==2 && ev.isSignal) continue;
	}
      if( duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) continue;
      top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);

      //preselect for KIN level: 2 jets, MET (no OS, no Z-veto)
      bool isSameFlavor(false);
      if(ev.cat==MUMU || ev.cat==EE) isSameFlavor=true;

      float minJetPt(30);
      if(ev.cat==ETAU) minJetPt=35;
    
      int njets(0),nbjets(0);
      for(size_t ijet=0; ijet<phys.jets.size(); ijet++)
	{
	  if(phys.jets[ijet].pt()<minJetPt || fabs(phys.jets[ijet].eta())>2.5) continue;
	  njets++;
	  nbjets += (phys.jets[ijet].btag1>1.7);
	}
      if(njets<2 || ((ev.cat==ETAU || ev.cat==MUTAU) && nbjets==0) ) continue;
    
      bool passMet( phys.met.pt()>30 );//(!isSameFlavor && phys.met.pt() > 30) || ( isSameFlavor && phys.met.pt()>30) );
      if(!passMet) continue;
   

      TString key(""); key+= ev.run; key+="-"; key += ev.lumi; key+="-"; key += ev.event;
      selEvents[key]=inum;
    }
  
  cout << "Selected : " << selEvents.size() << " events,  looking for kin results" << endl;

  //loop over kin results

  MassLikelihoodExtractor* massLikelihood = new MassLikelihoodExtractor(false);
  

  int nresults(0),neventsused(0);
  std::set<TString> usedEvents;
  for (int inum=0; inum < t->GetEntries(); ++inum)
    {
      t->GetEvent(inum);
      
      //get original event
      Int_t irun,ievent,ilumi;
      kinHandler.getEventInfo(irun,ievent,ilumi);
      
      TString key("");  key+= irun; key+="-"; key += ilumi;  key+="-"; key += ievent;
      if(selEvents.find(key)==selEvents.end())  continue;
      if(usedEvents.find(key)!=usedEvents.end()) continue;
      usedEvents.insert(key);

      nresults++;
      evTree->GetEntry( selEvents[key] );

      //get event summary
      EventSummary_t &ev = evSummaryHandler.getEvent();
   
      std::vector<TString> categs;
      categs.push_back("all");
      if(ev.cat==MUMU)  categs.push_back("mumu");
      if(ev.cat==EE)  categs.push_back("ee");
      if(ev.cat==EMU) categs.push_back("emu");
      if(ev.cat==ETAU) categs.push_back("etau");
      if(ev.cat==MUTAU) categs.push_back("mutau");


      //fill histos
      float weight = ev.weight;
      int normWeight=0;
      if(LumiWeights && ev.cat!=ETAU && ev.cat!=MUTAU) 
	{
	  weight = LumiWeights->weight( ev.ngenpu );
	  float rnd=gRandom->Uniform();
	  if(rnd > weight/maxPuWeight) normWeight=1;
	}
      //else if (isMC)  weight = ev.weight;


      PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      top::PhysicsObjectJetCollection prunedJets;
      int nbtags(0),nbtagscor(0);
      for(size_t ijet=0; ijet<phys.jets.size(); ijet++)
	{
	  if(phys.jets[ijet].pt()<30 || fabs(phys.jets[ijet].eta())>2.5) continue;
	  prunedJets.push_back(phys.jets[ijet]);
	  float btagDisc=phys.jets[ijet].btag7;
	  nbtags += (btagDisc>0.244);
	  if(isMC)
	    {
	      if(fabs(phys.jets[ijet].flavid)==5) nbtagscor += (btagDisc>0.207);
	      else                                nbtagscor += (btagDisc>0.233);
	    }
	  else nbtagscor += (btagDisc>0.244);
	}
      sort(prunedJets.begin(),prunedJets.end(),sortByBtag);
      sort(phys.leptons.begin(),phys.leptons.end(),sortByPt);

      //compute basic kinematics for the leptons / dileptons / dijet
      LorentzVector tauP4(0,0,0,0);
      if(fabs(phys.leptons[0].id)==15)  tauP4 = phys.leptons[0];
      if(fabs(phys.leptons[1].id)==15)  tauP4 = phys.leptons[1];
      LorentzVector dil = phys.leptons[0]+phys.leptons[1];
      float dilmass = dil.mass();
      LorentzVector dij = prunedJets[0]+prunedJets[1];
      float mjj=dij.M();
      double ptjet1(max(prunedJets[0].pt(),prunedJets[1].pt())), ptjet2(min(prunedJets[0].pt(),prunedJets[1].pt()));
      LorentzVector ptttbar=phys.leptons[0]+phys.leptons[1]+prunedJets[0]+prunedJets[1]+phys.met;
      double ptlep1(max(phys.leptons[0].pt(),phys.leptons[1].pt())), ptlep2(min(phys.leptons[0].pt(),phys.leptons[1].pt()));    

      //ht
      double ht(0);
      for(size_t ijet=0; ijet<prunedJets.size(); ijet++) ht += prunedJets[ijet].pt();
      double sumptlep(phys.leptons[0].pt()+phys.leptons[1].pt());
      double st(sumptlep+phys.met.pt());
      double htlep(st+ht);

      int nRecoBs( (fabs(prunedJets[0].flavid)==5) + (fabs(prunedJets[1].flavid)==5) );
      if(!ev.isSignal) nRecoBs=0;
      int iCorrectComb=0;
      if(nRecoBs>1)
	{
	  //the charge of the generator level matched particles must be opposite
	  int assignCode=(phys.leptons[0].genid*prunedJets[0].flavid);
	  if(assignCode<0) iCorrectComb=1;
	  else             iCorrectComb=2;
	  //  cout << iCorrectComb << " | " 
	  // 	       << phys.leptons[0].genid << " (" << phys.leptons[0].pt() << ") "
	  // 	       << phys.leptons[1].genid << " (" << phys.leptons[1].pt() << ") |"
	  // 	       << prunedJets[0].flavid << " (" << prunedJets[0].btag1 << ") "
	  // 	       << prunedJets[1].flavid << " (" << prunedJets[1].btag1 << ") |"
	  // 	       << endl;
	}
      
      //define the event category
      bool isZcand(fabs(dilmass-91)<15 && (ev.cat==EE || ev.cat==MUMU));
      bool isSS( phys.leptons[0].id*phys.leptons[1].id >0 );
      std::vector<TString> subcats;
      if(!isZcand && !isSS)
	{
	  subcats.push_back("");
	  if(nbtagscor==0) subcats.push_back("eq0btags");
	  else if(nbtagscor==1) subcats.push_back("eq1btags");
	  else if(nbtagscor>=2) subcats.push_back("geq2btags");
	}
      else if(isZcand && !isSS) subcats.push_back("zcands");
      else                      subcats.push_back("ss");
      for(std::vector<TString>::iterator cIt = categs.begin(); cIt != categs.end(); cIt++)
	{
	  for(std::vector<TString>::iterator scIt = subcats.begin(); scIt != subcats.end(); scIt++)
	    {
	      TString ctf=*cIt + *scIt;
	      controlHistos.fillHisto("evtflow",ctf,0,weight);
	    }
	}

      // get the top mass likelihood
      massLikelihood->processEvent(kinHandler);
      // massLikelihood.getEventLikelihood() // eventually store that

      //
      // get preferred combination and the top mass measurement from the MPV fit
      //
      TH1F *h1=kinHandler.getHisto("mt",1), *h2=kinHandler.getHisto("mt",2);
      //h1->Rebin(2); h2->Rebin(2);  //<- don't rebin you'll loose resolution
      Int_t icomb=getPreferredCombination(h1,h2);
      if(icomb<0) continue;
      TH1F *mpref=kinHandler.getHisto("mt",icomb);
      double mtop = kinHandler.getMPVEstimate(mpref) [1];
      if(mtop<=0) continue;
      TH1F *mttbarpref=kinHandler.getHisto("mttbar",icomb);
      double mttbar = kinHandler.getMPVEstimate(mttbarpref)[1];
      //      TH1F *afbpref=kinHandler.getHisto("afb",icomb);
      //      double afb = kinHandler.getMPVEstimate(afbpref)[1];
      LorentzVector lj1=phys.leptons[0]+prunedJets[icomb==1?0:1];
      LorentzVector lj2=phys.leptons[1]+prunedJets[icomb==1?1:0];
      
      //now fill the control plots
      for(std::vector<TString>::iterator cIt = categs.begin(); cIt != categs.end(); cIt++)
	{
	  for(std::vector<TString>::iterator scIt = subcats.begin(); scIt != subcats.end(); scIt++)
	    {
	      TString ctf=*cIt + *scIt;
	      controlHistos.fillHisto("evtflow",ctf,1,weight);
	      controlHistos.fillHisto("evtflow",ctf,2+(nbtagscor>2?2:nbtagscor),weight);
	      
	      controlHistos.fillHisto("taupt",ctf,tauP4.pt(),weight);
	      controlHistos.fillHisto("taueta",ctf,fabs(tauP4.eta()),weight);
	      
	      controlHistos.fillHisto("njets",ctf,prunedJets.size(),weight);
	      controlHistos.fillHisto("btags",ctf,nbtagscor,weight);
	      controlHistos.fillHisto("btagsraw",ctf,nbtags,weight);
	      controlHistos.fillHisto("leadjet",ctf,ptjet1,weight);
	      controlHistos.fillHisto("subleadjet",ctf,ptjet2,weight);
	      controlHistos.fillHisto("leadlepton",ctf,ptlep1,weight);
	      controlHistos.fillHisto("subleadlepton",ctf,ptlep2,weight);
	      controlHistos.fillHisto("met",ctf,phys.met.pt(),weight);
	      controlHistos.fillHisto("ht",ctf,ht,weight);
	      controlHistos.fillHisto("st",ctf,st,weight);
	      controlHistos.fillHisto("sumpt",ctf,sumptlep,weight);
	      controlHistos.fillHisto("htlep",ctf,htlep,weight);
	      controlHistos.fillHisto("ptttbar",ctf,ptttbar.Pt(),weight);
	      
	      controlHistos.fillHisto("mtop",ctf,mtop,weight);
	      controlHistos.fillHisto("mttbar",ctf,mttbar,weight);
	      
	      controlHistos.fillHisto("dilmass",ctf,dilmass,weight);
	      controlHistos.fill2DHisto("mtopvsdilmass",ctf,mtop,dilmass,weight);
	      controlHistos.fill2DHisto("mtopvsmlj",ctf,mtop,lj1.mass(),weight);
	      controlHistos.fill2DHisto("mtopvsmlj",ctf,mtop,lj2.mass(),weight);
	      controlHistos.fill2DHisto("mtopvsmttbar",ctf,mtop,mttbar,weight);
	      
	      //resolution plots
	      controlHistos.fill2DHisto("mtopvsnvtx",ctf,mtop,ev.nvtx,weight);
	      
	      TH2F *h= (TH2F *) controlHistos.getHisto("mtopvsmet",ctf);
	      if(h)
		{
		  for(int ibin=1; ibin<=h->GetYaxis()->GetNbins(); ibin++)
		    {
		      float metmin=h->GetYaxis()->GetBinLowEdge(ibin);
		      if(phys.met.pt()>metmin) h->Fill(mtop,metmin,weight);
		    }
		}
	      
	      h= (TH2F *) controlHistos.getHisto("mtopvsptjet",ctf);
	      if(h)
		{
		  for(int ibin=1; ibin<=h->GetYaxis()->GetNbins(); ibin++)
		    {
		      float ptmin=h->GetYaxis()->GetBinLowEdge(ibin);
		      if(prunedJets[0].pt()>ptmin && prunedJets[1].pt()>ptmin) h->Fill(mtop,ptmin,weight);
		    }
		}
	    }
	}
      neventsused++;
      
      //save for further study if required
      if(!isSS && spyEvents)
	{
	  ev.normWeight=normWeight; ev.xsecWeight=summaryWeight;
	  std::vector<float> measurements;
	  measurements.push_back(mtop);
	  measurements.push_back(mttbar);
	  measurements.push_back(isZcand);
	  measurements.push_back(ptttbar.Pt());
	  measurements.push_back(nbtags);
	  measurements.push_back(nbtagscor);
	  measurements.push_back(prunedJets.size());
	  spyEvents->fillTreeWithEvent( ev, measurements );
	}

      //for data only
      if (!isMC && mttbar>2000 && outf!=0) 
	*outf << "| " << irun << ":" << ilumi << ":" << ievent 
	      << " | " << categs[1] 
	      << " | " << mtop 
	      << " | " << mttbar 
	      << " | " << ptlep1 << ";" << ptlep2  << " | " << dilmass
	      << " | " << ptjet1 << ";" << ptjet2  << " | " << mjj
	      << " | " << phys.met.pt() << " | " << htlep << endl; 
	  
    }


  // store that
  //  controlHistos.addHistogram( &(massLikelihood.getFinalLikelihood()) );
  TH1F blah = massLikelihood->getFinalLikelihood();
  //  TH1F* blah2 = &blah;
  controlHistos.getHisto("totalMassLikelihood","all")->Add( &blah );

  kinHandler.end();
  
  //if MC: rescale to number of selected events and to units of pb
  cout << "From " << selEvents.size() << "original events found " << nresults << " kin results - used " << neventsused << endl; 
  
  if(!isMC && outf!=0) 
    {
      outf->close(); 
      delete outf;
    }

  if(isMC && nresults>0)
    {
      double scaleFactor=double(selEvents.size())/double(nresults);
      if(cnorm>0) 
	{
	  cout << "Will re-scale MC by: " << scaleFactor << "/" << cnorm << endl;
	  cnorm=scaleFactor/cnorm;
	}
    }


  //save to file
  TDirectory *baseOutDir=file->mkdir("localAnalysis");
  SelectionMonitor::StepMonitor_t &mons=controlHistos.getAllMonitors();
  std::map<TString, TDirectory *> outDirs;
  outDirs["all"]=baseOutDir->mkdir("all");
  outDirs["ee"]=baseOutDir->mkdir("ee");
  outDirs["emu"]=baseOutDir->mkdir("emu");
  outDirs["mumu"]=baseOutDir->mkdir("mumu");
  outDirs["etau"]=baseOutDir->mkdir("etau");
  outDirs["mutau"]=baseOutDir->mkdir("mutau");
  for(SelectionMonitor::StepMonitor_t::iterator it =mons.begin(); it!= mons.end(); it++)
    {
      TString icat("all");
      if(it->first.BeginsWith("ee")) icat="ee";
      if(it->first.BeginsWith("emu")) icat="emu";
      if(it->first.BeginsWith("mumu")) icat="mumu";
      if(it->first.BeginsWith("etau")) icat="etau";
      if(it->first.BeginsWith("mutau")) icat="mutau";
      outDirs[icat]->cd();
      for(SelectionMonitor::Monitor_t::iterator hit=it->second.begin(); hit!= it->second.end(); hit++)
	{
	  if(isMC && cnorm>0) hit->second->Scale(cnorm);
	  if( !((TClass*)hit->second->IsA())->InheritsFrom("TH2")
	      && !((TClass*)hit->second->IsA())->InheritsFrom("TGraph") )
	    fixExtremities(hit->second,true,true);
	  hit->second->Write();

	}
    }
  file->Close(); 

  //close the spy events
  if(spyEvents)
    {
      spyDir->cd();
      spyEvents->getTree()->Write();
      spyFile->Write();
      spyFile->Close();
    }
}  
