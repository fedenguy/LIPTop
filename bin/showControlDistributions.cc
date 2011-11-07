
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "LIP/Top/interface/EventSummaryHandler.h"

#include "CMGTools/HtoZZ2l2nu/interface/TMVAUtils.h"
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

using namespace std;
using namespace top;


//
bool isQQlike(LorentzVectorCollection &jets, LorentzVector &l1, LorentzVector &l2)
{
  int NJetsQQ(0);
  LorentzVector QQj1(0,0,0,0),QQj2(0,0,0,0);
  LorentzVectorCollection centralJets;
  for(size_t ijet=0; ijet<jets.size(); ijet++) 
    {
      float jpt  = jets[ijet].pt();
      float jeta = jets[ijet].eta();
      if(jpt<30)continue;
      if(fabs(jeta)<2.5)
	{
	  centralJets.push_back(jets[ijet]); 
	  continue;
	}
      
      NJetsQQ++;
      if(jets[ijet].pt()>QQj1.pt()) 
	{
	  QQj2=QQj1;    
	  QQj1=jets[ijet];
	}
      else if (jets[ijet].pt()>QQj2.pt())
	{
	  QQj2=jets[ijet];
	}
      
    }
  
  bool isQQ = false;
  if(NJetsQQ>=2)
    {
      LorentzVector QQSyst = QQj1+QQj2;
      double j1eta=QQj1.eta();
      double j2eta=QQj2.eta();
      double dEta = fabs(j1eta-j2eta);
	
      double MaxEta, MinEta;
      if(j1eta<j2eta) { MinEta=j1eta; MaxEta=j2eta;}
      else            { MinEta=j2eta; MaxEta=j1eta;}

      int NCentralLepton(0);
      if(l1.eta()>MinEta && l1.eta()<MaxEta) NCentralLepton++;
      if(l2.eta()>MinEta && l2.eta()<MaxEta) NCentralLepton++;

      int NCentralJet(0);
      for(size_t icentralJets=0; icentralJets<centralJets.size(); icentralJets++)
	if(centralJets[icentralJets].eta()>MinEta && centralJets[icentralJets].eta()<MaxEta)
	  NCentralJet++;
      
      
      isQQ = (dEta>3.5) && (QQSyst.M()>450) && (NCentralLepton==2) && (NCentralJet==2);
    }
  
  return isQQ;
}


float getArcCos(LorentzVector &a, LorentzVector &b)
{
  TVector3 mom1(a.px(),a.py(),a.pz());
  TVector3 mom2(b.px(),b.py(),b.pz());
  double cosine = mom1.Dot(mom2)/(mom1.Mag()*mom2.Mag());
  double arcCosine = acos(cosine)-TMath::Pi();
  return arcCosine;
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
  TString uncFile =  runProcess.getParameter<std::string>("jesUncFileName"); gSystem->ExpandPathName(uncFile);

  //pileup reweighter
  reweight::PoissonMeanShifter PShiftUp(+0.6);
  reweight::PoissonMeanShifter PShiftDown(-0.6);

  // b-tag working points: https://twiki.cern.ch/twiki/bin/view/CMS/BTagPerformanceOP
  std::map<TString,float> btagCuts;
  btagCuts["TCHEL"]=1.7;
  btagCuts["TCHEM"]=3.3;
  btagCuts["TCHET"]=10.2;
  btagCuts["TCHPM"]=1.93;
  btagCuts["TCHPT"]=3.41;
  btagCuts["JBPL"]=1.33;
  btagCuts["JBPM"]=2.55;
  btagCuts["JBPT"]=3.74;
  btagCuts["JPL"]=0.275;
  btagCuts["JPM"]=0.545;
  btagCuts["JPT"]=0.790;
  btagCuts["SSVHEM"]=1.74;
  btagCuts["SSVHET"]=3.05;
  btagCuts["SSVHPT"]=2.00;
  btagCuts["CSVL"]=0.244;
  btagCuts["CSVM"]=0.679;
  btagCuts["CSVT"]=0.898;
   
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
  double massAxis[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,400,500,1000,2000};
  int nMassBins=sizeof(massAxis)/sizeof(double)-1;
  controlHistos.addHistogram( new TH1D("mlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV/c^{2}];Lepton-jet pairs",nMassBins,massAxis) );
  controlHistos.addHistogram( new TH1D("drlj","Lepton-jet spectrum;#Delta R(l,j);Lepton-jet pairs",50,0,6) );
  controlHistos.addHistogram( new TH1D("mindrlj","Lepton-jet spectrum;min #Delta R(l,j);Events",50,0,6) );
  
  controlHistos.addHistogram( new TH1D("dilmass",";M(l,l') [GeV/c^{2}];Events",100,0,250) );
  TH1D *lepMult=new TH1D("nleptons",";Leptons;Events",3,0,3);
  lepMult->GetXaxis()->SetBinLabel(1,"=2 leptons");
  lepMult->GetXaxis()->SetBinLabel(2,"=3 leptons");
  lepMult->GetXaxis()->SetBinLabel(3,"#geq 4 leptons");
  controlHistos.addHistogram( lepMult );
  cout << lepMult << " " << lepMult->IsA()->InheritsFrom("TH1") << endl;
  TH1 *sslepMult = (TH1 *) lepMult->Clone("ssnleptons");

  controlHistos.addHistogram( sslepMult );
  controlHistos.addHistogram( new TH1D("dilarccosine",";arcCos(l,l');Events",100,-3.2,3.2) );
  controlHistos.addHistogram( new TH1D("dilcharge",";Charge;Events",3,-1.5,1.5) );
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
  controlHistos.addHistogram( new TH1F ("mt", "; M_{T}(l^{1},E_{T}^{miss}) [GeV/c^{2}]; Events / (25 GeV/c^{2})", 10, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 25, 0.,25.) );
  TH1F * h = new TH1F ("njets", "; Jet multiplicity; Events", 4, 0.,4.);
  h->GetXaxis()->SetBinLabel(1,"=2 jets");
  h->GetXaxis()->SetBinLabel(2,"=3 jets");
  h->GetXaxis()->SetBinLabel(3,"=4 jets");
  h->GetXaxis()->SetBinLabel(4,"#geq 5 jets");
  controlHistos.addHistogram( h );

  controlHistos.addHistogram( new TH1F ("mljmostiso",";M(l,j) [GeV/c^{2}]; Events",100,0,400) );
  controlHistos.addHistogram( new TH1F ("mljleastiso",";M(l,j) [GeV/c^{2}]; Events",100,0,400) );
  controlHistos.addHistogram( new TH1F ("ptmostisollep",";p_{T} [GeV/c]; Events ",20,0,200) );
  controlHistos.addHistogram( new TH1F ("ptleastisollep",";p_{T} [GeV/c]; Events ",20,0,200) );
  controlHistos.addHistogram( new TH1F ("mtmostisolep",";M_{T}(l,E_{T}^{miss}) [GeV/c^{2}]; Events",20,0,400) );
  controlHistos.addHistogram( new TH1F ("mtleastisolep",";M_{T}(l,E_{T}^{miss}) [GeV/c^{2}]; Events",20,0,400) );
  controlHistos.addHistogram( new TH1F ("dphimetmostisollep",";#Delta#phi(l,E_{T}^{miss}) [rad]; Events",100,0,6) );
  controlHistos.addHistogram( new TH1F ("dphimetleastisolep",";#Delta#phi(l,E_{T}^{miss}) [rad]; Events",100,0,6) );

  TString btagalgos[]={"TCHE","TCHP","SSVHE","SSVHP","JBP","JP","CSV"};
  float minba[]      ={-5    ,-5    ,0      ,0      ,0    ,0   ,0    };
  float maxba[]      ={25    ,25    ,6      ,6      ,8    ,3   ,1.2  };
  int nbtagalgos=sizeof(btagalgos)/sizeof(TString);
  TString btagwps[]={"L","M","T"};
  int nbtagwps=sizeof(btagwps)/sizeof(TString);
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
  controlHistos.addHistogram( new TH1F("btagdownvar",";b-tag multiplicity;Events", 3, 0.,3.) );
  controlHistos.addHistogram( new TH1F("btagupvar",";b-tag multiplicity;Events", 3, 0.,3.) );

  //event selection histogram
  enum SelectionSteps { SEL2LEPTONS, SELDILEPTON, SELJETS, SELMET, SELOS, SEL0BTAGS, SEL1BTAGS, SEL2BTAGS, SELMT, SELTIGHTLEPTONS, SELTIGHTJETS };
  TString labels[]={"2 leptons",
		    "M>12 #wedge |M-M_{Z}|>15",
		    "#geq 2 jets",
		    "E_{T}^{miss}>30,0",
		    "op. sign",
		    "=0 b-tags",
		    "=1 b-tag",
		    "#geq 2 b-tags",
		    "#Sigma M_{T}(l,E_{T}^{miss})>75 GeV/c^{2}",
		    "p_{T}^{l}>40",
		    "p_{T}^{jets}>40"
  };
  int nsteps=sizeof(labels)/sizeof(TString);
  TString cats[]={"","jer","jesdown","jesup","pudown","puup"};
  int nvarcats=sizeof(cats)/sizeof(TString); 
  if(!runSystematics) nvarcats=1; 
  else { cout << "Running sytematics: this will take a while..." << endl; }
  for(int ivar=0;ivar<nvarcats; ivar++) 
    {
      TH1D *cutflowH=new TH1D("evtflow"+cats[ivar],";Cutflow;Events",nsteps,0,nsteps);
      for(int ibin=0; ibin<nsteps; ibin++) cutflowH->GetXaxis()->SetBinLabel(ibin+1,labels[ibin]);
      controlHistos.addHistogram( cutflowH );
      controlHistos.addHistogram( new TH1D("dilmassctr"+cats[ivar],";Region;Events",2,0,2) );
      controlHistos.addHistogram( new TH1D("mtsum"+cats[ivar],";M_{T}(l^{(1)},E_{T}^{miss})+M_{T}(l^{(2)},E_{T}^{miss});Events",100,0,1000) );
      controlHistos.addHistogram( new TH1D("ptsum"+cats[ivar],";p_{T}(l^{(1)})+p_{T}(l^{(2)});Events",100,0,500) );


      for(int iba=0; iba<nbtagalgos; iba++)
	{
	  for(int ibwp=0; ibwp<nbtagwps; ibwp++)
	    {
	      TString key=btagalgos[iba]+btagwps[ibwp];
	      if(btagCuts.find(key)==btagCuts.end()) continue;
	      TH1F *bmultH=new TH1F (key+"full"+cats[ivar], ";b-tag multiplicity;Events", 15, 0.,15.) ;
	      for(int ibin=1; ibin<=5; ibin++)
		{
		  TString label(""); label += ibin-1;
		  if(ibin==bmultH->GetXaxis()->GetNbins()) label ="#geq" + label;
		  bmultH->GetXaxis()->SetBinLabel(ibin,label + " btags");
		  bmultH->GetXaxis()->SetBinLabel(ibin+5,label + " btags");
		  bmultH->GetXaxis()->SetBinLabel(ibin+10,label + " btags");
		}
	      controlHistos.addHistogram( bmultH );
	    }
	}
    }

  //TMVA configuration
  TMVA::Reader *tmvaReader = 0;
  bool useMVA                             = runProcess.getParameter<bool>("useMVA");
  std::vector<std::string> methodList,varsList;
  std::vector<int> evCategories;
  std::vector<Float_t> tmvaVars,discriResults;
  if(useMVA)
    {
      try{
	edm::ParameterSet tmvaInput             = runProcess.getParameter<edm::ParameterSet>("tmvaInput");
	methodList     = tmvaInput.getParameter<std::vector<std::string> >("methodList");
	varsList       = tmvaInput.getParameter<std::vector<std::string> >("varsList");
	evCategories           = tmvaInput.getParameter<std::vector<int> >("evCategories");
	discriResults.resize(methodList.size(),0);
	tmvaVars.resize(varsList.size()+1,0);
	std::string weightsDir                  = tmvaInput.getParameter<std::string>("weightsDir");
	std::string studyTag                    = tmvaInput.getParameter<std::string>("studyTag");
		
	std::cout << "==> Start TMVA Classification with " << methodList.size() << " methods and " << varsList.size() << " variables" << std::endl;
	
	//start the reader for the variables and methods
	tmvaReader = new TMVA::Reader( "!Color:!Silent" );
	for(size_t ivar=0; ivar<varsList.size(); ivar++)   tmvaReader->AddVariable( varsList[ivar], &tmvaVars[ivar] );
	tmvaReader->AddSpectator("eventCategory", &tmvaVars[varsList.size()]);
	
	//book the methods
	for(size_t imet=0; imet<methodList.size(); imet++)
	  {
	    //open the file with the method description
	    TString weightFile = weightsDir + "/" + studyTag + ( evCategories.size()>1 ? "_Category_" : "_" ) + methodList[imet] + TString(".weights.xml");
	    gSystem->ExpandPathName(weightFile);
	    
	    tmvaReader->BookMVA(methodList[imet], weightFile);
	    TH1 *h=tmva::getHistogramForDiscriminator( methodList[imet] );
	    controlHistos.addHistogram( h );
	    controlHistos.addHistogram( (TH1 *) h->Clone(methodList[imet]+TString("Tight")) );
	    controlHistos.addHistogram( (TH1 *) h->Clone(methodList[imet]+TString("0btagsTight")) );
	    controlHistos.addHistogram( (TH1 *) h->Clone(methodList[imet]+TString("1btagsTight")) );
	    controlHistos.addHistogram( (TH1 *) h->Clone(methodList[imet]+TString("2btagsTight")) );
	  }
      }catch(std::exception &e){
	useMVA=false;
	cout << "[Warning] disabling MVA - caught exception : " << e.what() << endl;
      }
    }
  
  controlHistos.initMonitorForStep("ee");
  controlHistos.initMonitorForStep("emu");
  controlHistos.initMonitorForStep("mumu");
  controlHistos.initMonitorForStep("qqlike");
  
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
  //efficiency (trigger,id+isolation) corrections for the different channels
  std::map<int,float> effCorr;
  effCorr[MUMU] = 0.92*pow(1.0,2);
  effCorr[EMU]  = 1.0*1.0*0.96;
  effCorr[EE]   = 1.0*pow(0.96,2);
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
      top::EventSummary_t &ev = evSummaryHandler.getEvent();
      if(isMC && mcTruthMode>0)
	{
	  if(mcTruthMode==1 && !ev.isSignal) continue;
	  if(mcTruthMode==2 && ev.isSignal)  continue;
	}
      
      if(duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event)){ 
	NumberOfDuplicated++; continue; 
      }
      float puweight = ev.weight/cnorm;    
      float normWeight = ev.normWeight;

      //get particles from event
      top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      
      bool isSameFlavor(false);
      TString ch("");
      if(ev.cat==MUMU)  { isSameFlavor=true; ch="mumu"; }
      if(ev.cat==EE)    { isSameFlavor=true; ch="ee"; }
      if(ev.cat==EMU)   ch="emu";
      if(ev.cat==ETAU)  ch="etau";
      if(ev.cat==MUTAU) ch="mutau";

      //efficiency correction for MC
      if(isMC) puweight *= effCorr[ev.cat];

      //b-tag analysis
      bcomp.compute(phys.nbjets,phys.nljets);
      std::vector<btag::Weight_t> wgt = bcomp.getWeights();
      double p0btags = wgt[0].first;  double p0btags_err=wgt[0].second;
      double p1btags = wgt[1].first;  double p1btags_err=wgt[1].second;
      double p2btags = wgt[2].first;  double p2btags_err=wgt[2].second;

      //order jets to apply variations
      top::PhysicsObjectJetCollection orderedJetColl = phys.jets;
      sort(orderedJetColl.begin(),orderedJetColl.end(),top::PhysicsEvent_t::sortJetsByBtag);
      LorentzVectorCollection jets;
      std::vector<int> baseJetsIdx;
      for(size_t ijet=0; ijet<orderedJetColl.size(); ijet++)
	{
	  float pt=orderedJetColl[ijet].pt();
	  float eta=orderedJetColl[ijet].eta();
	  if(pt<20 || fabs(eta)>2.5) continue;
	  jets.push_back(orderedJetColl[ijet]);
	  baseJetsIdx.push_back(ijet);
	}

      
      //jet/met variations
      std::vector<LorentzVectorCollection> jetsVar;
      LorentzVectorCollection metsVar;
      jet::computeVariation(jets,phys.met,jetsVar,metsVar,&stdPtResol,&stdEtaResol,&stdPhiResol,&jecUnc);
      LorentzVectorCollection jerVariedJets     = jetsVar[jet::JER];
      LorentzVector jerVariedMet                = metsVar[jet::JER];
      LorentzVectorCollection jesupVariedJets   = jetsVar[jet::JES_UP];
      LorentzVector jesUpVariedMet              = metsVar[jet::JES_UP];
      LorentzVectorCollection jesdownVariedJets = jetsVar[jet::JES_DOWN];
      LorentzVector jesDownVariedMet            = metsVar[jet::JES_DOWN];    

      //the cutflow with variations
      for(int ivar=0;ivar<(isMC ? nvarcats : 1); ivar++) 
	{
	  float weight=puweight;
	  //objects to use
	  LorentzVectorCollection jetColl = jets;
	  LorentzVector theMET           = phys.met;     
	  if(cats[ivar]=="jer")      { jetColl=jerVariedJets;      theMET=jerVariedMet;     }
	  if(cats[ivar]=="jesdown")  { jetColl=jesdownVariedJets;  theMET=jesUpVariedMet;   }
	  if(cats[ivar]=="jesup" )   { jetColl=jesupVariedJets;    theMET=jesDownVariedMet; }
	  if(cats[ivar]=="puup" )    { double TotalWeight_plus  = PShiftUp.ShiftWeight( ev.ngenpu );   weight *= TotalWeight_plus;  }
	  if(cats[ivar]=="pudown" )  { double TotalWeight_minus = PShiftDown.ShiftWeight( ev.ngenpu ); weight *= TotalWeight_minus; }

	  LorentzVector l1             = (phys.leptons[0].pt() > phys.leptons[1].pt() ? phys.leptons[0] : phys.leptons[1]);
	  LorentzVector l2             = (phys.leptons[0].pt() < phys.leptons[1].pt() ? phys.leptons[0] : phys.leptons[1]);	  
	  float dilcharge=(phys.leptons[0].id/fabs(phys.leptons[0].id)) *(phys.leptons[1].id/fabs(phys.leptons[1].id));
	  LorentzVector dileptonSystem = l1+l2;
	  LorentzVector jetSystem      = (jetColl.size() > 1 ? jetColl[0]+jetColl[1] : LorentzVector(0,0,0,0));
	  LorentzVector visible_t      = jetSystem+dileptonSystem;
	  LorentzVector ttbar_t        = jetSystem+dileptonSystem+theMET;

	  std::vector<TString> catsToFill;
	  catsToFill.push_back("all");
	  catsToFill.push_back(ch);
	  if(isQQlike(jetColl,l1,l2)) catsToFill.push_back("qqlike");
	  const size_t nCatsToFill=catsToFill.size();
	  	  
	  //kinematics with propagated variations
	  int nseljetsLoose(0), nseljetsTight(0);
	  LorentzVectorCollection prunedJetColl;
	  int nGoodassigments(0);
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
	  std::vector< std::pair<LorentzVector,double> > correctPairKin, wrongPairKin; 
	  for(size_t ijet=0; ijet<jetColl.size(); ijet++) 
	    {

	      float pt=jetColl[ijet].pt();
	      if(pt<30 || fabs(jetColl[ijet].eta())>2.5) continue;
	      nseljetsLoose ++;
	      prunedJetColl.push_back(jetColl[ijet]);
	      nseljetsTight += (pt>40 && fabs(jetColl[ijet].eta())<2.5);

	      //idx for the original jet
	      int baseIdx = baseJetsIdx[ijet];

	      //check MC truth flavor
	      bool hasBflavor(fabs(orderedJetColl[baseIdx].flavid)==5);
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

	      //check correct assignment
	      int assignCode1=(phys.leptons[0].genid*orderedJetColl[baseIdx].genid);
	      bool isCorrect1( (assignCode1<0) && hasBflavor);
	      LorentzVector mlj=phys.leptons[0]+orderedJetColl[baseIdx]; 
	      double dthetalj=getArcCos(phys.leptons[0],orderedJetColl[baseIdx]); 
	      if(isCorrect1) correctPairKin.push_back(std::pair<LorentzVector,double> (mlj,dthetalj) );
	      else           wrongPairKin.push_back(std::pair<LorentzVector,double> (mlj,dthetalj) );

	      int assignCode2=(phys.leptons[1].genid*orderedJetColl[baseIdx].genid);
	      bool isCorrect2 = ( (assignCode2<0) && hasBflavor);
	      mlj=phys.leptons[1]+orderedJetColl[baseIdx]; 
	      dthetalj=getArcCos(phys.leptons[1],orderedJetColl[baseIdx]); 
	      if(isCorrect2) correctPairKin.push_back(std::pair<LorentzVector,double> (mlj,dthetalj) );
	      else           wrongPairKin.push_back(std::pair<LorentzVector,double> (mlj,dthetalj) );
	      
	      nGoodassigments+=(isCorrect1||isCorrect2);
	      
	      //count b-tags
	      nbtags["TCHEL"]  += (orderedJetColl[baseIdx].btag1>btagCuts["TCHEL"]);
	      nbtags["TCHEM"]  += (orderedJetColl[baseIdx].btag1>btagCuts["TCHEM"]);
	      nbtags["TCHET"]  += (orderedJetColl[baseIdx].btag1>btagCuts["TCHET"]);
	      nbtags["TCHPM"]  += (orderedJetColl[baseIdx].btag2>btagCuts["TCHPM"]);
	      nbtags["TCHPT"]  += (orderedJetColl[baseIdx].btag2>btagCuts["TCHPT"]);
	      nbtags["SSVHEM"] += (orderedJetColl[baseIdx].btag3>btagCuts["SSVHEM"]);		
	      nbtags["SSVHET"] += (orderedJetColl[baseIdx].btag3>btagCuts["SSVHET"]);		
	      nbtags["SSVHPT"] += (orderedJetColl[baseIdx].btag6>btagCuts["SSVHPT"]);		
	      nbtags["JBPL"]   += (orderedJetColl[baseIdx].btag4>btagCuts["JBPL"]);
	      nbtags["JBPM"]   += (orderedJetColl[baseIdx].btag4>btagCuts["JBPM"]);
	      nbtags["JBPT"]   += (orderedJetColl[baseIdx].btag4>btagCuts["JBPT"]);
	      nbtags["JPL"]   += (orderedJetColl[baseIdx].btag5>btagCuts["JPL"]);
	      nbtags["JPM"]   += (orderedJetColl[baseIdx].btag5>btagCuts["JPM"]);
	      nbtags["JPT"]   += (orderedJetColl[baseIdx].btag5>btagCuts["JPT"]);
	      nbtags["CSVL"]   += (orderedJetColl[baseIdx].btag7>btagCuts["CSVL"]);
	      nbtags["CSVM"]   += (orderedJetColl[baseIdx].btag7>btagCuts["CSVM"]);
	      nbtags["CSVT"]   += (orderedJetColl[baseIdx].btag7>btagCuts["CSVT"]);
	  
	      if(isMC)
		{
		  //b-tag counting per jet flavor and pt-bin
		  for(size_t iptbin=0; iptbin<2; iptbin++)
		    {
		      nFlavBtags[ptbins[iptbin]]["TCHEL"+flavCat]  += (orderedJetColl[baseIdx].btag1>btagCuts["TCHEL"]);
		      nFlavBtags[ptbins[iptbin]]["TCHEM"+flavCat]  += (orderedJetColl[baseIdx].btag1>btagCuts["TCHEM"]);
		      nFlavBtags[ptbins[iptbin]]["TCHET"+flavCat]  += (orderedJetColl[baseIdx].btag1>btagCuts["TCHET"]);
		      nFlavBtags[ptbins[iptbin]]["TCHPM"+flavCat]  += (orderedJetColl[baseIdx].btag2>btagCuts["TCHPM"]);
		      nFlavBtags[ptbins[iptbin]]["TCHPT"+flavCat]  += (orderedJetColl[baseIdx].btag2>btagCuts["TCHPT"]);
		      nFlavBtags[ptbins[iptbin]]["SSVHEM"+flavCat] += (orderedJetColl[baseIdx].btag3>btagCuts["SSVHEM"]);		
		      nFlavBtags[ptbins[iptbin]]["SSVHET"+flavCat] += (orderedJetColl[baseIdx].btag3>btagCuts["SSVHET"]);		
		      nFlavBtags[ptbins[iptbin]]["SSVHPT"+flavCat] += (orderedJetColl[baseIdx].btag6>btagCuts["SSVHPT"]);		
		      nFlavBtags[ptbins[iptbin]]["JBPL"+flavCat]   += (orderedJetColl[baseIdx].btag4>btagCuts["JBPL"]);
		      nFlavBtags[ptbins[iptbin]]["JBPM"+flavCat]   += (orderedJetColl[baseIdx].btag4>btagCuts["JBPM"]);
		      nFlavBtags[ptbins[iptbin]]["JBPT"+flavCat]   += (orderedJetColl[baseIdx].btag4>btagCuts["JBPT"]);
		      nFlavBtags[ptbins[iptbin]]["JPL"+flavCat]   += (orderedJetColl[baseIdx].btag5>btagCuts["JPL"]);
		      nFlavBtags[ptbins[iptbin]]["JPM"+flavCat]   += (orderedJetColl[baseIdx].btag5>btagCuts["JPM"]);
		      nFlavBtags[ptbins[iptbin]]["JPT"+flavCat]   += (orderedJetColl[baseIdx].btag5>btagCuts["JPT"]);
		      nFlavBtags[ptbins[iptbin]]["CSVL"+flavCat]   += (orderedJetColl[baseIdx].btag7>btagCuts["CSVL"]);
		      nFlavBtags[ptbins[iptbin]]["CSVM"+flavCat]   += (orderedJetColl[baseIdx].btag7>btagCuts["CSVM"]);
		      nFlavBtags[ptbins[iptbin]]["CSVT"+flavCat]   += (orderedJetColl[baseIdx].btag7>btagCuts["CSVT"]);
		    }
		}
	    }
	  jetColl=prunedJetColl;

	  double acosine      = getArcCos(l1,l2);
	  double ptsum        = l1.pt()+l2.pt();
	  double drll         = deltaR(l1,l2);
	  double dphill       = deltaPhi(l1.phi(),l2.phi());
	  double st           = ptsum+theMET.pt();
	  double leadlepmt    = METUtils::transverseMass(l1,theMET,false);
	  double subleadlepmt = METUtils::transverseMass(l2,theMET,false);
	  double mtsum        = leadlepmt+subleadlepmt;
	  double mt           = METUtils::transverseMass(visible_t,theMET,false);
	  int btagbin(nbtags["JBPL"]);
	  if(btagbin>2) btagbin=2;
	  TString btagbinLabel(""); btagbinLabel += btagbin; btagbinLabel+="btags";
	  int njetbin(jetColl.size()-2);
	  if(njetbin>=2) njetbin=2;


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
	  
	  //lepton-jet pairs with leading b-ranked jets
	  LorentzVector l1ji[2]={LorentzVector(0,0,0,0),LorentzVector(0,0,0,0)};
	  LorentzVector l2ji[2]={LorentzVector(0,0,0,0),LorentzVector(0,0,0,0)};
	  if(nseljetsLoose>1) {
	    l1ji[0]=l1+jetColl[0]; l1ji[1]=l1+jetColl[1];
	    l2ji[0]=l2+jetColl[0]; l2ji[1]=l2+jetColl[1];
	  }
	  float minL1ji = min( l1ji[0].mass(), l1ji[1].mass() );
	  float minL2ji = min( l2ji[0].mass(), l2ji[1].mass() );
	  LorentzVector mostIsolPair,   leastIsolPair;
	  LorentzVector mostIsolLepton, leastIsolLepton;
	  float mostIsolDr(0), leastIsolDr(0);
	  if(nseljetsLoose)
	    {
	      if(minL1ji < minL2ji)
		{
		  mostIsolLepton=l1;
		  leastIsolLepton=l2;
		  int iComb1=(l1ji[0].mass()<l1ji[1].mass() ? 0 : 1);
		  int iComb2=(iComb1==0 ? 1 : 0); 
		  mostIsolPair=l1ji[iComb1]; 
		  leastIsolPair=l2ji[iComb2]; 
		  mostIsolDr = deltaR(l1,jetColl[iComb1]);
		  leastIsolDr = deltaR(l2,jetColl[iComb2]);
		}
	      else
		{
		  leastIsolLepton=l1;
		  mostIsolLepton=l2;
		  int iComb2=(l2ji[0].mass()<l2ji[1].mass() ? 0 : 1);
		  int iComb1=(iComb2==0 ? 1 : 0);
		  mostIsolPair=l2ji[iComb2];
		  leastIsolPair=l1ji[iComb1];
		  leastIsolDr = deltaR(l1,jetColl[iComb1]);
		  mostIsolDr = deltaR(l2,jetColl[iComb2]);
		}
	    }
	  float mostIsolMt = METUtils::transverseMass(mostIsolLepton,theMET,false);
	  float leastIsolMt = METUtils::transverseMass(leastIsolLepton,theMET,false);

	  //MVA evaluation                                                                                                                                       
	  if(useMVA)
	    {
	      for(size_t ivar=0; ivar<varsList.size(); ivar++)
		{
		  std::string variable=varsList[ivar];
		  if(variable=="mindrlj")              tmvaVars[ivar] = mindrlj;
		  if(variable=="mljmostiso")           tmvaVars[ivar] = mostIsolPair.mass();
		  if(variable=="mljleastiso")          tmvaVars[ivar] = leastIsolPair.mass();
		  if(variable=="ptmostisollep")        tmvaVars[ivar] = mostIsolLepton.pt();
		  if(variable=="ptleastisollep")       tmvaVars[ivar] = leastIsolLepton.pt();
		  if(variable=="mtmostisolep")         tmvaVars[ivar] = mostIsolMt;
		  if(variable=="mtleastisolep")        tmvaVars[ivar] = leastIsolMt;
		  if(variable=="drmostisol")           tmvaVars[ivar] = mostIsolDr;
		  if(variable=="drleastisol")          tmvaVars[ivar] = leastIsolDr;
		}
	      tmvaVars[varsList.size()] = btagbin;
	      for(size_t imet=0; imet<methodList.size(); imet++)
		discriResults[imet]=tmvaReader->EvaluateMVA( methodList[imet] );
	    }

	  bool isInQuarkoniaRegion( dileptonSystem.mass()<12 );
	  bool isZcand(isSameFlavor && fabs(dileptonSystem.mass()-91)<15);
	  bool passLooseJets(nseljetsLoose>1);
	  bool passMet( !isSameFlavor || (theMET.pt()>30) );
	  bool isOS(dilcharge<0);
	  bool passMtCut(mtsum>75);
	  bool passTightLepton(min(l1.pt(),l2.pt())>40);
	  bool passTightJets(nseljetsTight>1);
	  
	  //fill selection histograms
	  for(size_t ictf=0; ictf<nCatsToFill; ictf++)
	    {
	      TString ctf=catsToFill[ictf];
	      if(isInQuarkoniaRegion) continue;
	      if(!passLooseJets)      continue;
	      
	      // >= 2 jets
	      if(!isZcand)
		{
		  controlHistos.fillHisto("evtflow"+cats[ivar],ctf,SELJETS,weight);
		  if(ivar==0) 
		    {
		      controlHistos.fillHisto("met",ctf,theMET.pt(),weight);
		      controlHistos.fillHisto("st", ctf,st,weight);
		    }
		}

	      //MET
	      if(!passMet) continue;
	      if(!isZcand)
		{
		  controlHistos.fillHisto("evtflow"+cats[ivar],ctf,SELMET,weight);
		  if(ivar==0) 
		    {
		      controlHistos.fillHisto("dilcharge",ctf,dilcharge,weight);
		      controlHistos.fillHisto("dilarccosine",ctf,acosine,weight);
		    }
		}
 
	      // OS dilepton
	      if(!isOS) continue;
	      if(ivar==0)  controlHistos.fillHisto("dilmass",ctf,dileptonSystem.mass(),weight);
	      controlHistos.fillHisto("dilmassctr"+cats[ivar],ctf,isZcand,weight);
	      if(!isZcand)
		{
		  selEvents+=weight;
		  controlHistos.fillHisto("evtflow"+cats[ivar],ctf,SELOS,weight);
		  controlHistos.fillHisto("evtflow"+cats[ivar],ctf,SEL0BTAGS+btagbin,weight);
		  controlHistos.fillHisto("mtsum"+cats[ivar],ctf,mtsum,weight);
		  controlHistos.fillHisto("ptsum"+cats[ivar],ctf,ptsum,weight);
		  
		  //b-tag multiplicity 
		  for(std::map<TString,int>::iterator it = nbtags.begin(); it!= nbtags.end(); it++)
		    {
		      int btagCtr=it->second;
		      controlHistos.fillHisto(it->first+"full"+cats[ivar],ctf,njetbin*5+btagCtr,weight);
		      if(ivar==0)
			{
			  controlHistos.fillHisto(it->first,ctf,btagCtr,weight);
			  controlHistos.fillHisto(it->first+(ev.nvtx<6 ? "lowpu" :"highpu"),ctf,btagCtr,weight);
			}
		      
		    }
		  
		  
		  if(ivar==0)
		    {
		      if(isMC)
			{
			  //uncertainty on b-tag selection
			  double p0weight = p0btags+p0btags_err; 
			  double p1weight = p1btags+p1btags_err; 
			  double p2weight = p2btags+p2btags_err;   
			  controlHistos.fillHisto("btagupvar",ctf,0,weight*p0weight);
			  controlHistos.fillHisto("btagupvar",ctf,1,weight*p1weight);
			  controlHistos.fillHisto("btagupvar",ctf,2,weight*p2weight);
			  
			  p0weight = p0btags-p0btags_err; 
			  p1weight = p1btags-p1btags_err; 
			  p2weight = p2btags-p2btags_err;
			  controlHistos.fillHisto("btagdownvar",ctf,0,weight*p0weight);
			  controlHistos.fillHisto("btagdownvar",ctf,1,weight*p1weight);
			  controlHistos.fillHisto("btagdownvar",ctf,2,weight*p2weight);
			  
			  //b/mis-tag efficiency mc truth
			  for(size_t iptbin=0; iptbin<nptbins; iptbin++)	
			    for(std::map<TString,int>::iterator it = nFlavBtags[iptbin].begin(); it!= nFlavBtags[iptbin].end(); it++)
			      {
				int nmatchedflav=nbs[iptbin];
				if(it->first.EndsWith("udscg") )  nmatchedflav=nudscg[iptbin];
				int ntagged=it->second;
				int nnottagged=nmatchedflav-ntagged;
				controlHistos.fillHisto(it->first,ctf,iptbin*2,nnottagged*weight);
				controlHistos.fillHisto(it->first,ctf,iptbin*2+1,ntagged*weight);
			      }

			  //disriminators per jet flavor
			  for(size_t ijet=0; ijet<orderedJetColl.size(); ijet++)
			    {
			      bool hasBflavor(fabs(orderedJetColl[ijet].flavid)==5);
			      int jetFlavBin(2);
			      if(hasBflavor) jetFlavBin=0;
			      else if(fabs(orderedJetColl[ijet].flavid)==4) jetFlavBin=1;
			      TString flavCat(( hasBflavor ? "b" : "udscg" ) );
			      controlHistos.fillHisto("jetflav",ctf,jetFlavBin,puweight);
			      controlHistos.fillHisto("TCHE"+flavCat, ctf,orderedJetColl[ijet].btag1,puweight);
			      controlHistos.fillHisto("TCHP"+flavCat,ctf,orderedJetColl[ijet].btag2,puweight);
			      controlHistos.fillHisto("SSVHE"+flavCat,ctf,orderedJetColl[ijet].btag3,puweight);
			      controlHistos.fillHisto("JBP"+flavCat,ctf,orderedJetColl[ijet].btag4,puweight);
			      controlHistos.fillHisto("JP"+flavCat,ctf,orderedJetColl[ijet].btag5,puweight);
			      controlHistos.fillHisto("SSVHP"+flavCat,ctf,orderedJetColl[ijet].btag6,puweight);
			      controlHistos.fillHisto("CSV"+flavCat,ctf,orderedJetColl[ijet].btag7,puweight);
			    }
			}


		      //kinematics
		      controlHistos.fillHisto("nvertices",ctf,ev.nvtx,weight,true);
		      controlHistos.fillHisto("pttbar",ctf,ttbar_t.pt(),weight);
		      controlHistos.fillHisto("mjj",ctf,jetSystem.mass(),weight);
		      controlHistos.fillHisto("mt", ctf,mt,weight);
		      controlHistos.fillHisto("mindrlj",ctf,mindrlj,weight,true);
		      controlHistos.fillHisto("njets",ctf,nseljetsLoose-2,weight,true);
		      controlHistos.fillHisto("mljmostiso",ctf, mostIsolPair.mass(),weight);
		      controlHistos.fillHisto("mljleastiso",ctf,leastIsolPair.mass(),weight);
		      controlHistos.fillHisto("ptmostisollep",ctf,mostIsolLepton.pt(),weight);
		      controlHistos.fillHisto("ptleastisollep",ctf,leastIsolLepton.pt(),weight);
		      controlHistos.fillHisto("mtmostisolep",ctf,mostIsolMt,weight);
		      controlHistos.fillHisto("mtleastisolep",ctf,leastIsolMt,weight);
		      controlHistos.fillHisto("dphimetmostisollep",ctf,fabs(deltaPhi( mostIsolLepton.phi(),theMET.phi())),weight);
		      controlHistos.fillHisto("dphimetleastisolep",ctf,fabs(deltaPhi( leastIsolLepton.phi(),theMET.phi())),weight);
		      controlHistos.fillHisto("nleptons",ctf,phys.leptons.size()-2,weight);
		      controlHistos.fillHisto("alpha",ctf,min(nGoodassigments,2),weight);
		      controlHistos.fillHisto("leadjet",ctf,max(jetColl[0].pt(),jetColl[0].pt()),weight);
		      controlHistos.fillHisto("subleadjet",ctf,min(jetColl[0].pt(),jetColl[1].pt()),weight);
		      controlHistos.fill2DHisto("leadjetvssubleadjet",ctf,max(jetColl[0].pt(),jetColl[0].pt()),min(jetColl[0].pt(),jetColl[1].pt()),weight);
		      controlHistos.fillHisto("leadlepton",ctf,max(phys.leptons[0].pt(),phys.leptons[1].pt()),weight);
		      controlHistos.fillHisto("subleadlepton",ctf,min(phys.leptons[0].pt(),phys.leptons[1].pt()),weight);
		      controlHistos.fillHisto("drll",ctf,drll,weight);
		      controlHistos.fillHisto("dphill",ctf,dphill,weight);
		      for(size_t ijet=0; ijet<orderedJetColl.size(); ijet++)
			{
			  controlHistos.fillHisto("TCHE",ctf,orderedJetColl[ijet].btag1,puweight);
			  controlHistos.fillHisto("TCHP",ctf,orderedJetColl[ijet].btag2,puweight);
			  controlHistos.fillHisto("SSVHE",ctf,orderedJetColl[ijet].btag3,puweight);
			  controlHistos.fillHisto("JBP",ctf,orderedJetColl[ijet].btag4,puweight);
			  controlHistos.fillHisto("JP",ctf,orderedJetColl[ijet].btag5,puweight);
			  controlHistos.fillHisto("SSVHP",ctf,orderedJetColl[ijet].btag6,puweight);
			  controlHistos.fillHisto("CSV",ctf,orderedJetColl[ijet].btag7,puweight);
			}
		      for(size_t ilj=0; ilj<mljs.size(); ilj++)
			{
			  controlHistos.fillHisto("mlj",ctf,mljs[ilj],weight,true);
			  controlHistos.fillHisto("drlj",ctf,drljs[ilj],weight,true);
			}
		    }
		}
	      
	      if(isZcand) continue;
	      if(passMtCut) controlHistos.fillHisto("evtflow"+cats[ivar],ctf,SELMT,weight);
	      if(passTightLepton) controlHistos.fillHisto("evtflow"+cats[ivar],ctf,SELTIGHTLEPTONS,weight);
	      if(passTightJets) controlHistos.fillHisto("evtflow"+cats[ivar],ctf,SELTIGHTJETS,weight);
	      if(useMVA && useMVA)
		{
		  for(size_t imet=0; imet<methodList.size(); imet++)	      
		    {
		      controlHistos.fillHisto(methodList[imet],ctf,discriResults[imet],weight);
		      if(passMtCut)
			{
			  controlHistos.fillHisto(methodList[imet]+"Tight",ctf,discriResults[imet],weight);
			  controlHistos.fillHisto(methodList[imet]+btagbinLabel+"Tight",ctf,discriResults[imet],weight);
			}
		    }
		}
	      
	      //save summary if required
	      if(ictf==0 && ivar==0 && spyEvents && normWeight==1)
		{
		  std::vector<float> measurements;
		  ev.weight=summaryWeight;
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
  outDirs["qqlike"]=baseOutDir->mkdir("qqlike");
  for(SelectionMonitor::StepMonitor_t::iterator it =mons.begin(); it!= mons.end(); it++)
    {
      TString icat("all");
      if(it->first.BeginsWith("ee")) icat="ee";
      if(it->first.BeginsWith("emu")) icat="emu";
      if(it->first.BeginsWith("mumu")) icat="mumu";
      if(it->first.BeginsWith("qqlike")) icat="qqlike";
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
