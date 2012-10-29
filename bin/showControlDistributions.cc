#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuSummaryHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuPhysicsEvent.h"
#include "CMGTools/HtoZZ2l2nu/interface/MacroUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/METUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "CMGTools/HtoZZ2l2nu/interface/SmartSelectionMonitor.h"

#include "LIP/Top/interface/LeptonEfficiency.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "PhysicsTools/CondLiteIO/interface/RecordWriter.h"
#include "DataFormats/FWLite/interface/Record.h"
#include "DataFormats/FWLite/interface/EventSetup.h"
#include "DataFormats/FWLite/interface/ESHandle.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"


#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCanvas.h"
#include "TString.h"
#include "TDirectory.h"
#include "TEventList.h"
#include "TRandom.h"

#include <iostream>

using namespace std;

bool sortJetsByCSV(PhysicsObject_Jet a,PhysicsObject_Jet b) {   return (a.btag2>b.btag2); }
bool sortJetsByPt(PhysicsObject_Jet a,PhysicsObject_Jet b)  {   return (a.pt()>b.pt());   }

std::pair<float,float> getArcCos(LorentzVector &a, LorentzVector &b)
{
  TVector3 mom1(a.px(),a.py(),a.pz());
  TVector3 mom2(b.px(),b.py(),b.pz());
  double cosine = mom1.Dot(mom2)/(mom1.Mag()*mom2.Mag());
  double arcCosine = acos(cosine);
  return std::pair<float,float>(cosine,arcCosine);
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
  double sfMetCut = runProcess.getParameter<double>("sfMetCut");
  double ofMetCut = runProcess.getParameter<double>("ofMetCut");
  double jetPtCut = runProcess.getParameter<double>("jetPtCut");
  int jetIdToUse=JETID_LOOSE;
  int eIdToUse=EID_TIGHT;  float eIsoToUse=0.1;
  int mIdToUse=MID_TIGHT;  float mIsoToUse=0.2;//0.12;
  bool applyDYweight =  runProcess.getParameter<bool>("applyDYweight");
  TString uncFile =  runProcess.getParameter<std::string>("jesUncFileName");      gSystem->ExpandPathName(uncFile);

  //
  // check input file
  //
  TFile *evfile = TFile::Open(evurl);
  if(evfile==0) return -1;
  if(evfile->IsZombie()) return -1;
  TString proctag=gSystem->BaseName(evurl); 
  Ssiz_t pos=proctag.Index(".root");
  proctag.Remove(pos,proctag.Length());
  cout << "Processing: " << proctag << " @ " << evurl << endl;
  
  //
  // pileup reweighter
  //
  std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
  std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
  std::vector<float> mcPileupDistribution;
  if(isMC){
    //    TString puDist("evAnalyzer/h2zz/pileuptrue");
    TString puDist("evAnalyzer/h2zz/pileup");
    TH1F* histo = (TH1F *) evfile->Get(puDist);
    if(!histo)std::cout<<"pileup histogram is null!!!\n";
    for(int i=1;i<=histo->GetNbinsX();i++){mcPileupDistribution.push_back(histo->GetBinContent(i));}
    delete histo;
  }
  while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
  while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
  
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE                                                                                                                                          
  edm::LumiReWeighting *LumiWeights=0;
  PuShifter_t PuShifters;
  if(isMC)
    {
      LumiWeights= new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
      PuShifters=getPUshifters(dataPileupDistribution,0.05);
    }

  LeptonEfficiency lEff(2012);
  
  
  //
  // Instantiate uncertainty sources
  //
  TString systVars[]={"","jerdown","jerup","jesdown","jesup","pudown","puup"};//,"umetdown","umetup","leffup","leffdown"};
  size_t nSystVars=sizeof(systVars)/sizeof(TString); 
  if(!runSystematics || !isMC) { nSystVars=1; }
  else                         { cout << "Running sytematics: this will take a while..." << endl; }
  JetCorrectionUncertainty *totalUnc = new JetCorrectionUncertainty(uncFile.Data()); // *(new JetCorrectorParameters(uncFile.Data(), "Total")));

  //
  // control histograms
  //
  SmartSelectionMonitor controlHistos;
  TH1F* Hcutflow  = (TH1F*) controlHistos.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,5,0,5) ) ;
  TH1F* Hoptim_systs     =  (TH1F*) controlHistos.addHistogram( new TH1F ("optim_systs"    , ";syst;", nSystVars,0,nSystVars) );

  //vertex multiplicity
  controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 50, 0.,50.) );
  controlHistos.addHistogram( new TH1F ("nvertices_unwgt", "; Vertex multiplicity; Events", 50, 0.,50.) );
  
  ///lepton control
  controlHistos.addHistogram( new TH1D("dilmass",";M(l,l') [GeV];Events",50,0,250) );
  TH1D *lepMult=new TH1D("nleptons",";Leptons;Events",3,0,3);
  lepMult->GetXaxis()->SetBinLabel(1,"=2 leptons");
  lepMult->GetXaxis()->SetBinLabel(2,"=3 leptons");
  lepMult->GetXaxis()->SetBinLabel(3,"#geq 4 leptons");
  controlHistos.addHistogram( lepMult );
  TH1 *sslepMult = (TH1 *) lepMult->Clone("ssnleptons");
  controlHistos.addHistogram( sslepMult );
  controlHistos.addHistogram( new TH1F ("leadlepton", "; Leading lepton p_{T} [GeV]; Events / (10 GeV)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadlepton", "; Sub-leading lepton p_{T} [GeV]; Events / (10 GeV)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("leadleptoneta", "; Leading lepton #eta; Events", 25, 0.,2.5) );
  controlHistos.addHistogram( new TH1F ("subleadleptoneta", "; Sub-leading lepton #eta; Events", 25, 0.,2.5) );
  controlHistos.addHistogram( new TH1D("dilcharge",";Charge;Events",3,-1.5,1.5) );
  controlHistos.addHistogram( new TH1D("dphill",";#Delta#phi(l^{(1)},l^{(2)});Events",100,-3.2,3.2) );
  controlHistos.addHistogram( new TH1D("drll",";#Delta R(l^{(1)},l^{(2)});Events",100,0,6) );

  //jet control
  for(size_t ijet=1; ijet<=4; ijet++)
    {
      TString jetctr(""); jetctr += ijet;
      controlHistos.addHistogram( new TH1F ("jet"+jetctr, "; Jet #"+jetctr+" p_{T} [GeV]; Events / (10 GeV)", 50, 0.,500.) );
      controlHistos.addHistogram( new TH1F ("jet"+jetctr+"eta", "; Jet #"+jetctr+" #eta; Events", 30, 0.,3.) );
    }
  if(isMC)
    {
      TString jetFlavors[]={"b","c","udsg"};
      for(size_t iflav=0; iflav<3; iflav++)
	{
	  controlHistos.addHistogram( new TH1F ("jet"+jetFlavors[iflav], ";"+jetFlavors[iflav]+" jet p_{T} [GeV]; Events / (10 GeV)", 50, 0.,500.) );
	  controlHistos.addHistogram( new TH1F ("jet"+jetFlavors[iflav]+"eta", ";"+jetFlavors[iflav]+"jet #eta; Events", 30, 0.,3.) );
	}
      TH1F *h=new TH1F ("jetflavor", "; Jet flavor per event; Events", 10*3,0,10*3);
      for(int i=1; i<=h->GetXaxis()->GetNbins(); i++)
	{
	  TString label(""); label+=(i-1)%10;
	  h->GetXaxis()->SetBinLabel(i,label);
	}
      controlHistos.addHistogram( h );
    }
  
  ///////////////////////////////////
  // b-tagging                     //
  ///////////////////////////////////
  TString btagger[]  = {"csv", "jp", "tchp", "tche", "ivf"};
  float btaggerMin[] = {-0.2,  0.0,   -2,    -2,     -2 };
  float btaggerMax[] = {1.2,   2.5,   20,    20,     8};
  int idxForBtaggerToTemplate=1;
  TString btaggerToTemplate=btagger[idxForBtaggerToTemplate];
  TString btaggerWPs[]={"L","M","T"};
  TString jetRanges[]={"inc","30to50","50to80","80to120","120to200","200to300","300toInf","0to0.5","0.5to1.0","1.0to1.5","1.5to2.5","0to10","11to14","15to18","18toInf"};
  TString jetFlavors[]={"","b","udsg","c"};
  for(size_t ijf=0; ijf<4; ijf++)
    { 
      controlHistos.addHistogram( new TH1F (jetFlavors[ijf]+"jetlxy", "; SecVtx L_{xy} [cm]; Jets / (0.1 cm)", 50, 0.,5) );
      controlHistos.addHistogram( new TH1F (jetFlavors[ijf]+"jetlxysig", "; SecVtx #sigma_{L_{xy}}/ L_{xy}; Jets", 50, 0.,0.5) );
      controlHistos.addHistogram( new TH1F (jetFlavors[ijf]+"jetmass", "; SecVtx Mass [GeV]; Events / (0.2 GeV)", 50, 0.,10.) );
      controlHistos.addHistogram( new TH1F (jetFlavors[ijf]+"jetpt", "; Jet p_{T} [GeV]; Events / (10 GeV)", 50, 0.,500.) );
      controlHistos.addHistogram( new TH1F (jetFlavors[ijf]+"jeteta", "; Jet #eta; Events", 30, 0.,3.) );
      controlHistos.addHistogram( new TH2F (jetFlavors[ijf]+"jetptvseta", "; Jet p_{T} [GeV]; Jet #eta; Events / (10 GeV)", 50, 0.,500., 5, 0., 2.5) );
      controlHistos.addHistogram( new TH1F (jetFlavors[ijf]+"tagjetpt", "; Taggable jet p_{T} [GeV]; Events / (10 GeV)", 50, 0.,500.) );
      controlHistos.addHistogram( new TH1F (jetFlavors[ijf]+"tagjeteta", "; Taggable jet #eta; Events", 30, 0.,3.) );
    }
  TH1 *hflav=controlHistos.addHistogram( new TH1F("flavcount","Flavor at generator level",64,0,64) );
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      for(int k=0;k<4;k++)
	{
	  char buf[50]; sprintf(buf,"%db%dc%dl",i,j,k);
	  int bin(i | j<<2 | k<<4);
	  hflav->GetXaxis()->SetBinLabel(bin,buf);
	  cout << bin << " ";
	}
  cout << endl;
  const size_t nJetRanges=sizeof(jetRanges)/sizeof(TString);
  for(size_t i=0; i<sizeof(btagger)/sizeof(TString); i++)
    {
      TH1 *hl=controlHistos.addHistogram( new TH1F(btagger[i]+"Lbtagsextended",";b-tag multiplicity;Events", 3*3*5, 0.,3*3*5.) );
      TH1 *hm=controlHistos.addHistogram( new TH1F(btagger[i]+"Mbtagsextended",";b-tag multiplicity;Events", 3*3*5, 0.,3*3*5.) );
      TH1 *ht=controlHistos.addHistogram( new TH1F(btagger[i]+"Tbtagsextended",";b-tag multiplicity;Events", 3*3*5, 0.,3*3*5.) );
      for(int ibin=1; ibin<=hl->GetXaxis()->GetNbins(); ibin++)
	{
	  TString label(""); label += (ibin-1)%5;
	  hl->GetXaxis()->SetBinLabel(ibin,label);
	  hm->GetXaxis()->SetBinLabel(ibin,label);
	  ht->GetXaxis()->SetBinLabel(ibin,label);
	}
      
      for(size_t j=0; j<nSystVars; j++)
	{  
	  controlHistos.addHistogram( new TH1F("inc"+btagger[i]+systVars[j],";Discriminator;Jets", 50, btaggerMin[i],btaggerMax[i]) );

	  //use templates based on a given discriminator except in the case we're measuring the discriminator efficiency itself (switch to alternative)
	  int idxToUse=idxForBtaggerToTemplate;
	  TH2 *hinc  = (TH2 *)controlHistos.addHistogram( new TH2F(btagger[idxToUse]+systVars[j],";Discriminator;Jets", 50, btaggerMin[idxToUse],btaggerMax[idxToUse],nJetRanges,0,nJetRanges) );
	  TH2 *hb    = (TH2 *)controlHistos.addHistogram( new TH2F(btagger[idxToUse]+"b"+systVars[j],";b-tags;Range;Jets", 50,btaggerMin[idxToUse],btaggerMax[idxToUse],nJetRanges,0,nJetRanges) );
	  TH2 *hc    = (TH2 *)controlHistos.addHistogram( new TH2F(btagger[idxToUse]+"c"+systVars[j],";c-tags;Range;Jets", 50,btaggerMin[idxToUse],btaggerMax[idxToUse],nJetRanges,0,nJetRanges) );
	  TH2 *hudsg = (TH2 *)controlHistos.addHistogram( new TH2F(btagger[idxToUse]+"udsg"+systVars[j],";udsg-tags;Range;Jets", 50,btaggerMin[idxToUse],btaggerMax[idxToUse],nJetRanges,0,nJetRanges) );
	  for(int ybin=1; ybin<=hb->GetYaxis()->GetNbins(); ybin++)
	    {
	      hinc->GetYaxis()->SetBinLabel(ybin,jetRanges[ybin-1]);
	      hb->GetYaxis()->SetBinLabel(ybin,jetRanges[ybin-1]);
	      hc->GetYaxis()->SetBinLabel(ybin,jetRanges[ybin-1]);
	      hudsg->GetYaxis()->SetBinLabel(ybin,jetRanges[ybin-1]);
	    }
	  
	  for(size_t k=0; k<sizeof(btagger)/sizeof(TString); k++)
	    {
	      for(size_t iwp=0; iwp<sizeof(btaggerWPs)/sizeof(TString); iwp++)
		{
		  for(size_t ipf=0; ipf<2; ipf++)
		    {
		      TString pfix(ipf==0?"":"v");
		      controlHistos.addHistogram( (TH2 *)hinc->Clone( btagger[idxToUse]       +systVars[j]+btagger[k]+btaggerWPs[iwp]+pfix) );
		      controlHistos.addHistogram( (TH2 *)hb->Clone(   btagger[idxToUse]+"b"   +systVars[j]+btagger[k]+btaggerWPs[iwp]+pfix) );
		      controlHistos.addHistogram( (TH2 *)hc->Clone(   btagger[idxToUse]+"c"   +systVars[j]+btagger[k]+btaggerWPs[iwp]+pfix) );
		      controlHistos.addHistogram( (TH2 *)hudsg->Clone(btagger[idxToUse]+"udsg"+systVars[j]+btagger[k]+btaggerWPs[iwp]+pfix) );
		    }
		}
	    }
	}
    }
  
  //final control
  controlHistos.addHistogram( new TH1F ("pttbar", "; p_{T} (t#bar{t}) [GeV]; Events / (10 GeV)", 20, 0.,200.) );
  controlHistos.addHistogram( new TH1D("zmlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV];Lepton-jet pairs",50,0,500));
  controlHistos.addHistogram( new TH1D("mlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV];Lepton-jet pairs",50,0,500));
  controlHistos.addHistogram( new TH1D("correctmlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV];Lepton-jet pairs",50,0,500));
  controlHistos.addHistogram( new TH1D("wrongmlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV];Lepton-jet pairs",50,0,500));


  //event selection histogram
  enum SelectionSteps { SEL2LEPTONS, SELDILEPTON, SELJETS, SELMET, SELOS};
  TString labels[]={"2 leptons", "M>12 #wedge |M-M_{Z}|>15", "#geq 2 jets", "E_{T}^{miss}>40/0", "op. sign" };
  int nsteps=sizeof(labels)/sizeof(TString);

  //control histograms per variation
  TString thetallcats[]={"","eq0jets","eq1jets"};
  for(size_t ivar=0;ivar<nSystVars; ivar++) 
    {
      Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, systVars[ivar]);

      TH1D *cutflowH=new TH1D("evtflow"+systVars[ivar],";Cutflow;Events",nsteps,0,nsteps);
      for(int ibin=0; ibin<nsteps; ibin++) cutflowH->GetXaxis()->SetBinLabel(ibin+1,labels[ibin]);
      controlHistos.addHistogram( cutflowH );

      TH1D *finalCutflowH=new TH1D("finalevtflow"+systVars[ivar],";Category;Events",4,0,4);
      finalCutflowH->GetXaxis()->SetBinLabel(1,"=1 jets");
      finalCutflowH->GetXaxis()->SetBinLabel(2,"=2 jets");
      finalCutflowH->GetXaxis()->SetBinLabel(3,"=3 jets");
      finalCutflowH->GetXaxis()->SetBinLabel(4,"=4 jets");
      controlHistos.addHistogram( finalCutflowH );

      controlHistos.addHistogram( new TH1D("mtsum"+systVars[ivar],";M_{T}(l^{1},E_{T}^{miss})+M_{T}(l^{2},E_{T}^{miss});Events",100,0,1000) );
      controlHistos.addHistogram( new TH1D("eq1jetsmtsum"+systVars[ivar],";M_{T}(l^{1},E_{T}^{miss})+M_{T}(l^{2},E_{T}^{miss});Events",100,0,1000) );
      controlHistos.addHistogram( new TH1D("ptsum"+systVars[ivar],";p_{T}(l^{1})+p_{T}(l^{2});Events",100,0,500) );
      controlHistos.addHistogram( new TH1D("eq1jetsptsum"+systVars[ivar],";p_{T}(l^{1})+p_{T}(l^{2});Events",100,0,500) );
      controlHistos.addHistogram( new TH1D("mll"+systVars[ivar],";M(l,l') [GeV];Events",50,0,250) );
      controlHistos.addHistogram( new TH1D("eq1jetsmll",";M(l,l') [GeV];Events",50,0,250) );
      
      TH1F * h = new TH1F ("njets"+systVars[ivar], "; Jet multiplicity; Events", 6, 0.,6.);
      h->GetXaxis()->SetBinLabel(1,"=0 jets");
      h->GetXaxis()->SetBinLabel(2,"=1 jets");
      h->GetXaxis()->SetBinLabel(3,"=2 jets");
      h->GetXaxis()->SetBinLabel(4,"= 3 jets");
      h->GetXaxis()->SetBinLabel(5,"= 4 jets");
      h->GetXaxis()->SetBinLabel(6,"#geq 5 jets");
      controlHistos.addHistogram( h );
      controlHistos.addHistogram( (TH1F *)h->Clone("ntightjets") );
      controlHistos.addHistogram( (TH1F *)h->Clone("ncutbasedjets") );
      controlHistos.addHistogram( (TH1F *)h->Clone("prenjets") );
      controlHistos.addHistogram( (TH1F *)h->Clone("prentightjets") );
      controlHistos.addHistogram( (TH1F *)h->Clone("prencutbasedjets") );
       
      controlHistos.addHistogram( new TH1F ("met"+systVars[ivar], "; #slash{E}_{T} [GeV]; Events / (10 GeV)", 50, 0.,500.) );
      controlHistos.addHistogram( new TH1F ("rho"+systVars[ivar]  ,"; #rho [GeV]; Events / (10 GeV)", 50, 0., 25.) );
        
      for(size_t k=0; k<sizeof(thetallcats)/sizeof(TString); k++)
	{
	  controlHistos.addHistogram( new TH1D(thetallcats[k]+"dilarccosine"+systVars[ivar],";arcCos(l,l');Events",50,0,3.2) );
	  controlHistos.addHistogram( new TH1D(thetallcats[k]+"lowmetdilarccosine"+systVars[ivar],";arcCos(l,l');Events",50,0,3.2) );
	}
    }
  
  
  ///
  // process events file
  //
  DuplicatesChecker duplicatesChecker;
  ZZ2l2nuSummaryHandler evSummaryHandler;
  if( !evSummaryHandler.attachToTree( (TTree *)evfile->Get(dirname) ) )  { evfile->Close();  return -1; }  
  TTree *evTree=evSummaryHandler.getTree();

  int fType(0);
  if(evurl.Contains("DoubleEle")) fType=EE;
  if(evurl.Contains("DoubleMu"))  fType=MUMU;
  if(evurl.Contains("MuEG"))      fType=EMU;
  if(evurl.Contains("SingleMu"))  fType=MUMU;
  bool isSingleMuPD(!isMC && evurl.Contains("SingleMu"));
  bool isTauEmbedded(!isMC && evurl.Contains("DYTauEmbedded"));
  bool isDYJetsMC(isMC && evurl.Contains("DYJetsToLL"));

  //
  //total entries to process and normalization
  //
  const Int_t totalEntries=evTree->GetEntriesFast();
  if(evEnd<0 || evEnd>totalEntries) evEnd=totalEntries;
  if(evStart > evEnd || totalEntries==0)
    {
      evfile->Close();
      return -1;
    }
  float cnorm=1.0;
  if(isMC){
    TH1F* cutflowH = (TH1F *) evfile->Get("evAnalyzer/h2zz/cutflow");
    if(cutflowH) cnorm=cutflowH->GetBinContent(1);
    printf("cnorm = %f\n",cnorm);
  }
  Hcutflow->SetBinContent(1,cnorm);
  cout << " xSec x BR: " << xsec << endl
       << totalEntries << " events to be analyzed out of " << cnorm << " original events"<< endl;
  
  //prepare to save summaries
  ZZ2l2nuSummaryHandler *spyEvents=0;
  TFile *spyFile=0;
  TDirectory *spyDir=0;
  float summaryWeight(1);
  if(saveSummaryTree)
    {
      gSystem->Exec("mkdir -p " + outdir);
      TString summaryName = outdir + "/" + proctag;
      if(mcTruthMode!=0) { summaryName += "_filt"; summaryName += mcTruthMode; } 
      summaryName += "_summary.root";
      gSystem->ExpandPathName(summaryName);
      summaryWeight = xsec / cnorm;
      spyEvents = new ZZ2l2nuSummaryHandler;
      spyFile = TFile::Open(summaryName,"RECREATE");
      TString evtag=proctag;
      spyFile->rmdir(evtag);
      spyDir = spyFile->mkdir(evtag);
      TTree *outT = evTree->CloneTree(0);
      outT->SetTitle("Event summary");
      outT->SetDirectory(spyDir);
      spyEvents->initTree(outT,false);
      cout << "Creating event summary file:" << summaryName << endl;
    }
  
  //
  // analyze (puf...)
  //
  float selEvents(0);
  int NumberOfDuplicated(0);
  cout << evurl << " " << isMC << endl;
  for (int inum=evStart; inum < evEnd; ++inum)
    {
      if(inum%500==0) { printf("\r [ %d/100 ]",int(100*float(inum-evStart)/float(evEnd))); cout << flush; }
    
      evTree->GetEvent(inum);
      ZZ2l2nuSummary_t &ev = evSummaryHandler.getEvent();
      if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { NumberOfDuplicated++; continue; }
  
      //MC filter for specific processes
      int nleptons=(getNgenLeptons(ev.mccat,ELECTRON)+getNgenLeptons(ev.mccat,MUON)+getNgenLeptons(ev.mccat,TAU))/2;  //->bug: counted twice
      bool isTop=isTTbar(ev.mccat);
      bool isSTop=isSingleTop(ev.mccat);
      if(isMC && mcTruthMode==1 && (!isTop || nleptons<2) ) continue;
      if(isMC && mcTruthMode==2 && (isTop && nleptons>=2) ) continue;
      
      //require trigger to be consistent with the event to analyse
      if(ev.cat!=MUMU && ev.cat !=EE && ev.cat!=EMU) continue;
      bool hasEEtrigger = ev.triggerType & 0x1;
      bool hasMMtrigger = (ev.triggerType >> 1 ) & 0x1;
      bool hasEMtrigger = (ev.triggerType >> 2 ) & 0x1;
      bool hasMtrigger  = (ev.triggerType >> 3 ) & 0x1;
      if(!isMC && !isTauEmbedded){
	if(ev.cat!=fType) continue;
	if(ev.cat==EE   && !hasEEtrigger) continue;
	if(ev.cat==MUMU && !(hasMMtrigger||hasMtrigger) ) continue;
	if(ev.cat==EMU  && !hasEMtrigger) continue;

	//this is a safety veto for the single mu PD                                                                       
	if(isSingleMuPD) {
	  if(ev.cat!=MUMU) continue;
	  if(!hasMtrigger) continue;
	  if(hasMtrigger && hasMMtrigger) continue;
	}
      }
 
      //the channel to which the event belons
      bool isSameFlavor(false);
      TString ch("");
      if(ev.cat==MUMU)     { isSameFlavor=true; ch="mumu";  }
      else if(ev.cat==EE)  { isSameFlavor=true; ch="ee";    }
      else if(ev.cat==EMU) { isSameFlavor=false; ch="emu";  }
      
      //get the weights for this event
      float puweight=1;
      double TotalWeight_plus = 1.0;
      double TotalWeight_minus = 1.0;
      if(LumiWeights) 
	{
	  puweight = LumiWeights->weight( ev.ngenITpu );
	  TotalWeight_plus  = PuShifters[PUUP]->Eval(ev.ngenITpu);
	  TotalWeight_minus = PuShifters[PUDOWN]->Eval(ev.ngenITpu);
	}
      Hcutflow->Fill(1,1);
      Hcutflow->Fill(2,puweight);
      Hcutflow->Fill(3,puweight*TotalWeight_minus);
      Hcutflow->Fill(4,puweight*TotalWeight_plus);

      //
      // PHYSICS
      //
      PhysicsEvent_t phys=getPhysicsEventFrom(ev);
     
      //require ID and ISO for the dilepton
      bool passIdAndIso(true);
      float leffSF(1.0), leffSF_plus(1.0),leffSF_minus(1.0);
      for(size_t ilep=0; ilep<2; ilep++){
	int lpid=phys.leptons[ilep].pid;
	if(fabs(phys.leptons[ilep].id)==13){
	  float relIso = phys.leptons[ilep].pfRelIsoDbeta();
	  passIdAndIso &= hasObjectId(ev.mn_idbits[lpid], mIdToUse);
	  passIdAndIso &= (relIso<mIsoToUse);
	}else{
	  float relIso = phys.leptons[ilep].ePFRelIsoCorrected2012(ev.rho);
	  passIdAndIso &= hasObjectId(ev.en_idbits[lpid], eIdToUse);
	  passIdAndIso &= (relIso<eIsoToUse);
	}

	if(isMC)
	  {
	    std::pair<float,float> ieff=lEff.getLeptonEfficiency(phys.leptons[ilep].pt(),phys.leptons[ilep].eta(),fabs(phys.leptons[ilep].id));
	    leffSF       *= ieff.first;
	    leffSF_plus  *= (1.0+ieff.second/ieff.first);
	    leffSF_minus *= (1.0-ieff.second/ieff.first);
	  }
      }
      if(!passIdAndIso)continue;

      //count extra leptons
      int nextraleptons(0);
      PhysicsObjectLeptonCollection extraLeptons;
      for(size_t ilep=2; ilep<phys.leptons.size(); ilep++)
	{
	  bool isGood(false);
	  int lpid=phys.leptons[ilep].pid;
	  if(fabs(phys.leptons[ilep].id)==13)
	    {
	      isGood = (hasObjectId(ev.mn_idbits[lpid], MID_LOOSE) && phys.leptons[ilep].pfRelIsoDbeta()<0.2);
	      isGood |= (hasObjectId(ev.mn_idbits[lpid], MID_SOFT) && phys.leptons[ilep].pt()>3);
	    }
	  else
	    {
	      isGood = ( hasObjectId(ev.en_idbits[lpid],EID_VETO) && phys.leptons[ilep].ePFRelIsoCorrected2012(ev.rho)<0.15 && phys.leptons[ilep].pt()>10);
	    }
	  nextraleptons += isGood;
	  if(isGood) extraLeptons.push_back(phys.leptons[ilep]);
	}  

      //apply JER correction to jets (data only) and compute the variations
      std::vector<PhysicsObjectJetCollection> jetsVar;
      LorentzVectorCollection metsVar;
      //METUtils::computeVariation(phys.ajets, phys.leptons, phys.met[0], jetsVar, metsVar, totalUnc);   //AK5PF
      METUtils::computeVariation(phys.jets, phys.leptons, phys.met[0], jetsVar, metsVar, totalUnc);  //AK5PFchs

      //do the selection (taking into account the different varied collections)
      bool passBaseSelection(false);
      for(size_t ivar=0;ivar<(isMC ? nSystVars : 1); ivar++) 
	{
	  std::vector<TString> catsToFill;
	  catsToFill.push_back("all");
	  catsToFill.push_back(ch);

	  //weight to assign to event	
	  //float weight= (isMC ? puweight*leffSF : 1.0);
	  float weight= (isMC ? puweight : 1.0);

	  //trigger efficiencies
	  if(isMC)
	    {
	      if(ev.cat==EE)   weight *= 0.893;   //from AN-178/2012 (dilepton xsec)
	      if(ev.cat==MUMU) weight *= 0.927;  
	      if(ev.cat==EMU)  weight  *= 0.870; 
	    }
	  if(systVars[ivar]=="puup" )      { weight *= TotalWeight_plus;  }
	  if(systVars[ivar]=="pudown" )    { weight *= TotalWeight_minus; }
	  if(systVars[ivar]=="leffup" )    { weight *= leffSF_plus;  }
	  if(systVars[ivar]=="leffdown" )  { weight *= leffSF_minus; }

  	  //jets to use 
	  int ngoodJets(0),ngoodTightJets(0),ngoodCutBasedJets(0);
	  PhysicsObjectJetCollection prunedJetColl;
	  PhysicsObjectJetCollection jetColl = jetsVar[METUtils::JER];
	  LorentzVector theMET               = metsVar[METUtils::JER];
	  if(systVars[ivar]=="jerup")    { jetColl=jetsVar[METUtils::JER_UP];    theMET=metsVar[METUtils::JER_UP]; }
	  if(systVars[ivar]=="jerdown")  { jetColl=jetsVar[METUtils::JER_DOWN];  theMET=metsVar[METUtils::JER_DOWN]; }
	  if(systVars[ivar]=="jesdown")  { jetColl=jetsVar[METUtils::JES_DOWN];  theMET=metsVar[METUtils::JES_DOWN]; }
	  if(systVars[ivar]=="jesup" )   { jetColl=jetsVar[METUtils::JES_UP];    theMET=metsVar[METUtils::JES_UP]; }
	  if(systVars[ivar]=="umetdown") { jetColl=jetsVar[METUtils::UMET_DOWN]; theMET=metsVar[METUtils::UMET_DOWN]; }
	  if(systVars[ivar]=="umetup" )  { jetColl=jetsVar[METUtils::UMET_UP];   theMET=metsVar[METUtils::UMET_UP]; }
	  for(size_t ijet=0; ijet<jetColl.size(); ijet++) 
	    {
	      float pt=jetColl[ijet].pt();
	      float eta=fabs(jetColl[ijet].eta());
	      if(pt<=jetPtCut || fabs(eta)>=2.5) continue;
	      bool isGoodJet = hasObjectId(jetColl[ijet].pid,jetIdToUse);
	      if(!isGoodJet) continue;
	      ngoodTightJets += hasObjectId(jetColl[ijet].pid,JETID_MIN_LOOSE);
	      ngoodCutBasedJets += hasObjectId(jetColl[ijet].pid,JETID_CUTBASED_LOOSE);
	      ngoodJets++;
	      prunedJetColl.push_back( jetColl[ijet] );
	    }
	  sort(prunedJetColl.begin(),prunedJetColl.end(),sortJetsByPt);
	  PhysicsObjectJetCollection csvJetColl  = prunedJetColl; 
	  sort(csvJetColl.begin(),csvJetColl.end(),sortJetsByCSV);
	  TString jetCat          = "";
	  if(ngoodJets==0) jetCat = "eq0jets";
	  if(ngoodJets==1) jetCat = "eq1jets"; 

	  //lepton related quantities
	  PhysicsObject_Lepton l1       = (phys.leptons[0].pt() > phys.leptons[1].pt() ? phys.leptons[0] : phys.leptons[1]);
	  PhysicsObject_Lepton l2       = (phys.leptons[0].pt() < phys.leptons[1].pt() ? phys.leptons[0] : phys.leptons[1]);	  
	  float dilcharge              = (l1.id/fabs(l1.id))*(l2.id/fabs(l2.id));
	  LorentzVector dileptonSystem = l1+l2;
	  double acosine               = getArcCos(l1,l2).second;
	  double ptsum                 = l1.pt()+l2.pt();
	  double drll                  = deltaR(l1,l2);
	  double dphill                = deltaPhi(l1.phi(),l2.phi());
	  
	  //met related quantities
	  double leadlepmt    = METUtils::transverseMass(l1,theMET,false);
	  double subleadlepmt = METUtils::transverseMass(l2,theMET,false);
	  double mtsum        = leadlepmt+subleadlepmt;
 
	  //TTbar system candidate
	  LorentzVector ttbar_t = dileptonSystem+theMET;
	  if(csvJetColl.size()>=2) ttbar_t += csvJetColl[0]+csvJetColl[1];
	  double mt           = ttbar_t.mass();

	  //start filling selection histograms
	  bool isInQuarkoniaRegion( dileptonSystem.mass()<12 );
	  bool isZcand(isSameFlavor && fabs(dileptonSystem.mass()-91)<15);
	  //  bool isEmuInZRegion(!isSameFlavor  && fabs(dileptonSystem.mass()-91)<15);
	  bool passJet(ngoodJets>=2);
	  bool isOS(dilcharge<0);
	  bool passMet( (!isSameFlavor && theMET.pt()>ofMetCut) || (isSameFlavor && theMET.pt()>sfMetCut) );
	  if( isDYJetsMC && passMet && applyDYweight)
	    {
	      if(ngoodJets==1 && ev.cat==EE)   weight *= 1.58;
	      if(ngoodJets==1 && ev.cat==MUMU) weight *= 1.59;
	      if(ngoodJets==1 && ev.cat==EMU)  weight *= 1.04;
	      if(ngoodJets>1  && ev.cat==EE)   weight *= 1.69;
	      if(ngoodJets>1  && ev.cat==MUMU) weight *= 1.72;
	      if(ngoodJets>1  && ev.cat==EMU)  weight *= 1.13;
	    }
	  TString metCat("none");
	  if(passMet)              metCat="";
	  else if(theMET.pt()<=30) metCat="lowmet";

	  controlHistos.fillHisto("evtflow"+systVars[ivar],catsToFill,SEL2LEPTONS,weight);
	  if(ivar==0) 
	    {
	      controlHistos.fillHisto("nvertices",catsToFill,ev.nvtx,weight,true);
	      controlHistos.fillHisto("nvertices_unwgt",catsToFill,ev.nvtx,1,true);
	      passBaseSelection=(!isInQuarkoniaRegion && passJet && passMet && isOS);  
	    }
	  if(isInQuarkoniaRegion) continue;
	  if(!isZcand)   
	    {
	      controlHistos.fillHisto("evtflow"+systVars[ivar],catsToFill,SELDILEPTON,weight);
	      if(isOS)
		{
		  controlHistos.fillHisto(jetCat+metCat+"dilarccosine"+systVars[ivar],catsToFill,acosine,weight);
		  controlHistos.fillHisto("prenjets"+systVars[ivar],catsToFill,ngoodJets,weight);			      
		  controlHistos.fillHisto("prentightjets"+systVars[ivar],catsToFill,ngoodTightJets,weight);			      
		  controlHistos.fillHisto("prencutbasedjets"+systVars[ivar],catsToFill,ngoodCutBasedJets,weight);			      
		  if(passMet) 
		    {
		      controlHistos.fillHisto("njets"+systVars[ivar],catsToFill,ngoodJets,weight);			      
		      controlHistos.fillHisto("ntightjets"+systVars[ivar],catsToFill,ngoodTightJets,weight);			      
		      controlHistos.fillHisto("ncutbasedjets"+systVars[ivar],catsToFill,ngoodCutBasedJets,weight);			      
		    }
		}
	    }
	  if(!passJet)
	    {
	      if(!isZcand && passMet && isOS && ngoodJets==1)
		{
		  controlHistos.fillHisto("eq1jetsmtsum"+systVars[ivar],catsToFill,mtsum,weight);
		  controlHistos.fillHisto("eq1jetsmll"+systVars[ivar],catsToFill,dileptonSystem.mass(),weight);
		  controlHistos.fillHisto("e1jetsptsum"+systVars[ivar],catsToFill,ptsum,weight);
		  controlHistos.fillHisto("finalevtflow"+systVars[ivar],catsToFill,0,weight);
		}
	      continue;
	    }
	  if(!isZcand)
	    {
	      controlHistos.fillHisto("evtflow"+systVars[ivar],catsToFill,SELJETS,weight);
	      controlHistos.fillHisto("met"+systVars[ivar],catsToFill,theMET.pt(),weight);
	      controlHistos.fillHisto("rho"+systVars[ivar],catsToFill,ev.rho,weight);
	    }
	  if(!passMet) continue;
	  if(!isZcand)
	    {
	      controlHistos.fillHisto("evtflow"+systVars[ivar],catsToFill,SELMET,weight);
	      if(ivar==0) controlHistos.fillHisto("dilcharge",catsToFill,dilcharge,weight);
	    }
 	  if(!isOS) continue;
	  if(ivar==0) controlHistos.fillHisto("dilmass",catsToFill,dileptonSystem.mass(),weight);
	  if(isZcand) continue;

	  //this is the final selected dilepton sample
	  if(ngoodJets<=4)  controlHistos.fillHisto("finalevtflow"+systVars[ivar],catsToFill,ngoodJets-1,weight);

	  controlHistos.fillHisto("evtflow"+systVars[ivar],catsToFill,SELOS,weight);
	  controlHistos.fillHisto("mtsum"+systVars[ivar],catsToFill,mtsum,weight);
	  controlHistos.fillHisto("mll"+systVars[ivar],catsToFill,dileptonSystem.mass(),weight);
	  controlHistos.fillHisto("ptsum"+systVars[ivar],catsToFill,ptsum,weight);
	  if(ivar==0)
	    {
	      selEvents+=weight;
	      controlHistos.fillHisto("leadlepton",catsToFill,l1.pt(),weight);
	      controlHistos.fillHisto("subleadlepton",catsToFill,l2.pt(),weight);
	      controlHistos.fillHisto("leadleptoneta",catsToFill,fabs(l1.eta()),weight);
	      controlHistos.fillHisto("subleadleptoneta",catsToFill,fabs(l2.eta()),weight);
	      controlHistos.fillHisto("drll",catsToFill,drll,weight);
	      controlHistos.fillHisto("dphill",catsToFill,dphill,weight);
	      controlHistos.fillHisto("pttbar",catsToFill,ttbar_t.pt(),weight);
	      controlHistos.fillHisto("mt", catsToFill,mt,weight);
	      controlHistos.fillHisto("nleptons",catsToFill,extraLeptons.size(),weight);
	    }
	     
	  //jet kinematics, flavor templates, etc.
	  int nbs(0),ncs(0),nudsg(0);
	  std::map<TString,int> btagsCount; 
	  btagsCount["tcheL"]=0;   btagsCount["tcheM"]=0;
	  btagsCount["csvL"]=0;    btagsCount["csvM"]=0;   btagsCount["csvT"]=0;
	  btagsCount["jpL"]=0;     btagsCount["jpM"]=0;    btagsCount["jpT"]=0;
	  btagsCount["tchpT"]=0;   btagsCount["ivfL"]=0;   btagsCount["ivfM"]=0;
	  btagsCount["tcheL"]=0;   btagsCount["tcheM"]=0; 
	  TString leadingLxyJetFlav("");
	  float leadingLxy(-1), leadingLxySig(-1), leadingLxyMass(-1);
	  for(size_t ijet=0; ijet<prunedJetColl.size(); ijet++) 
	    {
	      TString jetctr(""); jetctr+=ijet;
	      
	      //jet flavor
	      TString jetFlav("udsg");
	      if(fabs(prunedJetColl[ijet].flavid)==5)      { nbs++;   jetFlav="b"; }
	      else if(fabs(prunedJetColl[ijet].flavid)==4) { ncs++;   jetFlav="c"; }
	      else                                         { nudsg++; jetFlav="udsg"; }

	      //kinematics
	      float pt  = prunedJetColl[ijet].pt();
	      float eta = fabs(prunedJetColl[ijet].eta());
	      std::vector<int> jetCategs;
	      jetCategs.push_back(0);
	      if(pt>30 && pt<=50)           jetCategs.push_back(1);
	      if(pt>50 && pt<=80)           jetCategs.push_back(2);
	      if(pt>80 && pt<=120)          jetCategs.push_back(3);
	      if(pt>120 && pt<=200)         jetCategs.push_back(4);
	      if(pt>200 && pt<=300)         jetCategs.push_back(5);
	      if(pt>300)                    jetCategs.push_back(6);
	      if(eta<=0.5)                  jetCategs.push_back(7);
	      if(eta>0.5 && eta<=1.0)       jetCategs.push_back(8);
	      if(eta>1.0 && eta<=1.5)       jetCategs.push_back(9);
	      if(eta>1.5)                   jetCategs.push_back(10);
	      if(ev.nvtx<=10)               jetCategs.push_back(11);
	      if(ev.nvtx>10 && ev.nvtx<=14) jetCategs.push_back(12);
	      if(ev.nvtx>14 && ev.nvtx<=18) jetCategs.push_back(13);
	      if(ev.nvtx>19)                jetCategs.push_back(14);

	      //analyse the leading lxy 
	      if(prunedJetColl[ijet].lxy>leadingLxy)
		{
		  leadingLxy=prunedJetColl[ijet].lxy;
		  leadingLxySig=prunedJetColl[ijet].lxyerr/leadingLxy;
		  leadingLxyMass=prunedJetColl[ijet].svmass;
		  leadingLxyJetFlav=jetFlav;
		}

	      //control distributions
	      if(ivar==0)
		{
		  controlHistos.fillHisto("jet"+jetctr,catsToFill,pt,weight);
		  controlHistos.fillHisto("jet"+jetctr+"eta",catsToFill,fabs(eta),weight);
		  controlHistos.fillHisto( "jetpt", catsToFill,pt,weight);
		  controlHistos.fillHisto( "jeteta",catsToFill,fabs(eta),weight);
	       
		  if(isMC)
		    {
		      controlHistos.fillHisto( jetFlav+"jetpt", catsToFill,pt,weight);
		      controlHistos.fillHisto( jetFlav+"jetptvseta", catsToFill,pt,fabs(eta),weight);
		      controlHistos.fillHisto( jetFlav+"jeteta",catsToFill,fabs(eta),weight);
		    }
		  if(prunedJetColl[ijet].btag3>0)
		    {
		      controlHistos.fillHisto( "tagjetpt", catsToFill,pt,weight);
		      controlHistos.fillHisto( "tagjeteta", catsToFill,fabs(eta),weight);
		      if(isMC)
			{
			  controlHistos.fillHisto( jetFlav+"tagjetpt", catsToFill,pt,weight);
			  controlHistos.fillHisto( jetFlav+"tagjeteta", catsToFill,fabs(eta),weight);
			}
		    }
		}

	      //Lxy analysis
	      if(ivar==0 && leadingLxy>0)
		{
		  controlHistos.fillHisto("jetlxy",catsToFill,leadingLxy,weight);
		  controlHistos.fillHisto("jetlxysig",catsToFill,leadingLxySig,weight);
		  controlHistos.fillHisto("jetmass",catsToFill,leadingLxyMass,weight);
		  if(isMC)
		    {
		      controlHistos.fillHisto(leadingLxyJetFlav+"jetlxy",catsToFill,leadingLxy,weight);
		      controlHistos.fillHisto(leadingLxyJetFlav+"jetmass",catsToFill,leadingLxyMass,weight);
		    }
		}

	      //count flavors
	      if(ivar==0)
		{
		  std::vector<LorentzVector> genP4;
		  int ngenb(0),ngenc(0),ngenl(0);
		  for(int ngen=0; ngen<ev.nmcparticles; ngen++)
		    {
		      if(abs(ev.mc_id[ngen])>5) continue;
		      LorentzVector p4(ev.mc_px[ngen],ev.mc_py[ngen],ev.mc_pz[ngen],ev.mc_en[ngen]);
		      if(fabs(p4.eta())>2.5 || p4.pt()<25) continue;
		      
		      bool hasOverlap(false);
		      for(size_t igp4=0; igp4<genP4.size(); igp4++)
			{
			  float dR = deltaR(p4,genP4[igp4]);
			  if(dR>0.1) continue;
			  hasOverlap=true;
			}
		      if(hasOverlap) continue;
		      
		      if(abs(ev.mc_id[ngen])==5)      ngenb++;
		      else if(abs(ev.mc_id[ngen])==4) ngenc++;
		      else                         ngenl++;
		    }
		  ngenb=min(4,ngenb);
		  ngenc=min(4,ngenc);
		  ngenl=min(4,ngenl);
		  int bin(ngenb | ngenc<<2 | ngenl<<4);
		  controlHistos.fillHisto("flavcount", catsToFill, hflav->GetBinCenter(bin),  weight);
		}
	      
	      //b-tagging
	      std::map<TString, bool> hasTagger;
	      float tche = prunedJetColl[ijet].btag1;
	      bool hasTCHEL(tche>1.7);  btagsCount["tcheL"] += hasTCHEL; hasTagger["tcheL"]=hasTCHEL;
	      bool hasTCHEM(tche>1.93); btagsCount["tcheM"] += hasTCHEM; hasTagger["tcheM"]=hasTCHEM; 
	      float csv  = prunedJetColl[ijet].btag2;
	      bool hasCSVL(csv>0.244);  btagsCount["csvL"] += hasCSVL;   hasTagger["csvL"]=hasCSVL;
	      bool hasCSVM(csv>0.679);  btagsCount["csvM"] += hasCSVM;   hasTagger["csvM"]=hasCSVM;
	      bool hasCSVT(csv>0.898);  btagsCount["csvT"] += hasCSVT;   hasTagger["csvT"]=hasCSVT;
	      float jp   = prunedJetColl[ijet].btag3;
	      bool hasJPL(jp>0.275);    btagsCount["jpL"] += hasJPL;     hasTagger["jpL"]=hasJPL; 
	      bool hasJPM(jp>0.545);    btagsCount["jpM"] += hasJPM;     hasTagger["jpM"]=hasJPM;
	      bool hasJPT(jp>0.790);    btagsCount["jpT"] += hasJPT;     hasTagger["jpT"]=hasJPT;
	      float tchp = prunedJetColl[ijet].btag4;
	      bool hasTCHPT(tchp>3.41); btagsCount["tchpT"] += hasTCHPT; hasTagger["tchpT"]=hasTCHPT;
	      float ivfHighEff = prunedJetColl[ijet].btag5;
	      float ivfHighPur = prunedJetColl[ijet].btag6;
	      bool hasIVFHighEff(ivfHighEff>0);       btagsCount["ivfL"] += hasIVFHighEff;
	      bool hasIVFHighPur(ivfHighPur>0);       btagsCount["ivfM"] += hasIVFHighPur;

	      controlHistos.fillHisto("inctche"+systVars[ivar], catsToFill, tche,  weight);
	      controlHistos.fillHisto("inctchp"+systVars[ivar], catsToFill, tchp,  weight);
	      controlHistos.fillHisto("inccsv"+systVars[ivar],  catsToFill, csv,   weight);
	      controlHistos.fillHisto("incjp"+systVars[ivar],   catsToFill, jp,    weight);
	      controlHistos.fillHisto("incivf"+systVars[ivar],   catsToFill, ivfHighEff,    weight);
	      float btaggerToTemplateVal=csv;
	      if(btaggerToTemplate=="jp")   btaggerToTemplateVal=jp;
	      if(btaggerToTemplate=="tche") btaggerToTemplateVal=tche;
	      if(btaggerToTemplate=="tchp") btaggerToTemplateVal=tchp;
	      if(btaggerToTemplate=="ivf")  btaggerToTemplateVal=ivfHighEff;
	      for(size_t ijcat=0; ijcat<jetCategs.size(); ijcat++)
		{
		  controlHistos.fillHisto(btaggerToTemplate+systVars[ivar],      catsToFill, btaggerToTemplateVal,      jetCategs[ijcat], weight);
		  if(isMC) controlHistos.fillHisto(btaggerToTemplate+jetFlav+systVars[ivar],      catsToFill, btaggerToTemplateVal,      jetCategs[ijcat], weight);
		  for(std::map<TString, bool>::iterator it=hasTagger.begin(); it!=hasTagger.end(); it++)
		    {
		      TString pfix(it->second?"":"v");
		      controlHistos.fillHisto(btaggerToTemplate        +systVars[ivar]+it->first+pfix,catsToFill, btaggerToTemplateVal, jetCategs[ijcat], weight);
		      if(isMC) controlHistos.fillHisto( btaggerToTemplate+jetFlav+systVars[ivar]+it->first+pfix,     catsToFill, btaggerToTemplateVal, jetCategs[ijcat], weight);
		    }
		}

	      //lepton-jet kinematics
	      LorentzVector l1j=l1+prunedJetColl[ijet];
	      LorentzVector l2j=l2+prunedJetColl[ijet];
	      if(ivar==0)
		{
		  controlHistos.fillHisto("mlj",catsToFill,l1j.mass(),weight);
		  controlHistos.fillHisto("mlj",catsToFill,l2j.mass(),weight);
		  if(isMC) 
		    {
		      bool isL1JCorrect( isMC && (isTop || isSTop) && l1.genid*prunedJetColl[ijet].genid<0 && jetFlav=="b");
		      bool isL2JCorrect( isMC && (isTop || isSTop) && l2.genid*prunedJetColl[ijet].genid<0 && jetFlav=="b");
		      controlHistos.fillHisto(isL1JCorrect ? "correctmlj" : "wrongmlj" ,catsToFill,l1j.mass(),weight);
		      controlHistos.fillHisto(isL2JCorrect ? "correctmlj" : "wrongmlj" ,catsToFill,l2j.mass(),weight);
		      
		      controlHistos.fillHisto("jet"+jetFlav,catsToFill,pt,weight);
		      controlHistos.fillHisto("jet"+jetFlav+"eta",catsToFill,eta,weight);
		    }
		}
	    }
	  if(isMC && ivar==0)
	    {
	      controlHistos.fillHisto("jetflavor",catsToFill,nudsg,weight);		      
	      controlHistos.fillHisto("jetflavor",catsToFill,ncs+10,weight);		      
	      controlHistos.fillHisto("jetflavor",catsToFill,nbs+20,weight);		      
	    } 
	  
	  //b-tag counting
	  if(ngoodJets<=4)
	    {
	      int addBin(0);
	      if(ngoodJets==3) addBin += 5;
	      if(ngoodJets==4) addBin += 10;
	      if(ev.cat==MUMU) addBin += 15;
	      if(ev.cat==EMU)  addBin += 2*15;
	      for(std::map<TString,int>::iterator it=btagsCount.begin(); it!= btagsCount.end(); it++)
		{ 
		  //if(it->second==0) continue;
		  controlHistos.fillHisto(it->first+"btagsextended"+systVars[ivar],catsToFill,it->second+addBin,weight);
		}
	    }
	}

      //save events which pass base selection (also under the Z mass window)
      if(spyEvents && passBaseSelection)
	{
	  //sample the MC according to the PU in data (generate unweighted sample for pseudo-exp.)
	  //  int normWeight=0;
	  // 	      if(isMC && maxPuWeight>0)
	  // 		{
	  // 		  float rnd=gRandom->Uniform();
	  // 		  if( rnd > puweight/maxPuWeight) normWeight=1;
	  // 		    }
	  ev.hptWeights[0]=summaryWeight;
	  ev.hptWeights[1]=puweight;
	  ev.hptWeights[2]=puweight*TotalWeight_plus;
	  ev.hptWeights[3]=puweight*TotalWeight_minus;
	
	  std::vector<float> measurements;
	  //FIXME: add measurements
	  spyEvents->fillTree();// ev, measurements );
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
  outUrl += proctag;
  if(mcTruthMode!=0) { outUrl += "_filt"; outUrl += mcTruthMode; }
  outUrl += ".root";
  TFile *file=TFile::Open(outUrl, "recreate");
  controlHistos.Write();
  file->Close();
  
  //that's all folks!
}  
