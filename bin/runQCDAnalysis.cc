#include "CMGTools/HtoZZ2l2nu/interface/SmartSelectionMonitor.h"
#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuPhysicsEvent.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h" 

#include <fstream>
#include <iostream>

using namespace std;

struct QCDevent_t 
{
  int Run, nPU, nPUtrue;
  int nPV, nJet, nMuon, nBFromGSplit, nTrkInc, BitTrigger;
  int  Jet_flavour[50], Jet_nFirstTrkInc[50], Jet_nLastTrkInc[50], Muon_IdxJet[50]; 
  float pthat, PVz;
  float Jet_pt[50], Jet_eta[50], Jet_phi[50];                //JET kinematics
  float Jet_Ip1P[50], Jet_Ip2P[50], Jet_Ip3P[50];            //leading IP sig, TCHE, TCHP
  float Jet_Svx[50], Jet_SvxHP[50];                          //SSVHE, SSVHP
  float Jet_ProbaP[50], Jet_ProbaN[50], Jet_Proba[50];       //JP, JPn, JP
  float Jet_BprobP[50], Jet_BprobN[50], Jet_Bprob[50];       //JBP, JPBn, JBP
  float Jet_CombSvxP[50], Jet_CombSvxN[50], Jet_CombSvx[50]; //CSV, CSVn, CSV
  float Muon_pt[50], Muon_eta[50], Muon_ptrel[50], Muon_IP[50], Muon_IPsig[50], Muon_Proba[50], Muon_chi2[50], Muon_chi2Tk[50]; 
  int Muon_nMuHit[50], Muon_nTkHit[50], Muon_nPixHit[50], Muon_nOutHit[50], Muon_isGlobal[50], Muon_nMatched[50];
  float TrkInc_pt[200], TrkInc_ptrel[200], TrkInc_IP[200], TrkInc_IPsig[200];
  float bFromGSplit_pT[50], bFromGSplit_eta[50], bFromGSplit_phi[50]; //to do the matching
};

QCDevent_t qcdev;

TTree *getQCDTreeFrom(TFile *inF);
PhysicsObjectJetCollection getQCDJetCollection();


//
TTree *getQCDTreeFrom(TFile *inF)
{
  if(inF==0) return 0;
  TTree *tchain = (TTree*)inF->Get("mistag/ttree"); 
  if(tchain==0) return 0;

  //attach to tree
  tchain->SetBranchAddress("Run", &qcdev.Run);
  tchain->SetBranchAddress("BitTrigger", &qcdev.BitTrigger);
  tchain->SetBranchAddress("nPU", &qcdev.nPU);
  tchain->SetBranchAddress("nPUtrue", &qcdev.nPUtrue);
  tchain->SetBranchAddress("nPV", &qcdev.nPV);
  tchain->SetBranchAddress("pthat", &qcdev.pthat);
  tchain->SetBranchAddress("PVz", &qcdev.PVz);
  tchain->SetBranchAddress("nJet", &qcdev.nJet);
  tchain->SetBranchAddress("Jet_pt", qcdev.Jet_pt);
  tchain->SetBranchAddress("Jet_eta", qcdev.Jet_eta);
  tchain->SetBranchAddress("Jet_phi", qcdev.Jet_phi);
  tchain->SetBranchAddress("Jet_flavour", qcdev.Jet_flavour);
  tchain->SetBranchAddress("Jet_Ip1P", qcdev.Jet_Ip1P);
  tchain->SetBranchAddress("Jet_Ip2P", qcdev.Jet_Ip2P);
  tchain->SetBranchAddress("Jet_Ip3P", qcdev.Jet_Ip3P);
  tchain->SetBranchAddress("Jet_Svx", qcdev.Jet_Svx);
  tchain->SetBranchAddress("Jet_SvxHP", qcdev.Jet_SvxHP);
  tchain->SetBranchAddress("Jet_ProbaP", qcdev.Jet_ProbaP);
  tchain->SetBranchAddress("Jet_ProbaN", qcdev.Jet_ProbaN);
  tchain->SetBranchAddress("Jet_Proba", qcdev.Jet_Proba);
  tchain->SetBranchAddress("Jet_BprobP", qcdev.Jet_BprobP);
  tchain->SetBranchAddress("Jet_BprobN", qcdev.Jet_BprobN);
  tchain->SetBranchAddress("Jet_Bprob", qcdev.Jet_Bprob);
  tchain->SetBranchAddress("Jet_CombSvxP", qcdev.Jet_CombSvxP);
  tchain->SetBranchAddress("Jet_CombSvxN", qcdev.Jet_CombSvxN);
  tchain->SetBranchAddress("Jet_CombSvx", qcdev.Jet_CombSvx);
  tchain->SetBranchAddress("Jet_nFirstTrkInc", qcdev.Jet_nFirstTrkInc);
  tchain->SetBranchAddress("Jet_nLastTrkInc", qcdev.Jet_nLastTrkInc);
  tchain->SetBranchAddress("nMuon", &qcdev.nMuon);
  tchain->SetBranchAddress("Muon_IdxJet", qcdev.Muon_IdxJet);
  tchain->SetBranchAddress("Muon_pt", qcdev.Muon_pt);
  tchain->SetBranchAddress("Muon_eta", qcdev.Muon_eta);
  tchain->SetBranchAddress("Muon_ptrel", qcdev.Muon_ptrel);
  tchain->SetBranchAddress("Muon_IP", qcdev.Muon_IP);
  tchain->SetBranchAddress("Muon_IPsig", qcdev.Muon_IPsig);
  tchain->SetBranchAddress("Muon_Proba", qcdev.Muon_Proba);
  tchain->SetBranchAddress("Muon_isGlobal", qcdev.Muon_isGlobal);
  tchain->SetBranchAddress("Muon_nMatched", qcdev.Muon_nMatched);
  tchain->SetBranchAddress("Muon_nTkHit", qcdev.Muon_nTkHit);
  tchain->SetBranchAddress("Muon_nPixHit", qcdev.Muon_nPixHit);
  tchain->SetBranchAddress("Muon_nMuHit", qcdev.Muon_nMuHit);
  tchain->SetBranchAddress("Muon_nOutHit", qcdev.Muon_nOutHit);
  tchain->SetBranchAddress("Muon_chi2", qcdev.Muon_chi2);
  tchain->SetBranchAddress("Muon_chi2Tk", qcdev.Muon_chi2Tk);
  tchain->SetBranchAddress("nBFromGSplit", &qcdev.nBFromGSplit);
  tchain->SetBranchAddress("bFromGSplit_pT", qcdev.bFromGSplit_pT);
  tchain->SetBranchAddress("bFromGSplit_eta", qcdev.bFromGSplit_eta);
  tchain->SetBranchAddress("bFromGSplit_phi", qcdev.bFromGSplit_phi);
  tchain->SetBranchAddress("nTrkInc", &qcdev.nTrkInc);
  tchain->SetBranchAddress("TrkInc_pt", qcdev.TrkInc_pt);
  tchain->SetBranchAddress("TrkInc_ptrel", qcdev.TrkInc_ptrel);
  tchain->SetBranchAddress("TrkInc_IP", qcdev.TrkInc_IP);
  tchain->SetBranchAddress("TrkInc_IPsig", qcdev.TrkInc_IPsig);
  
  return tchain;
}

//
PhysicsObjectJetCollection getQCDJetCollection(float minJetPt, bool isMC)
{
  PhysicsObjectJetCollection jets;
  for(int ijet=0; ijet<qcdev.nJet; ijet++) 
    {
      //check the kinematics
      float pt=qcdev.Jet_pt[ijet];
      float eta=qcdev.Jet_eta[ijet];
      float phi=qcdev.Jet_phi[ijet];
      if(pt<minJetPt || fabs(eta)>2.5) continue;
      
      //check for consistency with the generated ptHat? (check this with Luca) 
      bool hasGoodKinematics(true);
      if(isMC)
	{
	  if (qcdev.pthat<  30. && qcdev.Jet_pt[ijet]>  85.) hasGoodKinematics = false;
	  if (qcdev.pthat<  50. && qcdev.Jet_pt[ijet]> 100.) hasGoodKinematics = false;
	  if (qcdev.pthat<  80. && qcdev.Jet_pt[ijet]> 150.) hasGoodKinematics = false;
	  if (qcdev.pthat< 120. && qcdev.Jet_pt[ijet]> 250.) hasGoodKinematics = false;
	  if (qcdev.pthat< 170. && qcdev.Jet_pt[ijet]> 320.) hasGoodKinematics = false;
	  if (qcdev.pthat< 470. && qcdev.Jet_pt[ijet]> 620.) hasGoodKinematics = false;
	  if (qcdev.pthat< 600. && qcdev.Jet_pt[ijet]> 720.) hasGoodKinematics = false;
	  if (qcdev.pthat< 800. && qcdev.Jet_pt[ijet]> 920.) hasGoodKinematics = false;
	  if (qcdev.pthat> 800. && qcdev.Jet_pt[ijet]< 600.) hasGoodKinematics = false;
	}
      if (!hasGoodKinematics) continue;
      
      //wrap the jet info
      float en=pt*cosh(eta);
      LorentzVector p4(pt*TMath::Cos(phi), pt*TMath::Sin(phi), pt*sinh(eta), en);
      PhysicsObject_Jet pjet(p4,0,0,0,0);
      pjet.setBtagInfo(qcdev.Jet_Ip2P[ijet],qcdev.Jet_CombSvxP[ijet],qcdev.Jet_ProbaP[ijet],qcdev.Jet_Ip3P[ijet],0,0,0);
      pjet.setGenInfo(qcdev.Jet_flavour[ijet],0,qcdev.Jet_flavour[ijet],false,false);
      jets.push_back(pjet);
    }

  //order by pt and return
  sort(jets.begin(),jets.end(),&JetPtOrdering);
  return jets;
}



bool PassTriggerBit(int ThisBit, int ThisCode) {

  if (ThisCode==-1) return true;

  int ThisMagnitude = 0;
  for (int mg = 0; mg<10; mg++) 
    if (ThisCode/pow(10, mg)>=1.) ThisMagnitude = mg;

  int ThisBitRounded = ThisBit/pow(10, ThisMagnitude+1);
  int ThisBitCut = ThisBit/pow(10, ThisMagnitude);
  int ThisBitCode = ThisBitCut - 10*ThisBitRounded;

  int ThisCodeCode = ThisCode/pow(10, ThisMagnitude);

  return (ThisBitCode & ThisCodeCode);
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
  const edm::ParameterSet &runProcess       = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  TString evurl                             = runProcess.getParameter<std::string>("input");
  TString outdir                            = runProcess.getParameter<std::string>("outdir");
  bool isMC                                 = runProcess.getParameter<bool>("isMC");
  TString dirname                           = runProcess.getParameter<std::string>("dirName");
  std::vector<std::string> jetPtWeightsFile = runProcess.getParameter<std::vector<std::string> >("weightsFile");
  double jetPtCut                           = runProcess.getParameter<double>("jetPtCut");
  TString uncFile                           = runProcess.getParameter<std::string>("jesUncFileName");      gSystem->ExpandPathName(uncFile);


  //
  // check input file                                            
  //
  TFile *evfile = TFile::Open(evurl);
  if(evfile==0) return -1;
  if(evfile->IsZombie()) return -1;
  TString proctag=gSystem->BaseName(evurl);
  Ssiz_t pos=proctag.Index(".root");
  proctag.Remove(pos,proctag.Length());
  TTree *tchain=getQCDTreeFrom(evfile);
  cout << "Processing " << tchain->GetEntriesFast() << " events for " << proctag << " @ " << evurl << endl;

    //
  // pileup reweighter
  //
  std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
  std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
  std::vector<float> mcPileupDistribution;
  if(isMC){
    //    TString puDist("evAnalyzer/h2zz/pileuptrue");
    //     TString puDist("evAnalyzer/h2zz/pileup");
    //     TH1F* histo = (TH1F *) evfile->Get(puDist);
    //     if(!histo)std::cout<<"pileup histogram is null!!!\n";
    //     for(int i=1;i<=histo->GetNbinsX();i++){mcPileupDistribution.push_back(histo->GetBinContent(i));}
    //     delete histo;
  }
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
  edm::LumiReWeighting *LumiWeights=0;
  if(isMC && mcPileupDistribution.size()) LumiWeights= new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);


  //book histograms
  SmartSelectionMonitor controlHistos;
  TH1 *Hcutflow=(TH1 *)controlHistos.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;
  for(int ibin=1; ibin<=6; ibin++) Hcutflow->SetBinContent(ibin,1);

  controlHistos.addHistogram( new TH1F("pileup", ";Pileup; Events",100,-0.5,99.5) );
  controlHistos.addHistogram( new TH1F("pileup_unwgt", ";Pileup; Events",100,-0.5,99.5) );
  controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 50, 0.,50.) );
  controlHistos.addHistogram( new TH1F ("nvertices_unwgt", "; Vertex multiplicity; Events", 50, 0.,50.) );

  TH1 * h = controlHistos.addHistogram( new TH1F ("njets", "; Jet multiplicity; Events", 6, 0.,6.) );
  h->GetXaxis()->SetBinLabel(1,"=0 jets");
  h->GetXaxis()->SetBinLabel(2,"=1 jets");
  h->GetXaxis()->SetBinLabel(3,"=2 jets");
  h->GetXaxis()->SetBinLabel(4,"= 3 jets");
  h->GetXaxis()->SetBinLabel(5,"= 4 jets");
  h->GetXaxis()->SetBinLabel(6,"#geq 5 jets");
  controlHistos.addHistogram( new TH1D("dphijj",";#Delta#phi(j^{(1)},j^{(2)});Events",100,-3.2,3.2) );

  //jet control
  for(size_t ijet=1; ijet<=4; ijet++)
    {
      TString jetctr(""); jetctr += ijet;
      controlHistos.addHistogram( new TH1F ("jet"+jetctr, "; Jet #"+jetctr+" p_{T} [GeV/c]; Events / (10 GeV/c)", 50, 0.,500.) );
      controlHistos.addHistogram( new TH1F ("jet"+jetctr+"eta", "; Jet #"+jetctr+" #eta; Events", 30, 0.,3.) );
    }

  TString jetFlavors[]={"","b","udsg","c"};
  const size_t nJetFlavors=sizeof(jetFlavors)/sizeof(TString);
  if(isMC)
    {
      h = controlHistos.addHistogram( new TH1F ("recoilflav", "; Flavor; Events", nJetFlavors, 0,nJetFlavors) );
      TH1 *hbv = controlHistos.addHistogram( new TH1F ("recoilflavbveto", "; Flavor; Events", nJetFlavors, 0,nJetFlavors) );
      for(size_t iflav=0; iflav<nJetFlavors; iflav++)
	{
	  h->GetXaxis()->SetBinLabel(iflav+1,jetFlavors[iflav]);
	  hbv->GetXaxis()->SetBinLabel(iflav+1,jetFlavors[iflav]);
	  controlHistos.addHistogram( new TH1F (jetFlavors[iflav]+"jetpt", ";"+jetFlavors[iflav]+" jet p_{T} [GeV/c]; Events / (10 GeV/c)", 50, 0.,500.) );
	  controlHistos.addHistogram( new TH1F (jetFlavors[iflav]+"jeteta", ";"+jetFlavors[iflav]+"jet #eta; Events", 30, 0.,3.) );
	  controlHistos.addHistogram( new TH2F (jetFlavors[iflav]+"jetptvseta", "; Jet p_{T} [GeV/c]; Jet #eta; Events / (10 GeV/c)", 50, 0.,500., 5, 0., 2.5) );
	} 
    }


  ///////////////////////////////////
  // b-tagging                     //
  ///////////////////////////////////
  TString btagger[]={"csv","jp","tchp","tche"};
  float btaggerMin[]={-0.2, 0.0,    -2, -2 };
  float btaggerMax[]={1.2,  2.5,  20, 20};
  int idxForBtaggerToTemplate=1;
  TString btaggerToTemplate=btagger[idxForBtaggerToTemplate];
  int idxForBtaggerForTemplatedEff=0;
  TString btaggerForTemplatedEff=btagger[idxForBtaggerForTemplatedEff]; 
  TString btaggerWPs[]={"L","M","T"};
  TString jetRanges[]={"inc","30to50","50to80","80to120","120toInf","0to0.5","0.5to1.0","1.0to1.5","1.5to2.5","0to10","11to14","15to18","18toInf"};
  const size_t nJetRanges=sizeof(jetRanges)/sizeof(TString);
  for(size_t i=0; i<sizeof(btagger)/sizeof(TString); i++)
    {
      controlHistos.addHistogram( new TH1F("inc"+btagger[i],";Discriminator;Jets", 50, btaggerMin[i],btaggerMax[i]) );
      
      //use templates based on a given discriminator except in the case we're measuring the discriminator efficiency itself (switch to alternative)
      int idxToUse=idxForBtaggerToTemplate;
      if(btagger[i]==btaggerToTemplate) idxToUse=idxForBtaggerForTemplatedEff;
      TH2 *hinc=(TH2 *)controlHistos.addHistogram( new TH2F(btagger[idxToUse],";Discriminator;Jets", 50, btaggerMin[idxToUse],btaggerMax[idxToUse],nJetRanges,0,nJetRanges) );
      TH2 *hb=(TH2 *)controlHistos.addHistogram( new TH2F(btagger[idxToUse]+"b",";b-tags;Range;Jets", 50,btaggerMin[idxToUse],btaggerMax[idxToUse],nJetRanges,0,nJetRanges) );
      TH2 *hc=(TH2 *)controlHistos.addHistogram( new TH2F(btagger[idxToUse]+"c",";c-tags;Range;Jets", 50,btaggerMin[idxToUse],btaggerMax[idxToUse],nJetRanges,0,nJetRanges) );
      TH2 *hudsg= (TH2 *)controlHistos.addHistogram( new TH2F(btagger[idxToUse]+"udsg",";udsg-tags;Range;Jets", 50,btaggerMin[idxToUse],btaggerMax[idxToUse],nJetRanges,0,nJetRanges) );
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
		  controlHistos.addHistogram( (TH2 *)hinc->Clone( btagger[idxToUse]       +btagger[k]+btaggerWPs[iwp]+pfix) );
		  controlHistos.addHistogram( (TH2 *)hb->Clone(   btagger[idxToUse]+"b"   +btagger[k]+btaggerWPs[iwp]+pfix) );
		  controlHistos.addHistogram( (TH2 *)hc->Clone(   btagger[idxToUse]+"c"   +btagger[k]+btaggerWPs[iwp]+pfix) );
		  controlHistos.addHistogram( (TH2 *)hudsg->Clone(btagger[idxToUse]+"udsg"+btagger[k]+btaggerWPs[iwp]+pfix) );
		}
	    }
	}
    }


  //get the jet pt weights
  std::map<TString, TGraph *> jetWeights;
  for(size_t iwf=0; iwf<jetPtWeightsFile.size(); iwf++)
    {
      TString inUrl("${CMSSW_BASE}/src/LIP/Top/"); inUrl += jetPtWeightsFile[iwf].c_str();
      gSystem->ExpandPathName(inUrl);
      TFile *inF=TFile::Open(inUrl);
      cout << "Retrieving jet weights from " << inUrl << endl;
      bool isEta(inUrl.Contains("_eta"));
      if(inF)
	{
	  for(size_t iflav=0; iflav<nJetFlavors; iflav++)
	    {
	      TH1 *hwgt=(TH1 *)inF->Get(jetFlavors[iflav]+(isEta?"jetetawgt":"jetptwgt"));
	      if(hwgt==0) continue;
	      jetWeights[jetFlavors[iflav]+(isEta?"eta":"pt")]=new TGraph(hwgt);
	    }
	  inF->Close();
	}
    }
  if(jetWeights.size())  cout << " jet templates will be weighted independently for " << jetWeights.size() << " flavors" << endl;

    
  //
  // LOOP OVER THE EVENTS
  //
  for (Int_t i = 0; i<tchain->GetEntriesFast(); i++) 
    {
      tchain->GetEntry(i);
      if(i%500==0) { printf("\r [ %d/100 ]",int(100*float(i)/float(tchain->GetEntriesFast()))); cout << flush; }

      //check the trigger
      int HLTCode = -1;
      bool PassTrigger = PassTriggerBit(qcdev.BitTrigger, HLTCode);
      if(!PassTrigger) continue;

      //pu weight
      float weight = LumiWeights ? LumiWeights->weight(qcdev.nPUtrue) : 1.0;
	           
      std::vector<TString> catsToFill;
      catsToFill.push_back("all");

      controlHistos.fillHisto("pileup_unwgt",catsToFill,qcdev.nPUtrue,1.0);
      controlHistos.fillHisto("pileup",catsToFill,qcdev.nPUtrue,weight);


      //jet kinematics, flavor templates, etc.
      int nbs(0),ncs(0),nudsg(0);
      std::map<TString,int> btagsCount; 
      btagsCount["tcheL"]=0;   btagsCount["tcheM"]=0;
      btagsCount["csvL"]=0;    btagsCount["csvM"]=0;   btagsCount["csvT"]=0;
      btagsCount["jpL"]=0;     btagsCount["jpM"]=0;    btagsCount["jpT"]=0;
      btagsCount["tchpT"]=0;
      PhysicsObjectJetCollection jets=getQCDJetCollection(jetPtCut,isMC);
      for(size_t ijet=0; ijet<jets.size(); ijet++)
	{
	  //jet flavor
	  TString jetFlav("udsg");
	  if(fabs(jets[ijet].flavid)==5)      { nbs++;   jetFlav="b"; }
	  else if(fabs(jets[ijet].flavid)==4) { ncs++;   jetFlav="c"; }
	  else                                { nudsg++; jetFlav="udsg"; }

	  //kinematics
	  float pt=jets[ijet].pt();
	  float eta=fabs(jets[ijet].eta());
	  float jetWgt=1.0;
	  for(int iwgt=0; iwgt<2; iwgt++)
	    {
	      TString pfix(iwgt==0?"pt":"eta");
	      if(jetWeights.find(jetFlav+pfix)!=jetWeights.end())
		{
		  TGraph *wgtGr=jetWeights[jetFlav+pfix];
		  if(wgtGr)
		    {
		      jetWgt=wgtGr->Eval(iwgt==0?pt:eta);
		    }
		}
	    }
	  
	  //reset the weight
	  float iweight=weight*jetWgt;

	  TString jetctr(""); jetctr+=ijet;
	  std::vector<int> jetCategs;
	  jetCategs.push_back(0);
	  if(pt>30 && pt<=50)      jetCategs.push_back(1);
	  if(pt>50 && pt<=80)      jetCategs.push_back(2);
	  if(pt>80 && pt<=120)     jetCategs.push_back(3);
	  if(pt>120)               jetCategs.push_back(4);
	  if(eta<=0.5)             jetCategs.push_back(5);
	  if(eta>0.5 && eta<=1.0)  jetCategs.push_back(6);
	  if(eta>1.0 && eta<=1.5)  jetCategs.push_back(7);
	  if(eta>1.5)              jetCategs.push_back(8);
	  if(qcdev.nPV<=10)              jetCategs.push_back(9);
	  if(qcdev.nPV>10 && qcdev.nPV<=14)    jetCategs.push_back(10);
	  if(qcdev.nPV>14 && qcdev.nPV<=18)    jetCategs.push_back(11);
	  if(qcdev.nPV>19)               jetCategs.push_back(12);
	  
	  controlHistos.fillHisto("jet"+jetctr,catsToFill,pt,iweight);
	  controlHistos.fillHisto("jet"+jetctr+"eta",catsToFill,fabs(eta),iweight);
	  controlHistos.fillHisto( "jetpt", catsToFill,pt,iweight);
	  controlHistos.fillHisto( "jeteta",catsToFill,fabs(eta),iweight);
	  if(isMC)
	    {
	      controlHistos.fillHisto( jetFlav+"jetpt", catsToFill,pt,iweight);
	      controlHistos.fillHisto( jetFlav+"jeteta",catsToFill,fabs(eta),iweight);
	      controlHistos.fillHisto( jetFlav+"jetptvseta", catsToFill,pt,fabs(eta),iweight);
	    }
   
	  //b-tagging
	  std::map<TString, bool> hasTagger;
	  float tche = jets[ijet].btag1;
	  bool hasTCHEL(tche>1.7);  btagsCount["tcheL"] += hasTCHEL; hasTagger["tcheL"]=hasTCHEL;
	  bool hasTCHEM(tche>1.93); btagsCount["tcheM"] += hasTCHEM; hasTagger["tcheM"]=hasTCHEM; 
	  float csv  = jets[ijet].btag2;
	  bool hasCSVL(csv>0.244);  btagsCount["csvL"] += hasCSVL;   hasTagger["csvL"]=hasCSVL;
	  bool hasCSVM(csv>0.679);  btagsCount["csvM"] += hasCSVM;   hasTagger["csvM"]=hasCSVM;
	  bool hasCSVT(csv>0.898);  btagsCount["csvT"] += hasCSVT;   hasTagger["csvT"]=hasCSVT;
	  float jp   = jets[ijet].btag3;
	  bool hasJPL(jp>0.275);    btagsCount["jpL"] += hasJPL;     hasTagger["jpL"]=hasJPL; 
	  bool hasJPM(jp>0.545);    btagsCount["jpM"] += hasJPM;     hasTagger["jpM"]=hasJPM;
	  bool hasJPT(jp>0.790);    btagsCount["jpT"] += hasJPT;     hasTagger["jpT"]=hasJPT;
	  float tchp = jets[ijet].btag4;
	  bool hasTCHPT(tchp>3.41); btagsCount["tchpT"] += hasTCHPT; hasTagger["tchpT"]=hasTCHPT;
	  controlHistos.fillHisto("inctche", catsToFill, tche,  iweight);
	  controlHistos.fillHisto("inctchp", catsToFill, tchp,  iweight);
	  controlHistos.fillHisto("inccsv",  catsToFill, csv,   iweight);
	  controlHistos.fillHisto("incjp",   catsToFill, jp,    iweight);
	  float btaggerToTemplateVal=csv;
	  if(btaggerToTemplate=="jp")   btaggerToTemplateVal=jp;
	  if(btaggerToTemplate=="tche") btaggerToTemplateVal=tche;
	  if(btaggerToTemplate=="tchp") btaggerToTemplateVal=tchp;
	  float btaggerForTemplatedEffVal=jp;
	  if(btaggerForTemplatedEff=="csv")  btaggerForTemplatedEffVal=csv;
	  if(btaggerForTemplatedEff=="tche") btaggerForTemplatedEffVal=tche;
	  if(btaggerForTemplatedEff=="tchp") btaggerForTemplatedEffVal=tchp;
	  for(size_t ijcat=0; ijcat<jetCategs.size(); ijcat++)
	    {
	      controlHistos.fillHisto(btaggerToTemplate, catsToFill, btaggerToTemplateVal, jetCategs[ijcat], iweight);
	      controlHistos.fillHisto(btaggerForTemplatedEff, catsToFill, btaggerForTemplatedEffVal, jetCategs[ijcat], iweight);
	      if(isMC)
		{
		  controlHistos.fillHisto(btaggerToTemplate+jetFlav, catsToFill, btaggerToTemplateVal, jetCategs[ijcat], iweight);
		  controlHistos.fillHisto(btaggerForTemplatedEff+jetFlav, catsToFill, btaggerForTemplatedEffVal, jetCategs[ijcat], iweight);
		}
  
	      for(std::map<TString, bool>::iterator it=hasTagger.begin(); it!=hasTagger.end(); it++)
		{
		   TString pfix(it->second?"":"v");
		  float idisc=btaggerToTemplateVal;
		  TString idiscName=btaggerToTemplate;
		  if(it->first.Contains(btaggerToTemplate)) { idisc=btaggerForTemplatedEffVal; idiscName=btaggerForTemplatedEff; }
		  controlHistos.fillHisto( idiscName+it->first+pfix,catsToFill, idisc, jetCategs[ijcat], iweight);
		  if(isMC) controlHistos.fillHisto( idiscName+jetFlav+it->first+pfix,catsToFill, idisc, jetCategs[ijcat], iweight);
		}
	    }
	}

      controlHistos.fillHisto("njets",catsToFill,jets.size(),weight);
      if(jets.size()<2) continue;
      float dphijj=deltaPhi(jets[0].phi(),jets[1].phi());
      controlHistos.fillHisto("dphijj",catsToFill,dphijj,weight);
      if(fabs(dphijj)<2.6) continue;
      int rfBin(2);
      if(fabs(jets[1].flavid)==5) rfBin=1;
      if(fabs(jets[1].flavid)==4) rfBin=3;


      controlHistos.fillHisto("recoilflav",catsToFill,0,weight);
      controlHistos.fillHisto("recoilflav",catsToFill,rfBin,weight);
      if(jets[0].btag3>0.275) continue;
      controlHistos.fillHisto("recoilflavbveto",catsToFill,0,weight);
      controlHistos.fillHisto("recoilflavbveto",catsToFill,rfBin,weight);
	
    }
  evfile->Close();

  //
  // save histos to local file
  //
  TString outUrl(outdir);
  gSystem->ExpandPathName(outUrl);
  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  outUrl += proctag + ".root";
  TFile *file=TFile::Open(outUrl, "recreate");
  controlHistos.Write();
  file->Close();
 }
