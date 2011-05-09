#include <iostream>
#include <boost/shared_ptr.hpp>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/KinResultsHandler.h"
#include "LIP/Top/interface/KinAnalysis.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

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
  if ( argc < 5 ) {
    std::cout << "Usage : " << argv[0] << " kinDir eventsFile outputDir isMC" << std::endl;
    return 0;
  }

  std::map<TString, TH1 *> results;
  TString cats[]={"all","ee","mumu","emu"};
  size_t ncats=sizeof(cats)/sizeof(TString);
  for( size_t icat=0; icat<ncats; icat++)
    {
      results[cats[icat]+"_njets"] = (TH1*) formatPlot( new TH1F (cats[icat]+"_njets", ";Jets;Events", 6, 0.,6.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_btags"] = (TH1*) formatPlot( new TH1F (cats[icat]+"_btags", ";b-tag multiplicity;Events", 6, 0.,6.), 1,1,1,20,0,true,true,1,1,1);
      for(int ibin=1; ibin<=results[cats[icat]+"_njets"]->GetXaxis()->GetNbins(); ibin++)
	{
	  TString label(""); label += ibin-1; 
	  if(ibin==results[cats[icat]+"_njets"]->GetXaxis()->GetNbins()) label ="#geq" + label;
	  results[cats[icat]+"_njets"]->GetXaxis()->SetBinLabel(ibin,label + " jets");
	  results[cats[icat]+"_btags"]->GetXaxis()->SetBinLabel(ibin,label + " btags");
	}

      results[cats[icat]+"_leadjet"] = (TH1*) formatPlot( new TH1F (cats[icat]+"_leadjet", "; Leading jet p_{T} [GeV/c]; Events / (5 GeV/c)", 50, 0.,250.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_subleadjet"] = (TH1*) formatPlot( new TH1F (cats[icat]+"_subleadjet", "; Sub-leading jet p_{T} [GeV/c]; Events / (5 GeV/c)", 50, 0.,250.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_leadlepton"] = (TH1*) formatPlot( new TH1F (cats[icat]+"_leadlepton", "; Leading lepton p_{T} [GeV/c]; Events / (5 GeV/c)", 50, 0.,250.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_subleadlepton"] = (TH1*) formatPlot( new TH1F (cats[icat]+"_subleadlepton", "; Sub-leading lepton p_{T} [GeV/c]; Events / (5 GeV/c)", 50, 0.,250.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_met"] = (TH1*) formatPlot( new TH1F (cats[icat]+"_met", "; #slash{E}_{T} [GeV/c]; Events / (5 GeV/c)", 50, 0.,250.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_mtop"] = (TH1*) formatPlot( new TH1F (cats[icat]+"_mtop", "; m_{Top} [GeV/c^{2}]; Events / (5 GeV/c^{2})", 100, 0.,500.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_afb"] = (TH1*) formatPlot( new TH1F (cats[icat]+"_afb", "; #Delta #eta(t,#bar{t}); Events / (0.1)", 100, -5.,5.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_mttbar"] = (TH1*) formatPlot( new TH1F (cats[icat]+"_mttbar", "; Mass(t,#bar{t}) [GeV/c^{2}]; Events / (20 GeV/c^{2})", 100, 0.,2000.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_mtopvsdilmass"] = (TH1F*)formatPlot( new TH2F (cats[icat]+"_mtopvsdilmass", "; m_{Top} [GeV/c^{2}]; Mass(l,l') [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,500.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_mtopvsmlj"] = (TH1F*)formatPlot( new TH2F (cats[icat]+"_mtopvsmlj", "; m_{Top} [GeV/c^{2}]; Mass(l,j) [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,500.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_mtopvsmet"] = (TH1F*)formatPlot( new TH2F (cats[icat]+"_mtopvsmet", "; m_{Top} [GeV/c^{2}]; #slash{E}_{T} [GeV/c]; Events", 100, 0.,500.,50, 0.,250.), 1,1,1,20,0,true,true,1,1,1);

      results[cats[icat]+"_mtopvsmttbar"] = (TH1F*)formatPlot( new TH2F (cats[icat]+"_mtopvsmttbar", "; m_{Top} [GeV/c^{2}]; Mass(t,#bar{t}) [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,2000.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_mtopvsafb"] = (TH1F*)formatPlot( new TH2F (cats[icat]+"_mtopvsafb", "; m_{Top} [GeV/c^{2}]; #Delta #eta(t,#bar{t}); Events", 100, 0.,500.,100, -5.,5.), 1,1,1,20,0,true,true,1,1,1);
      results[cats[icat]+"_mttbarvsafb"] = (TH1F*)formatPlot( new TH2F (cats[icat]+"_mttbarvsafb", "; Mass(t,#bar{t}) [GeV/c^{2}];#Delta #eta(t,#bar{t}); Events", 100, 0.,2000.,100,-5.,5.), 1,1,1,20,0,true,true,1,1,1);
    }



  //fix entries flag
  TString isMCBuf(argv[4]);
  bool isMC=isMCBuf.Atoi();

  //process kin file

  TString url = argv[1];
  gSystem->ExpandPathName(url);
  KinResultsHandler kinHandler;
  kinHandler.init(url,false);
  TChain *t=kinHandler.getResultsChain();

  //process events file
  TString evurl = argv[2];
  gSystem->ExpandPathName(evurl);
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
  std::map<TString,int> selEvents;
  for (int inum=0; inum < evTree->GetEntriesFast(); ++inum)
  {
    evTree->GetEvent(inum);
    EventSummary_t &ev = evSummaryHandler.getEvent();
    TString key(""); key+= ev.run; key+="-"; key += ev.lumi; key+="-"; key += ev.event;
    selEvents[key]=inum;
  }

  //loop over kin results
  int nresults(0);
  for (int inum=0; inum < t->GetEntries(); ++inum){
    t->GetEvent(inum);
 
    //get original event
    Int_t irun,ievent,ilumi;
    kinHandler.getEventInfo(irun,ievent,ilumi);


    TString key("");  key+= irun; key+="-"; key += ilumi;  key+="-"; key += ievent;

    if(selEvents.find(key)==selEvents.end()) continue;
    nresults++;
    evTree->GetEntry( selEvents[key] );

    //get event summary
    EventSummary_t &ev = evSummaryHandler.getEvent();
    std::vector<TString> categs;
    categs.push_back("all");
    if(ev.cat==dilepton::MUMU)  categs.push_back("mumu");
    if(ev.cat==dilepton::EE)  categs.push_back("ee");
    if(ev.cat==dilepton::EMU)  categs.push_back("emu");


    //get particles from the event
    

    int njets(0),nbtags(0);
    KinCandidateCollection_t leptons, jets, mets;

    for(Int_t ipart=0; ipart<ev.nparticles; ipart++)
      {
	TLorentzVector p4(ev.px[ipart],ev.py[ipart],ev.pz[ipart],ev.en[ipart]);
	if(isnan(p4.Pt()) || isinf(p4.Pt())) continue;
	switch( ev.id[ipart] )
	  {
	  case 0:
	    mets.push_back( KinCandidate_t(p4,p4.Pt()) );
            break;
          case 1:
            jets.push_back( KinCandidate_t(p4, ev.info1[ipart]) );
	    njets++;
	    if(ev.info1[ipart]>1.74) nbtags++;
            break;
          default:
            leptons.push_back( KinCandidate_t(p4,ev.id[ipart]) );
            break;

	  }
      }

    sort(leptons.begin(),leptons.end(),KinAnalysis::sortKinCandidates);
    sort(jets.begin(),jets.end(),KinAnalysis::sortKinCandidates);
    sort(mets.begin(),mets.end(),KinAnalysis::sortKinCandidates);



    //get the combination preferred by KIN
    TH1F *h1=kinHandler.getHisto("mt",1), *h2=kinHandler.getHisto("mt",2);
    Int_t icomb=(h1->Integral()< h2->Integral())+1;
    TH1F *mpref=kinHandler.getHisto("mt",icomb);
    double mtop = kinHandler.getMPVEstimate(mpref) [1];
    TH1F *mttbarpref=kinHandler.getHisto("mttbar",icomb);
    double mttbar = kinHandler.getMPVEstimate(mttbarpref)[1];
    TH1F *afbpref=kinHandler.getHisto("afb",icomb);
    double afb = kinHandler.getMPVEstimate(afbpref)[1];
    
    //compute dilepton invariant mass
    TLorentzVector dil = leptons[0].first+leptons[1].first;

    float dilmass = dil.M();

    //get the lepton-jet pairs
    TLorentzVector lj1=leptons[0].first+jets[icomb==1?0:1].first;
    TLorentzVector lj2=leptons[1].first+jets[icomb==1?1:0].first;


    //fill histos
    float weight = ev.weight;
    for(std::vector<TString>::iterator cIt = categs.begin(); cIt != categs.end(); cIt++)
      {
	double ptjet1(jets[0].first.Pt()), ptjet2(jets[1].first.Pt());
	double ptlep1(leptons[0].first.Pt()), ptlep2(leptons[1].first.Pt());
	results[*cIt+"_njets"]->Fill(njets,weight);
	results[*cIt+"_btags"]->Fill(nbtags,weight);
	results[*cIt+"_leadjet"]->Fill(max(ptjet1,ptjet2),weight);
	results[*cIt+"_subleadjet"]->Fill(min(ptjet1,ptjet2),weight);
	results[*cIt+"_leadlepton"]->Fill(max(ptlep1,ptlep2),weight);
	results[*cIt+"_subleadlepton"]->Fill(min(ptlep1,ptlep2),weight);
	results[*cIt+"_met"]->Fill(mets[0].first.Pt(),weight);
	if(mtop>0)
	  {
	    results[*cIt+"_mtop"]->Fill(mtop,weight);
	    ((TH2F *)results[*cIt+"_mtopvsdilmass"])->Fill(mtop,dilmass,weight);
	    ((TH2F *)results[*cIt+"_mtopvsmlj"])->Fill(mtop,lj1.M(),weight);
	    ((TH2F *)results[*cIt+"_mtopvsmlj"])->Fill(mtop,lj2.M(),weight);
	    ((TH2F *)results[*cIt+"_mtopvsmet"])->Fill(mtop,mets[0].first.Pt(),weight);
	    ((TH2F *)results[*cIt+"_mtopvsmttbar"])->Fill(mtop,mttbar,weight);
	    ((TH2F *)results[*cIt+"_mtopvsafb"])->Fill(mtop,afb,weight);
	    ((TH2F *)results[*cIt+"_mttbarvsafb"])->Fill(mttbar,afb,weight);
	    results[*cIt+"_afb"]->Fill(afb);
	    results[*cIt+"_mttbar"]->Fill(mttbar);
	  }
      }



    if (!isMC && mtop>350) cout << irun << ":" << ilumi << ":" << ievent << " " << "(" << leptons[0].second <<" "<< leptons[1].second << ")" << flush;
   
  }
  kinHandler.end();
  cout << endl;


  //if MC: rescale to number of selected events and to units of pb
  if(isMC && nresults)

    {
      double scaleFactor=selEvents.size()/nresults;
      TString tag=gSystem->BaseName(evurl);
      tag.ReplaceAll(".root","");
      TH1F *cutflowH = (TH1F *) evfile->Get("evAnalyzer/"+tag+"/cutflow");
      if(cutflowH)
	{
	  float cnorm=cutflowH->GetBinContent(1);
	  if(cnorm>0) scaleFactor/=cnorm;
	}
      cout << scaleFactor << endl;
      for(std::map<TString,TH1 *>::iterator hIt = results.begin(); hIt != results.end(); hIt++) hIt->second->Scale(scaleFactor);
    }


  //save to file
  TString outUrl( argv[3] );
  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  outUrl += gSystem->BaseName(evurl);
  TFile *file=TFile::Open(outUrl, "recreate");

  TDirectory *baseOutDir=file->mkdir("massAnalyzer");
  for( size_t icat=0; icat<ncats; icat++)
    {
      baseOutDir->cd();
      if(icat) baseOutDir->mkdir( cats[icat] )->cd();
      for(std::map<TString,TH1 *>::iterator hIt = results.begin(); hIt != results.end(); hIt++) 
	{
	  if(!hIt->first.BeginsWith(cats[icat])) continue;
	  hIt->second->Write();
	}
    }

  file->Close(); 
}  
