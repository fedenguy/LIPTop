#include <iostream>
#include <sstream>
#include <boost/shared_ptr.hpp>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/KinResultsHandler.h"
#include "LIP/Top/interface/KinAnalysis.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "LIP/Top/interface/HistogramAnalyzer.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;
using namespace top;

typedef std::pair<Double_t,Double_t> minMax;
typedef std::vector<minMax> minMaxCollection;
typedef std::vector<Int_t> histogramEntries;
//


Double_t getBaselineFromHisto(TH1F& h){
  Double_t min;
  bool minLock(false);
  
  Double_t max;

  for(Int_t i=0; i< h.GetNbinsX(); i++){
    if(h.GetBinContent(i) != 0){

      if(minLock == false){
	min = h.GetBinCenter(i);
	minLock = true;
      }
      
      max = h.GetBinCenter(i);
    }
  }
  
  cout << "Max: " << max << ", Min: " << min << ", N: " << h.GetEntries() << endl;
  Double_t result = (max - min) / fabs(10000 - h.GetEntries()) ;
  return result;
}


TH1F processSingleEvent(TString url, TString ourl, Int_t event, TString& eventTag, histogramEntries& histEntries){
  //open the file
  KinResultsHandler kinHandler;
  kinHandler.init(url,false);
  TChain* t=kinHandler.getResultsChain();

  HistogramAnalyzer histoAnalyzer;
  const Int_t nEntries = t->GetEntries();
  minMaxCollection mMc;
    
  if(event>=0 && event<nEntries)
    {
      cout << "enter sandman" << endl;
      t->GetEntry(event);
      cout << "loaded event" << endl;
      Int_t irun,ievent,ilumi;
      kinHandler.getEventInfo(irun,ievent,ilumi);
      TString title("CMS preliminary, #sqrt{s}=7 TeV\\Run: "); title += irun; title += "\\"; 
      title += "Event: ";     title += ievent; title += "\\"; 
      title += "Lumi: ";      title += ilumi; title += "\\"; 

      stringstream evtTags;
      evtTags << ievent;
      TString evtTag= TString(evtTags.str());
      eventTag = evtTag;

      setStyle();
      TCanvas* c=getNewCanvas("kinres"+evtTag,"kinres"+evtTag,false);

      TH1F* h=kinHandler.getHisto("mt",1);
      TH1F* h1= kinHandler.getHisto("mt",1);

      histEntries.push_back(h->GetEntries());
      h->SetLineWidth(2);
      h->SetMarkerStyle(20);
      h->SetFillStyle(0);
      h->DrawClone("hist");
      std::map<TString,Double_t> res=histoAnalyzer.analyzeHistogram(h);
      cout << "Combination #1: " << h->Integral() << " has solutions" << endl;
      for(std::map<TString, Double_t>::iterator it = res.begin(); it != res.end(); it++) cout << "\t" << it->first << "=" << it->second << endl;

      h =kinHandler.getHisto("mt",2);
      TH1F* h2 = kinHandler.getHisto("mt",2);
      histEntries.push_back(h->GetEntries());
      h->SetFillStyle(3472);
      h->SetFillColor(2);
      h->SetMarkerStyle(21);
      h->DrawClone("histsame");
      res=histoAnalyzer.analyzeHistogram(h);
      cout << "Combination #2: " << h->Integral() << " has solutions" << endl;
      for(std::map<TString, Double_t>::iterator it = res.begin(); it != res.end(); it++) cout << "\t" << it->first << "=" << it->second << endl;

      TLegend *leg=c->BuildLegend();
      formatForCmsPublic(c,leg,title,2);
      c->Modified();
      c->Update();
      c->SaveAs(ourl+"kinres"+evtTag+".C");
      c->SaveAs(ourl+"kinres"+evtTag+".png");

      // add 
      Double_t baseline_1 = getBaselineFromHisto(*h1);
      Double_t baseline_2 = getBaselineFromHisto(*h2);
      
      cout << "baseline for combination 1: " << baseline_1 << endl;
      cout << "baseline for combination 2: " << baseline_2 << endl;

      for(Int_t i=0; i<h1->GetEntries(); i++)
	h1->AddBinContent(i, baseline_1);
      for(Int_t i=0; i<h2->GetEntries(); i++)
	h2->AddBinContent(i, baseline_2);
      
      h1->Add(h2);

      cout << "finishing single event" << endl;
      kinHandler.end();
      TH1F my = *h1;
      cout << "entries in kinideogram for event " << ievent << " are " << my.GetEntries() << endl;
      return *h1;    
    }
  TH1F blah;
  kinHandler.end();
  return blah;
  
}

void processNEvents(TString url, TString ourl, Int_t event){

  TH1F* kinideogram = new TH1F("kinideogram","M_{t} likelihood from KINb solutions;M_{t} [GeV/c^{2}];Likelihood / (5 GeV/c^{2})",400,0,2000);

  kinideogram->Reset();

  //open the file
  KinResultsHandler kinHandler;
  kinHandler.init(url,false);
  TChain* t=kinHandler.getResultsChain();

  const Int_t nEntries = t->GetEntries();
  kinHandler.end();

  if(event > nEntries)
    event = nEntries;

  for(Int_t nEvent=0; nEvent < event; nEvent++){

    histogramEntries histEntries;
    histEntries.clear(); // or erase?
    TString eventTag;
    TH1F eventideogram = processSingleEvent(url, ourl, nEvent, eventTag, histEntries);
    cout << "entries in kinideogram for event " << eventTag << " are " << eventideogram.GetEntries() << endl;    
    
    TCanvas* c=getNewCanvas("kinideogram"+eventTag,"kinideogram"+eventTag,false);
    eventideogram.DrawClone("hist");
    
    TLegend *leg=c->BuildLegend();
    formatForCmsPublic(c,leg,"Event "+eventTag+" M_{t} likelihood",2);
    c->Modified();
    c->Update();
    c->SaveAs(ourl+"kinideogram"+eventTag+".C");
    c->SaveAs(ourl+"kinideogram"+eventTag+".png");
    
    for(Int_t i=0; i<kinideogram->GetEntries(); i++)
      kinideogram->AddBinContent(i,TMath::Log(eventideogram.GetBinContent(i)));
  }
  

  setStyle();
  TCanvas* c=getNewCanvas("kinres","kinres",false);
  kinideogram->DrawClone("hist");

  TLegend *leg=c->BuildLegend();
  formatForCmsPublic(c,leg,"M_{t} likelihood",2);
  c->Modified();
  c->Update();
  c->SaveAs(ourl+"kinideogram.C");
  c->SaveAs(ourl+"kinideogram.png");

  kinideogram->Delete();
  return;  
}

int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  //  HistogramAnalyzer histoAnalyzer;

  //check arguments
  if ( argc < 4 ) {
    std::cout << "Usage : " << argv[0] << " [kinDir] [outDir] [single event || -(number of events)]" << std::endl;
    return 0;
  }

  TString url = argv[1];
  TString ourl= argv[2];
  gSystem->Exec("mkdir -p " + ourl);
  ourl = ourl+"/";
  Int_t event(-1);
  sscanf(argv[3],"%d",&event);
  event--;


  if(event>=0){
    histogramEntries histEntries;
    histEntries.clear(); // or erase?
    cout << "single event start" << endl;
    TString eventTag;
    TH1F blah = processSingleEvent(url, ourl, event, eventTag, histEntries);
    cout << "single event stop" << endl;
  }
  else{
    event = -1*event;
    cout << "multiple events start" << endl;
    processNEvents(url, ourl, event);    
  }
  

}  
