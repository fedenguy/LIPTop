#include <iostream>
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

//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  HistogramAnalyzer histoAnalyzer;

  //check arguments
  if ( argc < 3 ) {
    std::cout << "Usage : " << argv[0] << " [kinDir] [event]" << std::endl;
    return 0;
  }

  TString url = argv[1];
  int event(-1);
  sscanf(argv[2],"%d",&event);
  event--;

  //open the file
  KinResultsHandler kinHandler;
  kinHandler.init(url,false);
  TChain *t=kinHandler.getResultsChain();
  const Int_t nEntries = t->GetEntries();
  if(event>=0 && event<nEntries)
    {
      t->GetEntry(event);

      Int_t irun,ievent,ilumi;
      kinHandler.getEventInfo(irun,ievent,ilumi);
      TString title("CMS preliminary, #sqrt{s}=7 TeV\\Run: "); title += irun; title += "\\"; 
      title += "Event: ";     title += ievent; title += "\\"; 
      title += "Lumi: ";      title += ilumi; title += "\\"; 

      setStyle();
      TCanvas *c=getNewCanvas("kinres","kinres",false);

      TH1F *h=kinHandler.getHisto("mt",1);
      h->SetLineWidth(2);
      h->SetMarkerStyle(20);
      h->SetFillStyle(0);
      h->DrawClone("hist");
      std::map<TString,Double_t> res=histoAnalyzer.analyzeHistogram(h);
      cout << "Combination #1: " << h->Integral() << " has solutions" << endl;
      for(std::map<TString, Double_t>::iterator it = res.begin(); it != res.end(); it++) cout << "\t" << it->first << "=" << it->second << endl;

      h =kinHandler.getHisto("mt",2);
      h->SetFillStyle(3472);
      h->SetFillColor(1);
      h->SetMarkerStyle(21);
      h->DrawClone("histsame");
      res=histoAnalyzer.analyzeHistogram(h);
      cout << "Combination #2: " << h->Integral() << " has solutions" << endl;
      for(std::map<TString, Double_t>::iterator it = res.begin(); it != res.end(); it++) cout << "\t" << it->first << "=" << it->second << endl;

      TLegend *leg=c->BuildLegend();
      formatForCmsPublic(c,leg,title,2);
      c->Modified();
      c->Update();
      c->SaveAs("kinres.C");
    }

  kinHandler.end();
}  
