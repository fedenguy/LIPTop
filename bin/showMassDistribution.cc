#include <iostream>
#include <boost/shared_ptr.hpp>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/KinResultsHandler.h"
#include "LIP/Top/interface/KinAnalysis.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCanvas.h"

using namespace std;

//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  //check arguments
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [kinDir] " << std::endl;
    return 0;
  }
  TH1F *graph=(TH1F*)formatPlot( new TH1F ("top_mass", "; m_{Top} [GeV/c^{2}]; Events", 100, 0.,500.), 1,1,1,20,0,true,true,1,1,1);
  TString url = argv[1];
  
  //process the file
  KinResultsHandler kinHandler;
  kinHandler.init(url,false);
  TChain *t=kinHandler.getResultsChain();
  const Int_t nEntries = t->GetEntries();
  cout << nEntries << endl;
  for (int inum=0; inum < nEntries; ++inum){
    t->GetEvent(inum);

    TH1F *h1=kinHandler.getHisto("mt",1), *h2=kinHandler.getHisto("mt",2);
    TH1F *hprev=h1->Integral()> h2->Integral()? h1:h2;
    vector<double> res=kinHandler.getMPVEstimate(hprev);
    double mtop=res[1];
      
    Int_t irun,ievent,ilumi;
    kinHandler.getEventInfo(irun,ievent,ilumi);
      
    if(inum%100==0)
      {
	cout << inum << " "
	     << irun << " " << ievent << " " << ilumi << endl;
	//	       << h1->Integral() << " " << h2->Integral() << endl
	//	       << mtop << " " << endl;
      }

    graph->Fill(mtop);
  }

  kinHandler.end();

  
  TString title("CMS preliminary, #sqrt{s}=7 TeV "); 
  setStyle();

  TCanvas *c=getNewCanvas("kinresprev","kinresprev",false);
  graph->SetLineWidth(2);
  graph->SetMarkerStyle(20);
  graph->SetFillStyle(0);
  graph->Draw("E1P");
    
  TLegend *leg=c->BuildLegend();
  formatForCmsPublic(c,leg,title,2);
  c->Modified();
  c->Update();
  c->SaveAs("kinresprev.C");
}  
