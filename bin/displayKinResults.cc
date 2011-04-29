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
#include "TTree.h"
#include "TCanvas.h"

using namespace std;

//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

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
  KinResultsHandler kinHandler(url,false);
  TTree *t=kinHandler.getResultsTree();
  const Int_t nEntries = t->GetEntriesFast();
  if(event>=0 && event<nEntries)
    {
      t->GetEntry(event);
      cout << "event retrieved" << endl;

      Int_t irun,ievent,ilumi;
      kinHandler.getEventInfo(irun,ievent,ilumi);
      TString title("Run: "); title += irun; title += "\\"; 
      title += "Event: ";     title += ievent; title += "\\"; 
      title += "Lumi: ";      title += ilumi; title += "\\"; 
      cout << title << endl;

      setStyle();
      TCanvas *c=getNewCanvas("kinres","kinres",false);
      TH1D *h =kinHandler.getHisto("mt",1);
      cout << h << " " << dynamic_cast<TH1D *>(h) << endl;
      h->Print("all");
      h->Draw("hist");
      h=kinHandler.getHisto("mt",2);
      cout << h << endl;
      h->Draw("histsame");

      TLegend *leg = c->BuildLegend();
      formatForCmsPublic(c,leg,title,2);
      c->SaveAs("kinres.C");
    }

  kinHandler.end();
}  
