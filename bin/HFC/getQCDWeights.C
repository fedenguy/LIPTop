#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>
#include <vector>

using namespace std;

//
void getQCDWeights(TString llUrl,TString qcdUrl)
{
  TString flavors[]={"b","udsg","c"};
  const size_t nflavs=sizeof(flavors)/sizeof(TString);
  TString histos[]={"jetpt","jeteta","jetptvseta"};
  const size_t nhistos=sizeof(histos)/sizeof(TString);
  std::vector<TH1 *> wgtColl;

  //get histos from files
  TFile *inLLF=TFile::Open(llUrl);
  TFile *inQCDF=TFile::Open(qcdUrl);
  if(inLLF==0 || inQCDF==0) 
    {
      cout << "Error opening files: " << llUrl << " / " << qcdUrl << endl;
      return;
    }
  for(size_t iflav=0; iflav<nflavs; iflav++)
    {
      for(size_t ih=0; ih<nhistos; ih++)
	{
	  TH1 *qcdH = (TH1 *)inQCDF->Get("Multijets/"+flavors[iflav]+histos[ih]); qcdH->SetDirectory(0);
	  TH1 *ttH  = (TH1 *)inLLF->Get("t#bar{t}/"+flavors[iflav]+histos[ih]);   ttH->SetDirectory(0);
	  if(qcdH==0 || ttH==0) { cout << flavors[iflav]+histos[ih] << " not found" << endl; continue; }
	  
	  //divide the normalized flavor templates to get the weights
	  TH1 *qcdWgtH= (TH1 *)ttH->Clone(flavors[iflav]+histos[ih]+"wgt");
	  qcdWgtH->SetDirectory(0);
	  qcdWgtH->Reset("ICE");
	  qcdWgtH->Divide(ttH,qcdH,1./ttH->Integral(),1./qcdH->Integral());
	  wgtColl.push_back(qcdWgtH);

	  if(ih==0)
	    {
	      TCanvas *c=new TCanvas(flavors[iflav],flavors[iflav]);
	      c->cd();
	      ttH->SetTitle("t#bar{t}"); ttH->SetFillColor(kGray); ttH->SetFillStyle(1001); ttH->SetName(flavors[iflav]+"ttbar"); ttH->SetMarkerStyle(0);
	      qcdH->SetLineWidth(2);  qcdH->SetName(flavors[iflav]+"mj");  qcdH->SetTitle("Multijets"); qcdH->SetMarkerStyle(0);
	      ttH->GetYaxis()->SetTitle(flavors[iflav]+"-jets");
	      ttH->Draw("hist");   
	      qcdH->Scale(ttH->Integral()/qcdH->Integral());
	      qcdH->Draw("hist same");
	      ttH->GetYaxis()->SetRangeUser(0,qcdH->GetMaximum()*1.1);
	      TLegend *leg=c->BuildLegend();
	      leg->SetBorderSize(0);
	      leg->SetTextFont(42);
	      leg->SetFillStyle(0);
	      leg->SetHeader("CMS simulation");
	    }
	}
    }
  inLLF->Close();
  inQCDF->Close();

  //save weights to file
  TFile *fOut=TFile::Open("QCDweights.root","RECREATE");
  for(size_t ih=0; ih<wgtColl.size(); ih++) wgtColl[ih]->Write();
  fOut->Close();
  cout << "Weights available @ QCDweights.root" << endl;
}
