/*
  gSystem->Load("libCMGToolsHtoZZ2l2nu.so");
  .L getCalibrationResults.C+
*/

#include "TString.h"
#include "TH1F.h"
#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TLine.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TObjArray.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "TArrow.h"

#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/plotter.h"

void getCalibrationResults(TString url="TopMassCalibration.root",int nCategs=4)
{
  TString mpts[]={"161.5","163.5","166.5","169.5","172.5","175.5","178.5","181.5","184.5"};
  const size_t nmpts=sizeof(mpts)/sizeof(TString);

  std::vector<TString> labels;
  if(nCategs==2)
    {
      labels.push_back("1 b-tags");
      labels.push_back("#geq 2 b-tags");
    }
  else
    {
      labels.push_back("1 b-tags (ee+#mu#mu)");
      labels.push_back("#geq 2 b-tags (ee+#mu#mu)");
      labels.push_back("1 b-tags (e#mu)");
      labels.push_back("#geq 2 b-tags (e#mu)");
    }
  
  TString outDir=gSystem->DirName(url);
  cout << "Plots will be stored at @ " << outDir << endl;

  TGraphErrors *lingr = new TGraphErrors;  
  lingr->SetMarkerStyle(20);
  TGraphErrors *biasgr = new TGraphErrors; 
  biasgr->SetMarkerStyle(20);
  std::vector<TGraphErrors *> categbiasgr;
  for(int icat=0; icat<nCategs; icat++) categbiasgr.push_back( new TGraphErrors );
  TObjArray statUncResults;
  TObjArray pullResults;

  //get results from file
  TFile *fin=TFile::Open(url);
  for(size_t i=0; i<nmpts; i++)
    {
      Float_t trueMass=mpts[i].Atof();
      
      TH1D *massFitH=(TH1D *)fin->Get(mpts[i]+"/massfit");
      TF1 *ffunc=massFitH->GetFunction("gaus");
      float avgMass=ffunc->GetParameter(1);
      float avgMassErr=ffunc->GetParError(1);
      int ipt=lingr->GetN();
      lingr->SetPoint(ipt,trueMass,avgMass);
      lingr->SetPointError(ipt,0,avgMassErr);

      for(int icat=0; icat<=nCategs;icat++)
	{
	  TString postfix("");
	  TGraphErrors *gr=biasgr;
	  if(icat<nCategs){ postfix ="_"; postfix +=( icat+1); gr=categbiasgr[icat]; }
	  TH1D *biasH=(TH1D *)fin->Get(mpts[i]+"/bias"+postfix);
	  ffunc=biasH->GetFunction("gaus");
	  float avgBias=ffunc->GetParameter(1);
	  float avgBiasErr=ffunc->GetParError(1);   
	  ipt=gr->GetN();
	  gr->SetPoint(ipt,trueMass,avgBias);
	  gr->SetPointError(ipt,0,avgBiasErr);
	
	  if(trueMass==172.5) 
	    { 
	      TH1D *h=(TH1D *)fin->Get(mpts[i]+"/statunc"+postfix); 
	      h->SetDirectory(0);
	      h->SetLineColor(icat==nCategs? 1 : icat+2); 
	      h->SetMarkerColor(icat==nCategs? 1 : icat+2); 
	      h->SetMarkerStyle(1);
	      h->SetLineWidth(2);
	      statUncResults.Add(h);
	    }  
	}



      TH1D *pullH=(TH1D *)fin->Get(mpts[i]+"/pull");
      pullH->SetDirectory(0);
      pullResults.Add(pullH);
    }
  fin->Close();


  //plot results
  setStyle();
  gStyle->SetOptFit(0);

  TCanvas *c=new TCanvas("linc","linc",600,600);
  c->SetWindowSize(600,600);
  c->Divide(1,2);
  TPad *p=(TPad *)c->cd(1);
  p->SetPad(0,0.3,1.0,1.0);
  lingr->Draw("ap");
  lingr->GetXaxis()->SetTitle("Generated m_{top} [GeV/c^{2}]");
  lingr->GetXaxis()->SetTitle("Fitted m_{top} [GeV/c^{2}]");

  TLine *l=new TLine(lingr->GetXaxis()->GetXmin(),lingr->GetYaxis()->GetXmin(),
		     lingr->GetXaxis()->GetXmax(),lingr->GetYaxis()->GetXmax());
  l->SetLineColor(kGray);
  l->SetLineStyle(6);
  l->Draw("same");
  TLegend *leg = p->BuildLegend();
  formatForCmsPublic(p,leg,"CMS simulation",3);
  leg->Delete();


  p=(TPad *)c->cd(2);
  p->SetPad(0,0.0,1.0,0.28);
  p->SetTopMargin(0);
  p->SetBottomMargin(0.5);
  float yscale = (1.0-0.3)/(0.28-0);
  biasgr->Draw("ap");
  biasgr->GetXaxis()->SetTitle("Generated m_{top} [GeV/c^{2}]");
  biasgr->GetYaxis()->SetTitle("Bias");
  biasgr->GetYaxis()->SetRangeUser(-5.2,5.2);
  biasgr->GetXaxis()->SetTitleOffset(0.85);
  biasgr->GetXaxis()->SetLabelSize(0.04 * yscale);
  biasgr->GetXaxis()->SetTitleSize(0.05 * yscale);
  biasgr->GetXaxis()->SetTickLength( 0.03 * yscale );
  biasgr->GetYaxis()->SetTitleOffset(0.5);
  biasgr->GetYaxis()->SetLabelSize(0.04 * yscale);
  biasgr->GetYaxis()->SetTitleSize(0.04 * yscale);
  biasgr->Fit("pol1");
  float a=biasgr->GetFunction("pol1")->GetParameter(1);
  float b=biasgr->GetFunction("pol1")->GetParameter(0);
  cout << "[Inclusive bias correction] resslope:" << (a+1.) << " resbias:" << b << endl;
  c->Modified();
  c->Update();
  c->SaveAs(outDir+"/TopMassCalibrationResults.C");
  c->SaveAs(outDir+"/TopMassCalibrationResults.png");
  c->SaveAs(outDir+"/TopMassCalibrationResults.pdf");
  //  delete c;
  
  //biases per category
  c=new TCanvas("biasc","biasc",600,600);
  c->SetWindowSize(600,600);
  c->Divide(1,nCategs+1);
  c->cd(nCategs+1)->Delete();
  float ystep=(1.0)/(nCategs+1);
  for(int icat=0; icat<nCategs; icat++) 
    {
      TGraphErrors *gr=categbiasgr[icat];
      gr->SetMarkerStyle(20);
      TPad *p=(TPad *)c->cd(icat+1);
      p->SetFillStyle(1001);
      p->SetFillColor(0);
      float ymin=(icat+1)*ystep-0.05;
      float ymax=(icat+2)*ystep-0.05;
      if(icat==0)
	{
	  ymin=icat*ystep;
	}
      p->SetPad(0,ymin,1.0,ymax);
      p->SetTopMargin(0.0);
      p->SetBottomMargin(icat==0?0.4:0.0);
      gr->Draw("ap");
      float yscale = (0.95-ystep)/(ystep-0.0);
      if(icat==0) yscale=(0.95-2*ystep)/(1.2*ystep-0);
      gr->GetXaxis()->SetTitle("Generated m_{top} [GeV/c^{2}]");
      gr->GetYaxis()->SetTitle("Bias");
      gr->GetYaxis()->SetRangeUser(-5.2,5.2);
      gr->GetXaxis()->SetTitleOffset(0.85);
      gr->GetXaxis()->SetLabelSize(0.04 * yscale);
      gr->GetXaxis()->SetTitleSize(0.05 * yscale);
      gr->GetXaxis()->SetTickLength( 0.03 * yscale );
      gr->GetYaxis()->SetTitleOffset(icat==0?0.5:0.3);
      gr->GetYaxis()->SetLabelSize(0.04 * yscale);
      gr->GetYaxis()->SetTitleSize(0.04 * yscale);
      gr->Fit("pol1");
      float a=gr->GetFunction("pol1")->GetParameter(1);
      float b=gr->GetFunction("pol1")->GetParameter(0);
      cout << "[Bias correction #" << icat+1 << "] resslope:" << (a+1.) << " resbias:" << b << endl;
      
      //prepare label
      TPaveText *pave = new TPaveText(0.65,0.75,0.9,0.92,"NDC");
      pave->SetTextFont(42);
      pave->SetFillStyle(0);
      pave->SetBorderSize(0);
      pave->AddText( labels[icat] )->SetTextAlign(11);
      pave->Draw();

      if(icat==nCategs-1)
	{
	  pave = new TPaveText(0.2,0.75,0.45,0.92,"NDC");
	  pave->SetTextFont(42);
	  pave->SetFillStyle(0);
	  pave->SetBorderSize(0);
	  pave->AddText( "CMS simulation" )->SetTextAlign(11);
	  pave->Draw();
	}
    }
  c->Modified();
  c->Update();
  c->SaveAs(outDir+"/TopMassCalibrationBiases.C");
  c->SaveAs(outDir+"/TopMassCalibrationBiases.png");
  c->SaveAs(outDir+"/TopMassCalibrationBiases.pdf");
  //delete c;

  //stat uncertainty
  c=new TCanvas("statc","statc",600,600);
  c->SetWindowSize(600,600);
  for(int i=0; i<statUncResults.GetEntriesFast(); i++)
    {
      int idx=statUncResults.GetEntriesFast()-1-i;
      statUncResults.At(idx)->Draw(i==0?"hist":"histsame");
      ((TH1*)statUncResults.At(idx))->SetTitle( i==0 ? "combined" :labels[idx] );
    }
  leg = c->BuildLegend();
  formatForCmsPublic(c,leg,"CMS simulation",3);  c->Modified();
  c->Update();
  c->SaveAs(outDir+"/TopMassCalibrationStat.C");
  c->SaveAs(outDir+"/TopMassCalibrationStat.png");
  c->SaveAs(outDir+"/TopMassCalibrationStat.pdf");
  //delete c;
  
  //pulls
  TGraphErrors *avgPull = new TGraphErrors;   avgPull->SetMarkerStyle(20);
  TGraphErrors *sigmaPull = new TGraphErrors; sigmaPull->SetMarkerStyle(20);
  c=new TCanvas("pullsc","pullsc",1600,(pullResults.GetEntriesFast()/4+1)*400);
  c->SetWindowSize(1600,(pullResults.GetEntriesFast()/4+1)*400);
  c->Divide(4,pullResults.GetEntriesFast()/4+1);
  for(int i=0; i<pullResults.GetEntriesFast(); i++)
    {
      TPad *p=(TPad *)c->cd(i+1);
      TH1 *pullH=(TH1 *)pullResults.At(i);
      pullH->Draw("e1");

      Float_t trueMass=mpts[i].Atof();
      TF1 *gaus=pullH->GetFunction("gaus");
      avgPull->SetPoint(i,trueMass,gaus->GetParameter(1));
      avgPull->SetPointError(i,0,gaus->GetParError(1));
      sigmaPull->SetPoint(i,trueMass,gaus->GetParameter(2));
      sigmaPull->SetPointError(i,0,gaus->GetParError(2));

      if(i==0)
	{
	  leg = p->BuildLegend();
	  formatForCmsPublic(p,leg,"CMS simulation",3);
	  leg->Delete();
	}
      
      TPaveText *pave = new TPaveText(0.2,0.8,0.35,0.88,"NDC");
      pave->SetTextFont(42);
      pave->SetFillStyle(0);
      pave->SetBorderSize(0);
      pave->AddText( mpts[i] )->SetTextAlign(11);
      pave->Draw();

    }
  c->Modified();
  c->Update();
  c->SaveAs(outDir+"/TopMassCalibrationPulls.C");
  c->SaveAs(outDir+"/TopMassCalibrationPulls.png");
  c->SaveAs(outDir+"/TopMassCalibrationPulls.pdf");
  //delete c;

  //pull momenta
  c=new TCanvas("pullsmomc","pullsmomc",600,600);
  c->SetWindowSize(600,600);
  c->Divide(1,2);
  p=(TPad *)c->cd(1);
  avgPull->Draw("ap");
  avgPull->GetYaxis()->SetRangeUser(-3,3);
  avgPull->GetXaxis()->SetTitle("Generated m_{top} [GeV/c^{2}]");
  avgPull->GetYaxis()->SetTitle("Average pull");
  leg = p->BuildLegend();
  formatForCmsPublic(p,leg,"CMS simulation",3);
  leg->Delete();
  c->cd(2);
  sigmaPull->GetYaxis()->SetRangeUser(0.9,1.2);
  sigmaPull->Draw("ap");
  sigmaPull->GetXaxis()->SetTitle("Generated m_{top} [GeV/c^{2}]");
  sigmaPull->GetYaxis()->SetTitle("Pull width");
  c->Modified();
  c->Update();
  c->SaveAs(outDir+"/TopMassCalibrationPullsMomenta.C");
  c->SaveAs(outDir+"/TopMassCalibrationPullsMomenta.png");
  c->SaveAs(outDir+"/TopMassCalibrationPullsMomenta.pdf");
  //delete c;




}
