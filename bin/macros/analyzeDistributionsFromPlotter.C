/*
  gSystem->Load("libCMGToolsHtoZZ2l2nu.so");
  .L analyzeDistributionsFromPlotter.C+
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

//
TGraphAsymmErrors *getEfficiencyFrom(TH1F *h)
{
  if(h==0) return 0;

  // compute efficiency: N(x>cut) / N as function of x
  // we assume overflows / underflows are included
  
  TGraphAsymmErrors *effgr = new TGraphAsymmErrors;
  effgr->SetName(h->GetName()+TString("eff"));
  effgr->SetTitle(TString("#varepsilon_{")+h->GetTitle()+TString("}"));
  effgr->SetMarkerStyle(h->GetMarkerStyle());
  effgr->SetMarkerSize(1.3);
  effgr->SetMarkerColor(h->GetLineColor());
  effgr->SetLineColor(h->GetLineColor());
  effgr->SetFillStyle(0);
  int nbins=h->GetXaxis()->GetNbins();
  Double_t tot=h->Integral();
  for(int ibin=1; ibin<=nbins; ibin++)
    {
      Double_t cut = h->GetBinLowEdge(ibin);
      Double_t err;
      Double_t tot_cut=h->IntegralAndError(ibin,nbins,err);
      effgr->SetPoint(ibin-1,cut,tot_cut/tot);
      effgr->SetPointError(ibin-1,0,0,err/tot,err/tot);
    }
  
  return effgr;
}


//
TObjArray getDistributionFromPlotter(TString plot,TString baseURL="~/scratch0/top-nosyst/plotter.root")
{
  using namespace std;

  setStyle();
  
  //the relevant processes
  TString procs[]={"Di-bosons","Single top", "W+jets", "Z-#gamma^{*}+jets#rightarrow ll", "other t#bar{t}", "t#bar{t} dileptons", "data"};
  const size_t nprocs=sizeof(procs)/sizeof(TString);

  //containers for the histos
  TList *data     = new TList;
  TH1F *totalData=0;
  TList *mc       = new TList;
  TH1F *totalMC=0;
  TH1F *ttbarMC=0;
  TList *spimpose = new TList;
  
  //open file with histograms
  TFile *f=TFile::Open(baseURL);
  for(size_t iproc=0; iproc<nprocs; iproc++)
    {
      TH1F *histo = (TH1F *) f->Get(procs[iproc]+"/"+plot);
      if(histo==0) { cout << "[Warning] " << plot << " not found for " << procs[iproc] << " ...skipping" << endl; continue; }
      histo->SetDirectory(0);
      if(procs[iproc].Contains("data")) 
	{
	  data->Add(histo);
	  if(totalData==0) { totalData = (TH1F *) histo->Clone(plot+"totaldata"); totalData->SetDirectory(0); }
	  else             { totalData->Add(histo); }
	}
      else
	{
	  mc->Add(histo);
	  if(totalMC==0) { totalMC = (TH1F *) histo->Clone(plot+"totalmc"); totalMC->SetDirectory(0); }
	  else           { totalMC->Add(histo); }
	  if(procs[iproc]=="t#bar{t} dileptons") { ttbarMC = (TH1F *)histo->Clone(plot+"ttbar"); ttbarMC->SetDirectory(0); }
	}
    }
  f->Close();

  cout << "Drawing now" << endl;
  
  //draw
  TCanvas *cnv = getNewCanvas(plot+"c",plot+"c",false);
  cnv->Clear();
  cnv->SetWindowSize(600,600);
  cnv->cd(); 
  TLegend *leg=showPlotsAndMCtoDataComparison(cnv,*mc,*spimpose,*data,false);
  formatForCmsPublic(cnv,leg,"CMS preliminary", 4);
  cnv->SaveAs(plot+".C");

  cout << "Stat comparison for " << plot << endl
       << "Sample \t Average \t\t RMS " << endl
       << "----------------------------------------" << endl;
  if(totalData) cout << "Data \t " << totalData->GetMean() << " +/- " << totalData->GetMeanError() << " \t " << totalData->GetRMS() << " +/- " << totalData->GetRMSError() << endl;
  if(totalMC)   cout << "MC \t "   << totalMC->GetMean() << " +/- " << totalMC->GetMeanError() << " \t " << totalMC->GetRMS() << " +/- " << totalMC->GetRMSError() << endl;
  if(ttbarMC)   cout << "Signal \t "   << ttbarMC->GetMean() << " +/- " << ttbarMC->GetMeanError() << " \t " << ttbarMC->GetRMS() << " +/- " << ttbarMC->GetRMSError() << endl;

  TObjArray res;
  res.Add(totalData);
  res.Add(totalMC);
  res.Add(ttbarMC);
  return res;
}


//
void determineWorkingPoint(TString algo="csv",TString baseURL="~/scratch0/top-nosyst/plotter.root")
{
  float workPoint(0.244);
  float sfb(1.020),    sfberr(0.04);      //from BTV-11-003
  float sflight(1.08), sflighterr(0.13);  //from BTV-11-002
  //  float sfb(0.99),        sfberr(0.099);      //from BTV-11-001
  //  float sflight(1.07882), sflighterr(0.244);  //from BTV-11-001
  if(algo=="tche")
    {
      workPoint=1.7;
      sfb=0.95; sfberr=0.095;
      sflight=1.08018; sflighterr=0.1125;
    }

  TObjArray bjets=getDistributionFromPlotter(algo+"b",baseURL);
  TObjArray lightjets=getDistributionFromPlotter(algo+"light",baseURL);
  
  TH1F *bDisc=(TH1F *) bjets.At(1);
  TH1F *lightDisc=(TH1F *) lightjets.At(1);
  
  setStyle();
  gStyle->SetOptFit(0);

  TCanvas *cnv = getNewCanvas("c","c",false);
  cnv->Clear();
  cnv->SetCanvasSize(1200,1200);
  cnv->SetWindowSize(1200,1200);
  cnv->Divide(2,2);

  //draw the discriminator
  TPad *p = (TPad *)cnv->cd(1); 
  bDisc->SetTitle("b");
  bDisc->SetLineColor(1);
  bDisc->SetMarkerColor(1);
  bDisc->SetMarkerStyle(20);
  bDisc->SetFillStyle(0);
  bDisc->DrawNormalized("histe1");
  lightDisc->SetLineColor(1);
  lightDisc->SetMarkerColor(1);
  lightDisc->SetMarkerStyle(24);
  lightDisc->SetFillStyle(0);
  lightDisc->SetTitle("udcsg");
  lightDisc->DrawNormalized("histe1same");
  TLegend *leg=p->BuildLegend();
  formatForCmsPublic(p,leg,"CMS simulation",2);

  //draw the b/light efficiencies
  p=(TPad *)cnv->cd(2);
  p->SetLogy();
  TGraphAsymmErrors *bEff     = getEfficiencyFrom(bDisc);
  TGraphAsymmErrors *lightEff = getEfficiencyFrom(lightDisc);
  bEff->SetMarkerStyle(20);
  bEff->SetFillStyle(0);
  bEff->Draw("ap");
  bEff->GetXaxis()->SetTitle(bDisc->GetXaxis()->GetTitle());
  bEff->GetYaxis()->SetTitle("Efficiency");
  lightEff->SetMarkerStyle(24);
  lightEff->SetFillStyle(0);
  lightEff->Draw("p");

  //draw relatively to a given working point
  p=(TPad *)cnv->cd(3);
  Double_t baseBEff=bEff->Eval(workPoint);
  TGraphAsymmErrors *relBEff=new TGraphAsymmErrors;
  relBEff->SetMarkerStyle(20);
  relBEff->SetFillStyle(0);
  Double_t baseLightEff=lightEff->Eval(workPoint);
  TGraphAsymmErrors *relLightEff=new TGraphAsymmErrors;
  relLightEff->SetMarkerStyle(24);
  relLightEff->SetFillStyle(0);
  for(int ip=0; ip<bEff->GetN(); ip++)
    {
      Double_t cut, y,ey; 
      bEff->GetPoint(ip,cut,y);     ey = bEff->GetErrorY(ip);
      Double_t relEff(y/baseBEff);
      if(relEff<sfb+2*sfberr && relEff>sfb-2*sfberr)
	{
	  int ipt=relBEff->GetN();
	  relBEff->SetPoint(ipt,cut,relEff);
	  relBEff->SetPointError(ipt,0,0,ey/baseBEff,ey/baseBEff);
	}

      lightEff->GetPoint(ip,cut,y);     ey = lightEff->GetErrorY(ip);
      relEff=y/baseLightEff;
      //      if(relEff<sflight+7*sflighterr && relEff>sflight-7*sflighterr)
      if(relEff<sflight+3*sflighterr && relEff>sflight-3*sflighterr)
	{
	  int ipt=relLightEff->GetN();
	  relLightEff->SetPoint(ipt,cut,relEff);
	  relLightEff->SetPointError(ipt,0,0,ey/baseLightEff,ey/baseLightEff);
	}
    }

  relLightEff->Draw("ap");
  relLightEff->GetXaxis()->SetTitle( bDisc->GetXaxis()->GetTitle() );
  relLightEff->GetYaxis()->SetTitle( "#varepsilon/#varepsilon_{0}" );
  relLightEff->Fit("expo","Q+");
  TF1 *ffunc=relLightEff->GetFunction("expo");
  float newLightCut=(TMath::Log(sflight)-ffunc->GetParameter(0))/ffunc->GetParameter(1);
  float newLightCutErrPlus=(TMath::Log(sflight+sflighterr)-ffunc->GetParameter(0))/ffunc->GetParameter(1)-newLightCut;
  float newLightCutErrMinus=(TMath::Log(sflight-sflighterr)-ffunc->GetParameter(0))/ffunc->GetParameter(1)-newLightCut;
  TArrow *lightArrow = new TArrow(newLightCut, sflight-3*sflighterr, newLightCut, sflight-3*sflighterr*0.8, 0.02, "<|");
  lightArrow->SetLineColor(kGray+2);
  lightArrow->SetFillColor(kGray+2);
  lightArrow->Draw("SAME <|");

  relBEff->Draw("p");
  relBEff->Fit("expo","Q+");
  ffunc=relBEff->GetFunction("expo");
  float newBCut=(TMath::Log(sfb)-ffunc->GetParameter(0))/ffunc->GetParameter(1);
  float newBCutErrPlus=(TMath::Log(sfb+sfberr)-ffunc->GetParameter(0))/ffunc->GetParameter(1)-newBCut;
  float newBCutErrMinus=(TMath::Log(sfb-sfberr)-ffunc->GetParameter(0))/ffunc->GetParameter(1)-newBCut;
  cout <<  sfb << " " << sfb+sfberr << " " << sfb-sfberr << endl;
  TArrow *bArrow    = new TArrow(newBCut, sflight-3*sflighterr, newBCut, sflight-3*sflighterr*0.8, 0.02, "<|");
  bArrow->Draw("SAME <|");

  //draw epsilon_b vs epsilon_q
  p=(TPad *)cnv->cd(4);
  p->SetLogy();

  TGraphAsymmErrors *perf=new TGraphAsymmErrors;
  perf->SetName("algoperformance");
  perf->SetMarkerStyle(20);
  perf->SetFillStyle(0);
  for(int ip=0; ip<bEff->GetN(); ip++)
    {
      Double_t cut;
      Double_t x,ex; bEff->GetPoint(ip,cut,x);     ex = bEff->GetErrorY(ip);
      Double_t y,ey; lightEff->GetPoint(ip,cut,y); ey = lightEff->GetErrorY(ip);
      perf->SetPoint(ip,x,y);
      perf->SetPointError(ip,ex,ex,ey,ey);
    }
  perf->Draw("ap");
  perf->GetXaxis()->SetTitle( bEff->GetTitle() );
  perf->GetYaxis()->SetTitle( lightEff->GetTitle() );

  //new performance expected after applying the new cuts
  TGraphAsymmErrors *newperf=new TGraphAsymmErrors;
  newperf->SetName("newalgoperformance");
  newperf->SetMarkerStyle(24);
  newperf->SetFillStyle(0);
  newperf->SetLineWidth(2);
  newperf->SetLineColor(kRed);
  newperf->SetLineColor(kRed);
  newperf->SetPoint(0,baseBEff*sfb,baseLightEff*sflight);
  newperf->SetPointError(0,baseBEff*sfberr,baseBEff*sfberr,baseLightEff*sflighterr,baseLightEff*sflighterr);
  newperf->Draw("p");

  cnv->Modified();
  cnv->Update();
  cnv->SaveAs("discFlavor.C");  
  cnv->SaveAs("discFlavor.png");  


  cout << "[determineWorkingPoint]" << endl
       << "To emulate the measured scale-factors you can use the following new cuts per jet flavor" << endl
       << "SF-b : "     << newBCut << " +" << newBCutErrPlus << " " << newBCutErrMinus << endl
       << "SF-light : " << newLightCut << " +" << newLightCutErrPlus << " " << newLightCutErrMinus << endl;
}


