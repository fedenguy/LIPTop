/*
  root -l fitDYdistribution.C
*/

{
#include "TString.h"
#include "TH1F.h"
#include "TMath.h"
#include "TFile.h"
#include "TObjArray.h"
  
  gSystem->Load("libCMGToolsHtoZZ2l2nu.so");
  setStyle();

  using namespace RooFit;
  using namespace std;

  TString baseURL="${HOME}/scratch0/top/";
  gSystem->ExpandPathName(baseURL);
  TString stdUrl=baseURL+"/std_plotter.root";

  //histograms of interest

  //dy->tautau
  TString dyReplacementUrl=baseURL+"/syst_plotter.root";
  TString title("Z#rightarrow #tau#tau");
  TString histCompare[]={"emu_mtsum"};//"emu_leadlepton","emu_subleadlepton","emu_ptsum", "emu_dilmass","emu_met","emu_mtsum"};
  int histoFit=0;
  bool rebin(false);
  Float_t fitmin(0), fitmax(150);
  TString replacementHisto("data/emu_mtsum");


  //dy->ee
  //TString dyReplacementUrl=baseURL+"/std_plotter.root";
  //TString title("Z#rightarrow ee");
  //TString histCompare[]={"ee_dilarccosine","ee_dilarccosinelowmet"};
  //Float_t fitmin(0.), fitmax(3.3);
  //bool rebin(false);
  //int histoFit=0;
  //TString replacementHisto("data/ee_lowmetdilarccosine");
  
  //dy->mumu
  //TString dyReplacementUrl=baseURL+"/std_plotter.root";
  //TString title("Z#rightarrow #mu#mu");
  //TString histCompare[]={"mumu_dilarccosine","mumu_dilarccosinelowmet"};
  //Float_t fitmin(0.), fitmax(3.3);
  //bool rebin(false);
  //int histoFit=0;
  //TString replacementHisto("data/mumu_lowmetdilarccosine");

  const size_t nHistos=sizeof(histCompare)/sizeof(TString);

  //
  //get histograms from files
  //
  TFile *dyReplacementFile=TFile::Open(dyReplacementUrl);
  TObjArray dyReplacementHistos;
  TFile *stdFile=TFile::Open(stdUrl);
 
  TObjArray mcHistos,dymcHistos,dataHistos;
  TString procs[]={"Di-bosons","Single top", "W+jets", "Z-#gamma^{*}+jets#rightarrow ll", "other t#bar{t}", "t#bar{t} dileptons", "data"};
  const size_t nprocs=sizeof(procs)/sizeof(TString);
  for(size_t ihisto=0; ihisto<nHistos; ihisto++)
    {
      //save copy
      TH1F *histo=(TH1F *) dyReplacementFile->Get(replacementHisto);
 
      if(rebin) histo->Rebin(2);
      histo->SetDirectory(0);
      histo->SetName("dymodel");
      histo->SetTitle(title);
      dyReplacementHistos.Add(histo);

      for(size_t iproc=0; iproc<nprocs; iproc++)
	{
	  TString theHisto=histCompare[ihisto];
	  if(ihisto==histoFit && !procs[iproc].Contains("data"))  //scale variations
	    {
	      //theHisto += "pudown";
	      theHisto += "puup";
	      //theHisto += "jesup";
	      //theHisto += "jesdown";
	      //theHisto += "jer";
	    }
	  TH1F *histo = (TH1F *) stdFile->Get(procs[iproc]+"/"+theHisto);
	  if(rebin) histo->Rebin(2);
	  //if(!procs[iproc].Contains("data")) histo.Scale(1.0-0.045);   //lumi
	  //if(!procs[iproc].Contains("data")) histo.Scale(1.0+0.045);
	  //if(!procs[iproc].Contains("data")) histo.Scale(1.0+0.04);    //selection
	  //if(!procs[iproc].Contains("data")) histo.Scale(1.0-0.04);
	  //if(procs[iproc].Contains("t#bar{t}")) histo.Scale(1.0-0.04);  //ttbar xsec
	  //if(procs[iproc].Contains("t#bar{t}")) histo.Scale(1.0+0.04);  

	  //check category
	  TString key("Other processes"), keyName("other");
	  TObjArray *container = &mcHistos;
	  if(procs[iproc].Contains("Z-#gamma"))   
	    {
	      keyName="dymc";
	      key=title;
	      container = &dymcHistos;
	    }
	  else if(procs[iproc].Contains("data"))  
	    { 
	      keyName="data";
	      key="data";
	      container = &dataHistos;
	    }
	  
	  //add the current histo
	  if(container->GetEntriesFast() <= ihisto )
	    {
	      histo = (TH1F *) histo->Clone(keyName);
	      histo->SetDirectory(0);
	      histo->SetTitle(key);
	      container->Add(histo);
	    }
	  else
	    {
	      ((TH1F *)container->At(ihisto))->Add(histo);
	    }
	}
    }
  dyReplacementFile->Close();
  stdFile->Close();
  
  cout << dyReplacementHistos.GetEntriesFast() << " kinematic distributions retrieved from tau replacement" << endl;

  //compare data and MC
  TCanvas *cnv = getNewCanvas("compc","compc",false);
  cnv->Clear();
  cnv->SetWindowSize(1200,800);
  cnv->Divide(3,2);
  for(size_t ihisto=0; ihisto<nHistos; ihisto++)
    {
      TPad *p=(TPad *) cnv->cd(ihisto+1);
      TH1F *kindy=(TH1F *) dyReplacementHistos.At(ihisto);
      TH1F *kindymc=(TH1F *) dymcHistos.At(ihisto);
      kindymc->GetYaxis()->SetTitle("Events (a.u.)");
      kindymc->SetFillStyle(3001);
      kindymc->SetFillColor(1);
      kindymc->SetMarkerColor(1);
      kindymc->SetMarkerStyle(1);
      kindymc->DrawNormalized("e4");
      kindy->DrawNormalized("e1same");
      
      if(ihisto==0)
	{
	  TLegend *leg=p->BuildLegend();
	  leg->SetBorderSize(0);
	  leg->SetHeader("CMS preliminary");
	  leg->SetFillStyle(0);
	  leg->SetFillColor(0);
	  leg->SetTextFont(42);
	}
    }
  cnv->SaveAs("dydistcomp.C");

  //
  // Now fit
  //
  cnv = getNewCanvas("fitc","fitc",false);
  cnv->Clear();
  cnv->SetWindowSize(600,600);
  TH1* kindata = (TH1F *) dataHistos.At(histoFit);
  TH1F *kindy=(TH1F *) dyReplacementHistos.At(histoFit);
  TH1F *kinmc=(TH1F *) mcHistos.At(histoFit);
  TH1F *kindymc=(TH1F *) dymcHistos.At(histoFit);
  
  RooRealVar x("x","x", fitmin, kindata->GetXaxis()->GetXmax());
  x->setBins(kindata->GetXaxis()->GetNbins());
      
  RooDataHist* mcdyTemplate = new RooDataHist("mcdyTemplate", "mcdyTemplate", RooArgList(x), kindymc );
  
  std::cout << " Preparing data template" << std::endl;
  RooDataHist* dataTemplate = new RooDataHist("dataTemplate", "dataTemplate", RooArgList(x), kindy );
  RooHistPdf modelDataTemplate("modelDataTemplate", "modelDataTemplate", RooArgSet(x), *dataTemplate);
  Double_t dymcyields=kindymc->Integral();
  RooRealVar ndyexp("<N>_{"+title+"}","dyyieldsexp",dymcyields);
  RooRealVar ndysf("SF_{DY}","dyyieldssfactor",1.0,0.,5.0);
  RooFormulaVar ndy("N_{"+title+"}","@0*@1",RooArgSet(ndyexp,ndysf));
  
  std::cout << " Preparing mc template" << std::endl;
  RooDataHist *mcTemplate = new RooDataHist("mcTemplate", "mcTemplate", RooArgList(x), kinmc);
  RooHistPdf modelMcTemplate("modelMcTemplate", "modelMcTemplate", RooArgSet(x), *mcTemplate);
  Double_t otheryields(0), otheryields_err(0);
  otheryields=kinmc->IntegralAndError(0,kinmc->GetXaxis()->GetNbins()+1,otheryields_err);
  RooRealVar nother("N_{other}","otheryields",otheryields,TMath::Max(otheryields-2*otheryields_err,0.),otheryields+2*otheryields_err);
  RooRealVar nother_mean("meanother","meanother",otheryields);  
  RooRealVar nother_sigma("sigmaother","sigmaother",otheryields_err);
  RooGaussian other_constraint("otherconstraintpdf","otherconstraintpdf", nother, nother_mean, nother_sigma);
  
  std::cout << " Preparing data to fit" << std::endl;
  RooDataHist* sumData = new RooDataHist("sumData", "sumData", RooArgList(x), kindata);
  
  std::cout << " Fitting now ..." << std::endl;
  RooAddPdf shapeModel("shapemodel","signal+background",RooArgList(modelDataTemplate,modelMcTemplate),RooArgList(ndy,nother));
  //  RooProdPdf model("model","(signal+background)*evconstraint*bkgconstraint",RooArgSet(other_constraint,shapeModel));
  //  model.fitTo(*sumData,Extended(kTRUE), Constrain(nother),Minos(),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(fitmin,fitmax));
  RooProdPdf model("model","(signal+background)*evconstraint*bkgconstraint",RooArgSet(shapeModel));
  model.fitTo(*sumData,Minos(),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(fitmin,fitmax));
  
  RooPlot *genericFrame = x.frame();
  genericFrame->GetXaxis()->SetTitle( kindata->GetXaxis()->GetTitle() );
  genericFrame->GetYaxis()->SetTitle( "Events" );
  sumData->plotOn(genericFrame);
  // mcdyTemplate->plotOn(genericFrame,RooFit::FillStyle(1001));
  // model.plotOn(genericFrame,Range(fitmin,kindata->GetXaxis()->GetXmax()),RooFit::LineStyle(kDashed));
  // model.plotOn(genericFrame,RooFit::Components(modelMcTemplate),Range(fitmin,kindata->GetXaxis()->GetXmax()));

  model.plotOn(genericFrame,Range(fitmin,kindata->GetXaxis()->GetXmax()));
  model.plotOn(genericFrame,RooFit::Components(modelDataTemplate),Range(fitmin,kindata->GetXaxis()->GetXmax()),RooFit::LineStyle(kDashed));
  genericFrame->Draw();
  
  //prepare label
  TPaveText *pave = new TPaveText(0.7,0.75,0.95,0.93,"NDC");
  pave->SetFillStyle(0);
  pave->SetBorderSize(0);
  pave->AddText("CMS preliminary, #int L=2.2 fb^{-1}, #sqrt{s}=7 TeV");
  char buf[100];
  sprintf(buf,"SF_{DY}=%3.1f #pm %3.1f",ndysf.getVal(),ndysf.getError());
  pave->AddText(buf)->SetTextAlign(11);
  pave->Draw();
  pave->SetTextFont(42);
  
  //likelihoods
  TPad *npad = new TPad("llpad","ll", 0.6, 0.6, 0.9, 0.9);
  npad->Draw();
  npad->cd();
  RooNLLVar *nll = (RooNLLVar*) model.createNLL(*sumData,RooFit::CloneData(kFALSE),Extended(kTRUE),Constrain(nother),Range(fitmin,fitmax));
  RooMinuit minuit(*nll); 
  minuit.migrad();
  minuit.hesse();
  
  RooPlot *frame2 = ndysf.frame();
  nll->plotOn(frame2,ShiftToZero(),Name("ll"));
  frame2->GetXaxis()->SetTitle("SF_{DY}");
  frame2->GetXaxis()->SetTitleOffset(0.8);
  frame2->GetYaxis()->SetTitle("-log(L/L_{max})");
  frame2->GetYaxis()->SetTitleOffset(1);
  frame2->Draw();

  cnv->SaveAs("dydistfit.C");
  cnv->SaveAs("dydistfit.png");

}
