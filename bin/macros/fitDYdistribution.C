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

  TString baseURL="~/scratch0/";

  //histograms of interest
  TString histCompare[]={"emu_leadlepton","emu_subleadlepton","emu_ptsum",
			 "emu_dilmass","emu_met","emu_mtsum"};
  const size_t nHistos=sizeof(histCompare)/sizeof(TString);
  int histoFit=5;

  //
  //get histograms from files
  //
  TString tautauReplacementUrl=baseURL+"taurep/plots/plotter.root";
  TFile *tautauReplacementFile=TFile::Open(tautauReplacementUrl);
  TObjArray tautauReplacementHistos;
  
  TString stdUrl=baseURL+"syst/plots/plotter.root";
  TFile *stdFile=TFile::Open(stdUrl);
 
  TObjArray mcHistos,dymcHistos,dataHistos;
  TString procs[]={"Di-bosons","Single top", "W+jets", "Z-#gamma^{*}+jets#rightarrow ll", "other t#bar{t}", "t#bar{t} dileptons", "data"};
  const size_t nprocs=sizeof(procs)/sizeof(TString);
  for(size_t ihisto=0; ihisto<nHistos; ihisto++)
    {
      //save copy
      TH1F *histo=(TH1F *) tautauReplacementFile->Get("data/"+histCompare[ihisto]);
      histo->SetDirectory(0);
      histo->SetName("tautaureplace");
      histo->SetTitle("Z#rightarrow #tau#tau (data)");
      tautauReplacementHistos.Add(histo);
      
      for(size_t iproc=0; iproc<nprocs; iproc++)
	{

	  TString theHisto=histCompare[ihisto];

	  //if(theHisto.Contains("mtsum") && !procs[iproc].Contains("data")) theHisto += "pudown";
	  TH1F *histo = (TH1F *) stdFile->Get(procs[iproc]+"/"+theHisto);
	  if(!procs[iproc].Contains("data")) histo.Scale(1.0+0.04);
	  //if(procs[iproc].Contains("t#bar{t}")) histo.Scale(1.0+0.04);

	  //check category
	  TString key("Other processes"), keyName("other");
	  TObjArray *container = &mcHistos;
	  if(procs[iproc].Contains("Z-#gamma"))   
	    {
	      keyName="dytautaumc";
	      key="Z#rightarrow #tau#tau (MC)";
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
  tautauReplacementFile->Close();
  stdFile->Close();
  
  cout << tautauReplacementHistos.GetEntriesFast() << " kinematic distributions retrieved from tau replacement" << endl;

  //compare data and MC
  TCanvas *cnv = getNewCanvas("compc","compc",false);
  cnv->Clear();
  cnv->SetWindowSize(1200,800);
  cnv->Divide(3,2);
  for(size_t ihisto=0; ihisto<nHistos; ihisto++)
    {
      TPad *p=(TPad *) cnv->cd(ihisto+1);
      TH1F *kintautau=(TH1F *) tautauReplacementHistos.At(ihisto);
      TH1F *kindymc=(TH1F *) dymcHistos.At(ihisto);
      kindymc->GetYaxis()->SetTitle("Events (a.u.)");
      kindymc->SetFillStyle(3001);
      kindymc->SetFillColor(1);
      kindymc->SetMarkerColor(1);
      kindymc->SetMarkerStyle(1);
      kindymc->DrawNormalized("e4");
      kintautau->DrawNormalized("e1same");

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

  //
  // Now fit
  //
  cnv = getNewCanvas("fitc","fitc",false);
  cnv->Clear();
  cnv->SetWindowSize(600,600);
  TH1* kindata = (TH1F *) dataHistos.At(histoFit);
  TH1F *kintautau=(TH1F *) tautauReplacementHistos.At(histoFit);
  TH1F *kinmc=(TH1F *) mcHistos.At(histoFit);
  TH1F *kindymc=(TH1F *) dymcHistos.At(histoFit);
  
  RooRealVar x("x","x", kindata->GetXaxis()->GetXmin(), kindata->GetXaxis()->GetXmax());
  x->setBins(kindata->GetXaxis()->GetNbins());
      
  RooDataHist* mcdyTemplate = new RooDataHist("mcdyTemplate", "mcdyTemplate", RooArgList(x), kindymc );
  
  std::cout << " Preparing data template" << std::endl;
  RooDataHist* dataTemplate = new RooDataHist("dataTemplate", "dataTemplate", RooArgList(x), kintautau );
  RooHistPdf modelDataTemplate("modelDataTemplate", "modelDataTemplate", RooArgSet(x), *dataTemplate);
  Double_t dymcyields=kindymc->Integral();
  RooRealVar ndytautauexp("<N>_{Z#rightarrow #tau#tau}","dytautauyieldsexp",dymcyields);
  RooRealVar ndytautausf("SF_{DY}","dytautauyieldssfactor",1.0,0.2,3.0);
  RooFormulaVar ndytautau("N_{Z#rightarrow #tau#tau}","@0*@1",RooArgSet(ndytautauexp,ndytautausf));
  
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
  RooAddPdf shapeModel("shapemodel","signal+background",RooArgList(modelDataTemplate,modelMcTemplate),RooArgList(ndytautau,nother));
  RooProdPdf model("model","(signal+background)*evconstraint*bkgconstraint",RooArgSet(other_constraint,shapeModel));
  model.fitTo(*sumData,Extended(kTRUE), Constrain(nother),Minos(),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(0,100));
  
  RooPlot *genericFrame = x.frame();
  genericFrame->GetXaxis()->SetTitle( kindata->GetXaxis()->GetTitle() );
  genericFrame->GetYaxis()->SetTitle( "Events" );
  sumData->plotOn(genericFrame);
  // mcdyTemplate->plotOn(genericFrame,RooFit::FillStyle(1001));
  model.plotOn(genericFrame,Range(kindata->GetXaxis()->GetXmin(), kindata->GetXaxis()->GetXmax()));
  model.plotOn(genericFrame,RooFit::Components(modelDataTemplate),RooFit::LineStyle(kDashed));
  genericFrame->Draw();
  
  //prepare label
  TPaveText *pave = new TPaveText(0.7,0.75,0.95,0.93,"NDC");
  pave->SetFillStyle(0);
  pave->SetBorderSize(0);
  pave->AddText("CMS preliminary");
  char buf[100];
  sprintf(buf,"SF_{DY}=%3.1f #pm %3.1f",ndytautausf.getVal(),ndytautausf.getError());
  pave->AddText(buf)->SetTextAlign(11);
  pave->Draw();
  pave->SetTextFont(42);
  
  //likelihoods
  TPad *npad = new TPad("llpad","ll", 0.6, 0.6, 0.9, 0.9);
  npad->Draw();
  npad->cd();
  RooNLLVar *nll = (RooNLLVar*) model.createNLL(*sumData,RooFit::CloneData(kFALSE),Extended(kTRUE),Constrain(nother),Range(0,100));
  RooMinuit minuit(*nll); 
  minuit.migrad();
  minuit.hesse();
  
  RooPlot *frame2 = ndytautausf.frame();
  nll->plotOn(frame2,ShiftToZero(),Name("ll"));
  frame2->GetXaxis()->SetTitle("SF_{DY}");
  frame2->GetXaxis()->SetTitleOffset(0.8);
  frame2->GetYaxis()->SetTitle("-log(L/L_{max})");
  frame2->GetYaxis()->SetTitleOffset(1);
  frame2->Draw();
}
