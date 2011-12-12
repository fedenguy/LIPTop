#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooSimPdfBuilder.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooNLLVar.h"
#include "RooProdPdf.h"
#include "RooProfileLL.h"

#include "TList.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TSystem.h"
#include "TString.h"
#include "TH1F.h"
#include "TMath.h"
#include "TFile.h"
#include "TObjArray.h"

#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include<vector>
#include<sstream>

using namespace RooFit;
using namespace std;

int main(int argc,char *argv[])
{
  stringstream report;

  if(argc<2)
    {
      cout << endl
	   << "fitDYcontribution path/std_plotter.root <path/syst_plotter.root>" << endl
	   << endl;
      return -1;
    }

  setStyle();

  TString stdUrl=argv[1];
  gSystem->ExpandPathName(stdUrl);
  TString dyReplacementUrl("");
  if(argc>2) 
    {
      dyReplacementUrl=argv[2];
      gSystem->ExpandPathName(dyReplacementUrl);
    }

  TString vars[]={"","puup","pudown","jesup","jesdown","jer"};
  const size_t nvars=sizeof(vars)/sizeof(TString);

  TString ch[]={"ee","mumu","emu"};
  const size_t nchs=sizeof(ch)/sizeof(TString);

  TString signalRegionHisto[]={"ee_dilarccosine","mumu_dilarccosine","emu_mtsum"};
  TString templateHisto[]={"ee_lowmetdilarccosine","mumu_lowmetdilarccosine","emu_mtsum"};
  TString templateTitle[]={"Z#rightarrow ee","Z#rightarrow #mu#mu","Z#rightarrow e#mu"};
  TString templateName[]={"dytoee","dytomumu","dytoemu"};
  Float_t fitmin[]={0.,0.,0.};
  Float_t fitmax[]={3.3,3.3,400.};
  bool rebin[]={true,true,false};

  for(size_t ich=0; ich<nchs; ich++)
    {
      
      //get the DY template
      if(ich==2 && dyReplacementUrl.IsNull()) continue;
      TFile *dyReplacementFile=TFile::Open(ich==2 ? dyReplacementUrl : stdUrl);
      TH1F *dyReplacementHisto=(TH1F *) dyReplacementFile->Get("data/"+templateHisto[ich]);
      dyReplacementHisto->SetDirectory(0); 
      if(rebin[ich]) dyReplacementHisto->Rebin(2);
      dyReplacementHisto->SetName(templateName[ich]+"model");
      dyReplacementHisto->SetTitle(templateTitle[ich]);
      dyReplacementFile->Close();
  
      //get the histograms for data, DY (simulation) and other processes (simulation)  in the signal region
      TFile *stdFile=TFile::Open(stdUrl);
      
      TH1F *dataHisto=(TH1F *) stdFile->Get("data/"+signalRegionHisto[ich])->Clone("data");
      dataHisto->SetDirectory(0);
      if(rebin[ich]) dataHisto->Rebin(2);
      dataHisto->SetTitle("data");
      
      std::vector<TH1F *> dyMCHisto, otherProcsHisto;
      TString procs[]={"Di-bosons","Single top", "W+jets", "Z-#gamma^{*}+jets#rightarrow ll", "other t#bar{t}", "t#bar{t} dileptons"};
      const size_t nprocs=sizeof(procs)/sizeof(TString);
      for(size_t ivar=0; ivar<nvars; ivar++)
	{
	  for(size_t iproc=0; iproc<nprocs; iproc++)
	    {
	      TH1F *histo = (TH1F *) stdFile->Get(procs[iproc]+"/"+signalRegionHisto[ich]+vars[ivar]);
	      if(rebin[ich]) histo->Rebin(2);
	      if(procs[iproc].Contains("Z-#gamma"))   
		{
		  dyMCHisto.push_back( (TH1F *)histo->Clone(templateName[ich]+vars[ivar]+"mc") );
		  dyMCHisto[ivar]->SetDirectory(0);
		  dyMCHisto[ivar]->SetTitle(templateTitle[ich]+" (MC)");
		}
	      else if(otherProcsHisto.size()==ivar)
		{
		  otherProcsHisto.push_back( (TH1F *) histo->Clone("otherprocs"+vars[ivar]) );
		  otherProcsHisto[ivar]->SetDirectory(0);
		  otherProcsHisto[ivar]->SetTitle("Other processes");
		}
	      else
		{
		  otherProcsHisto[ivar]->Add(histo);
		}
	    }
	}
      stdFile->Close();

      //
      // Now fit
      //
      RooRealVar x("x","x", fitmin[ich], dataHisto->GetXaxis()->GetXmax());
      
      RooDataHist* sumData = new RooDataHist("sumData", "sumData", RooArgList(x), dataHisto);
      float totalData=dataHisto->Integral(0,dataHisto->GetXaxis()->GetNbins()+1);
      
      RooDataHist* dataTemplate = new RooDataHist("dataTemplate", "dataTemplate", RooArgList(x), dyReplacementHisto );
      RooHistPdf modelDataTemplate("modelDataTemplate", "modelDataTemplate", RooArgSet(x), *dataTemplate);

      float dyExpected=dyMCHisto[0]->Integral(0,dyMCHisto[0]->GetXaxis()->GetNbins()+1);
      float totalOthers=otherProcsHisto[0]->Integral(0,otherProcsHisto[0]->GetXaxis()->GetNbins()+1);

      RooRealVar ndyexp("<N>_{DY-MC}","dyyieldsexp",dyExpected);
      //float maxdysf=(totalData-totalOthers)/dyExpected;
      RooRealVar ndysf("SF_{DY}","dyyieldssfactor",1.0,0.,5.0);
      RooFormulaVar ndy("N_{DY}","@0*@1",RooArgSet(ndyexp,ndysf));

      RooRealVar nother("N_{other}","otheryields",totalOthers,0,2*totalData);
      RooRealVar nother_mean("meanother","meanother",totalOthers);
      RooRealVar nother_sigma("sigmaother","sigmaother",0.06*totalOthers); //use lumi uncertainty (conservative 6%)
      RooGaussian other_constraint("otherconstraintpdf","otherconstraintpdf", nother, nother_mean, nother_sigma);
      
      report << "***** Fit results for " << templateTitle[ich] << " *********" << endl;
      for(size_t ivar=0; ivar<nvars; ivar++)
	{
	  RooDataHist *mcTemplate = new RooDataHist("mcTemplate", "mcTemplate", RooArgList(x), otherProcsHisto[ivar]);
	  RooHistPdf modelMcTemplate("modelMcTemplate", "modelMcTemplate", RooArgSet(x), *mcTemplate);
	  
	  RooAddPdf shapeModel("shapemodel","signal+background",RooArgSet(modelDataTemplate,modelMcTemplate),RooArgSet(ndy,nother));
	  RooProdPdf model("model","(signal+background)*otherconstrain",RooArgSet(shapeModel,other_constraint)); 
	  model.fitTo(*sumData,Minos(),Extended(kTRUE),Constrain(nother),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(fitmin[ich],fitmax[ich]));

	  report << "\t" << vars[ivar] << "\t" << ndy.getVal() << " +/- " << ndysf.getError()*ndyexp.getVal() << flush;
	  if(ivar==0)   report << "\t SF-DY: "<< ndysf.getVal() << " +/- " << ndysf.getError() << "\t\t (" << nother.getVal() << " +/- " << nother.getError() << " )" << flush;
	  report << endl;

	  //
	  // Display results
	  //
	  if(ivar>0) continue;
	  //compare data template and MC
	  TCanvas *cnv = getNewCanvas("compc","compc",false);
	  cnv->SetWindowSize(600,600);
	  cnv->cd();
	  dyMCHisto[0]->DrawNormalized("hist");
	  dyReplacementHisto->DrawNormalized("e1same");
	  TLegend *leg=cnv->BuildLegend();
	  formatForCmsPublic(cnv,leg,"CMS preliminary",2);
	  cnv->SaveAs(ch[ich]+"dydistcomp.C");
	  cnv->SaveAs(ch[ich]+"dydistcomp.png");
	  
	  
	  //fit results
	  cnv = getNewCanvas("fitc","fitc",false);
	  cnv->SetWindowSize(600,600);
	  cnv->cd();
	  
	  RooPlot *genericFrame = x.frame(Range(fitmin[ich],fitmax[ich]));
	  genericFrame->GetXaxis()->SetTitle( dataHisto->GetXaxis()->GetTitle() );
	  genericFrame->GetYaxis()->SetTitle( "Events" );
	  sumData->plotOn(genericFrame);
	  shapeModel.plotOn(genericFrame,Range(fitmin[ich],dataHisto->GetXaxis()->GetXmax()));
	  shapeModel.plotOn(genericFrame,RooFit::Components(modelDataTemplate),Range(fitmin[ich],dataHisto->GetXaxis()->GetXmax()),RooFit::LineStyle(kDashed));
	  genericFrame->Draw();
	  
	  TString caption("CMS preliminary\\");
	  char buf[100];
	  sprintf(buf,"SF_{DY}=%3.1f #pm %3.1f\\",ndysf.getVal(),ndysf.getError());
	  caption += buf;
	  sprintf(buf,"N_{DY}=%3.1f #pm %3.1f",ndy.getVal(),ndysf.getError()*ndyexp.getVal());
	  caption += buf;
	  leg=cnv->BuildLegend();
	  formatForCmsPublic(cnv,leg,caption,2);
	  
	  //display likelihood as inset
	  TPad *npad = new TPad("llpad","ll", 0.6, 0.6, 0.9, 0.9);
	  npad->Draw();
	  npad->cd();
	  RooNLLVar *nll = (RooNLLVar*) model.createNLL(*sumData,RooFit::CloneData(kFALSE),Extended(kTRUE),Constrain(nother),Range(fitmin[ich],fitmax[ich]));
	  RooAbsReal* profileLL = nll->createProfile(ndysf);
	  //RooMinuit minuit(*nll); 
	  //minuit.migrad();
	  //minuit.hesse();
	  
	  RooPlot *frame2 = ndysf.frame();
	  profileLL->plotOn(frame2,ShiftToZero(),Name("ll"));
	  frame2->GetXaxis()->SetTitle("SF_{DY}");
	  frame2->GetXaxis()->SetTitleOffset(0.8);
	  frame2->GetYaxis()->SetTitle("-log(L/L_{max})");
	  frame2->GetYaxis()->SetTitleOffset(1);
	  frame2->GetYaxis()->SetRangeUser(0,5);
	  frame2->Draw();
	  
	  cnv->Modified();
	  cnv->Update();
	  
	  cnv->SaveAs(ch[ich]+"dydistfit.C");
	  cnv->SaveAs(ch[ich]+"dydistfit.png");
	}
      report << "******************************************" << endl;
    }

  //all done... print report
  cout << report.str() << endl;
}
