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

//
void printHelp()
{
  printf("--in      --> input file with standard control plots from showControlPlots\n");
  printf("--ttrep   --> input file with control plots from showControlPlots for the syst. json (optional)\n");
  printf("--syst    --> will run systematics also\n");
  printf("command line example: fitDYcontribution --in std_plotter.root --ttrep syst_plotter.root\n");
}

//
int main(int argc,char *argv[])
{

  bool doSyst(false);
  TString stdUrl(""),dyReplacementUrl("");
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos) { printHelp();    return 0; }
      if(arg.find("--in")!=string::npos)                  { stdUrl=argv[i+1];           gSystem->ExpandPathName(stdUrl);            i++;  printf("in      = %s\n", stdUrl.Data()); }
      if(arg.find("--ttrep")!=string::npos)               { dyReplacementUrl=argv[i+1]; gSystem->ExpandPathName(dyReplacementUrl);  i++;  printf("ttRep   = %s\n", dyReplacementUrl.Data()); }
      if(arg.find("--syst")!=string::npos)                { doSyst=true;                                                            i++;  printf("Will run systematics\n"); }
    }
  if(stdUrl=="") { printHelp(); return 0; }


  stringstream report;
  setStyle();

  //systematics to consider for the templates
  std::vector<TString> vars;
  vars.push_back("");
  if(doSyst)
    {
      vars.push_back("puup");
      vars.push_back("pudown");
      vars.push_back("jesup");
      vars.push_back("jesdown");
      vars.push_back("jer");
    }
  const size_t nvars=vars.size();

  TString ch[]               ={"ee1jets",                      "ee",                    "mumu1jets",                      "mumu",                    "emu"};
  size_t runNvars[]          ={1,                              nvars,                   1,                                 nvars,                    nvars};
  TString signalRegionHisto[]={"ee_eq1jetsdilarccosine",       "ee_dilarccosine",       "mumu_eq1jetsdilarccosine",       "mumu_dilarccosine",       "emu_mtsum"};
  TString templateHisto[]    ={"ee_eq1jetslowmetdilarccosine", "ee_lowmetdilarccosine", "mumu_eq1jetslowmetdilarccosine", "mumu_lowmetdilarccosine", "emu_mtsum"};
  TString templateTitle[]    ={"Z#rightarrow ee (=1 jets)",    "Z#rightarrow ee",       "Z#rightarrow #mu#mu (=1 jets)",  "Z#rightarrow #mu#mu",     "Z#rightarrow e#mu"};
  TString templateName[]     ={"dytoee",                       "dytoee",                "dytomumu",                       "dytomumu",                "dytoemu"};
  Float_t fitmin[]           ={0.,                             0.,                      0.,                               0.,                        0.};
  Float_t fitmax[]           ={3.3,                            3.3,                     3.3,                              3.3,                       400.};
  bool rebin[]               ={false,                          false,                   false,                            false,                     false};
  //  bool rebin[]               ={true,                           true,                    true,                             true,                      false};

  const size_t nchs=sizeof(ch)/sizeof(TString);
  for(size_t ich=0; ich<nchs; ich++)
    {
      //get the DY template
      bool getTemplateFromSyst(false);
      if(ch[ich]=="emu")
	{
	  if(dyReplacementUrl.IsNull()) continue;
	  getTemplateFromSyst=true;
	}
      TFile *dyReplacementFile=TFile::Open(getTemplateFromSyst ? dyReplacementUrl : stdUrl);
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
      
      std::vector<TH1F *> dyReplacementMCHisto, dyMCHisto, otherProcsHisto, mcSumHisto;
      TString procs[]={"Di-bosons","Single top", "W+jets", "Z-#gamma^{*}+jets#rightarrow ll", "other t#bar{t}", "t#bar{t} dileptons"};
      const size_t nprocs=sizeof(procs)/sizeof(TString);
      for(size_t ivar=0; ivar<runNvars[ich]; ivar++)
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

		  TH1F *repHistoMC = (TH1F *) stdFile->Get(procs[iproc]+"/"+templateHisto[ich]+vars[ivar]);
		  if(rebin[ich]) repHistoMC->Rebin(2);
		  dyReplacementMCHisto.push_back( (TH1F *) repHistoMC->Clone(templateName[ich]+vars[ivar]+"repmc") );
		  dyReplacementMCHisto[ivar]->SetDirectory(0);
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
	      
	      //total mc
	      if(iproc==0) {
		mcSumHisto.push_back( (TH1F *)histo->Clone("mcsum"+vars[ivar]) );
		mcSumHisto[ivar]->SetDirectory(0);
		mcSumHisto[ivar]->SetTitle("Total MC");
	      }
	      else
		{
		  mcSumHisto[ivar]->Add(histo);
		}
	    }
	}
      stdFile->Close();

      //
      // Now fit
      //
      RooRealVar x("x","x", fitmin[ich], dataHisto->GetXaxis()->GetXmax());
      
      //data
      RooDataHist* sumData = new RooDataHist("sumData", "sumData", RooArgList(x), dataHisto);
      float totalData=dataHisto->Integral(0,dataHisto->GetXaxis()->GetNbins()+1);

      //data driven template from control region
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

      //closure test in MC
      RooDataHist* sumMC = new RooDataHist("sumMC", "sumMC", RooArgList(x), mcSumHisto[0]);
      RooDataHist* replacementMCTemplateHist = new RooDataHist("replacementmctemplatehist","replacementmctemplatehist",RooArgList(x), dyReplacementMCHisto[0]);
      RooHistPdf replacementMCTemplate("replacementmctemplate","replacementmctemplate",RooArgSet(x),*replacementMCTemplateHist);
      RooDataHist *closureMCTemplate = new RooDataHist("closureMCTemplate", "mcTemplate", RooArgList(x), otherProcsHisto[0]);
      RooHistPdf modelClosureMcTemplate("modelClosureMcTemplate", "modelMcTemplate", RooArgSet(x), *closureMCTemplate);
      RooRealVar ndymcsf("SF_{DY}^{MC}","dyyieldssfactormc",1.0,0.,5.0);
      RooFormulaVar ndymc("N_{DY}^{MC}","@0*@1",RooArgSet(ndyexp,ndymcsf));
      RooAddPdf mcShapeModel("mcshapemodel","signal+background",RooArgSet(replacementMCTemplate,modelClosureMcTemplate),RooArgSet(ndymc,nother));
      RooProdPdf mcModel("mcmodel","(signal+background)*otherconstrain",RooArgSet(mcShapeModel,other_constraint)); 
      mcModel.fitTo(*sumMC,Minos(),Extended(kTRUE),SumW2Error(kTRUE),Constrain(nother),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(fitmin[ich],fitmax[ich]));
      
      report << "***** Fit results for " << templateTitle[ich] << " *********" << endl;
      report << "MC closure SF-DY: "<< ndymcsf.getVal() << " +/- " << ndymcsf.getError() << "\t\t (" << ndymc.getVal() << " vs " << ndyexp.getVal() << " )" << endl;
      for(size_t ivar=0; ivar<runNvars[ich]; ivar++)
	{
	  //template for non-DY processes
	  RooDataHist *mcTemplate = new RooDataHist("mcTemplate", "mcTemplate", RooArgList(x), otherProcsHisto[ivar]);
	  RooHistPdf modelMcTemplate("modelMcTemplate", "modelMcTemplate", RooArgSet(x), *mcTemplate);
	  	  
	  //fit to data
	  RooAddPdf shapeModel("shapemodel","signal+background",RooArgSet(modelDataTemplate,modelMcTemplate),RooArgSet(ndy,nother));
	  RooProdPdf model("model","(signal+background)*otherconstrain",RooArgSet(shapeModel,other_constraint)); 
	  model.fitTo(*sumData,Minos(),Extended(kTRUE),Constrain(nother),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(fitmin[ich],fitmax[ich]));

	  report << "\t" << vars[ivar] << "\t" << ndy.getVal() << " +/- " << ndysf.getError()*ndyexp.getVal() << flush;
	  if(ivar<=0)   report << "\t SF-DY: "<< ndysf.getVal() << " +/- " << ndysf.getError() << "\t\t (" << nother.getVal() << " +/- " << nother.getError() << " )" << flush;
	  report << endl;

	  //
	  // Display results
	  //
	  if(ivar>0) continue;
	  TString cnvTitle[]={"closurec","datac"};
	  for(size_t ifit=0; ifit<2; ifit++)
	    {
	      TCanvas *cnv = getNewCanvas(cnvTitle[ifit],cnvTitle[ifit],false);
	      cnv->SetCanvasSize(600,600);
	      cnv->SetWindowSize(600,600);
	      cnv->cd();

	      TString caption("CMS preliminary");
	      TString fitResults("#splitline");
	      
	      RooPlot *genericFrame = x.frame(Range(fitmin[ich],fitmax[ich]));
	      genericFrame->GetXaxis()->SetTitle( dataHisto->GetXaxis()->GetTitle() );
	      genericFrame->GetYaxis()->SetTitle( "Events" );
	      if(ifit==0)
		{
		  sumMC->plotOn(genericFrame,Name("data"),SumW2Error(kTRUE)); 
		  mcShapeModel.plotOn(genericFrame,Range(fitmin[ich],dataHisto->GetXaxis()->GetXmax()),Name("total"));
		  mcShapeModel.plotOn(genericFrame,
				      RooFit::Components(replacementMCTemplate),
				      Range(fitmin[ich],dataHisto->GetXaxis()->GetXmax()),
				      RooFit::LineStyle(kDashed),
				      Name("dytemplate"));
		  caption="CMS simulation";
		  char buf[100];
		  sprintf(buf,"{SF_{DY}=%3.2f #pm %3.2f}",ndymcsf.getVal(),ndymcsf.getError());
		  fitResults += buf;
		  sprintf(buf,"{N_{DY}=%3.1f #pm %3.1f}",ndymc.getVal(),ndymcsf.getError()*ndyexp.getVal());
		  fitResults += buf;
		}
	      else
		{
		  sumData->plotOn(genericFrame,Name("data"));
		  shapeModel.plotOn(genericFrame,Range(fitmin[ich],dataHisto->GetXaxis()->GetXmax()),Name("total"));
		  shapeModel.plotOn(genericFrame,
				    RooFit::Components(modelDataTemplate),
				    Range(fitmin[ich],dataHisto->GetXaxis()->GetXmax()),
				    RooFit::LineStyle(kDashed),
				    Name("dytemplate"));
		  char buf[100];
		  sprintf(buf,"{SF_{DY}=%3.2f #pm %3.2f}",ndysf.getVal(),ndysf.getError());
		  fitResults += buf;
		  sprintf(buf,"{N_{DY}=%3.1f #pm %3.1f}",ndy.getVal(),ndysf.getError()*ndyexp.getVal());
		  fitResults += buf;
		}
	      genericFrame->Draw();
	      
	      TLegend *leg=new TLegend;
	      formatForCmsPublic(cnv,leg,caption,1);
	      leg->Delete();
	      leg = new TLegend(0.3,0.12,0.5,0.22);
	      leg->SetBorderSize(0);
	      leg->SetFillStyle(0);
	      leg->SetTextFont(42);
	      if(ifit==0) leg->AddEntry("data","MC","p");
	      else        leg->AddEntry("data","data","p");
	      leg->AddEntry("total","Total","l");
	      leg->AddEntry("dytemplate","DY template","l");
	      leg->Draw("same");

	      TPaveText *pave = new TPaveText(0.3,0.22,0.5,0.3,"NDC");
	      pave->SetBorderSize(0);
	      pave->SetFillStyle(0);
	      pave->SetTextFont(42);
	      pave->AddText(fitResults);
	      pave->Draw("same");
	  
	      //display likelihood as inset
	      TPad *npad = new TPad("llpad","ll", 0.65, 0.65, 0.95, 0.92);
	      npad->Draw();
	      npad->cd();
	      RooNLLVar *nll = 0;
	      RooAbsReal* profileLL = 0;
	      if(ifit==0)
		{
		  nll = (RooNLLVar*) mcModel.createNLL(*sumMC,RooFit::CloneData(kFALSE),SumW2Error(kTRUE),Extended(kTRUE),Constrain(nother),Range(fitmin[ich],fitmax[ich]));
		  profileLL = nll->createProfile(ndymcsf);
		}
	      else   
		{
		  nll = (RooNLLVar*) model.createNLL(*sumData,RooFit::CloneData(kFALSE),Extended(kTRUE),Constrain(nother),Range(fitmin[ich],fitmax[ich]));
		  profileLL = nll->createProfile(ndysf);
		}

	      //no need to minimize
	      //RooMinuit minuit(*nll); 
	      //minuit.migrad();
	      //minuit.hesse();
	      
	      RooPlot *frame2 = (ifit==0 ? ndymcsf.frame(Range(0.5,3)) : ndysf.frame(Range(0.5,3)));
	      profileLL->plotOn(frame2,ShiftToZero(),Name("ll"));
	      frame2->GetXaxis()->SetTitle("SF_{DY}");
	      frame2->GetXaxis()->SetTitleOffset(0.8);
	      frame2->GetYaxis()->SetTitle("-log(L/L_{max})");
	      frame2->GetYaxis()->SetTitleOffset(1);
	      frame2->GetYaxis()->SetRangeUser(0,5);
	      frame2->Draw();
	      
	      //update and save
	      cnv->Modified();
	      cnv->Update();
	      cnv->SaveAs(ch[ich] + TString(ifit==0? "mc":"") + "dydistfit.C");
	      cnv->SaveAs(ch[ich] + TString(ifit==0? "mc":"") + "dydistfit.png");
	      cnv->SaveAs(ch[ich] + TString(ifit==0? "mc":"") + "dydistfit.pdf");
	    }
	}
      
      report << "******************************************" << endl;
    }
  
  //all done... print report
  cout << report.str() << endl;
}
