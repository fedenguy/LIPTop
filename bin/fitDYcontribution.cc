#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooConstVar.h"
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
#include "RooKeysPdf.h"

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

void printHelp();
void showDYFitResults(RooRealVar &x, RooDataHist &data, RooAbsPdf &model,RooAbsPdf &dyModel, RooRealVar &ndysf,RooNLLVar &nll,
		      TString tag, TString caption, TString xvar);

//
void printHelp()
{
  printf("--in      --> input file with standard control plots from showControlPlots\n");
  printf("--ttrep   --> input file with control plots from showControlPlots for the syst. json (optional)\n");
  printf("--syst    --> will run systematics also\n");
  printf("--smooth  --> smooth histograms before templating\n");
  printf("command line example: fitDYcontribution --in std_plotter.root --ttrep syst_plotter.root\n");
}

//
int main(int argc,char *argv[])
{

  bool doSyst(false),doSmoothing(false);
  TString stdUrl(""),dyReplacementUrl("");
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos) { printHelp();    return 0; }
      if(arg.find("--in")!=string::npos)                  { stdUrl=argv[i+1];           gSystem->ExpandPathName(stdUrl);            i++;  printf("in      = %s\n", stdUrl.Data()); }
      if(arg.find("--ttrep")!=string::npos)               { dyReplacementUrl=argv[i+1]; gSystem->ExpandPathName(dyReplacementUrl);  i++;  printf("ttRep   = %s\n", dyReplacementUrl.Data()); }
      if(arg.find("--syst")!=string::npos)                { doSyst=true;                                                            printf("Will run systematics\n"); }
      if(arg.find("--smooth")!=string::npos)              { doSmoothing=true;                                                       printf("Will smooth templates\n"); }
    }
  if(stdUrl=="") { printHelp(); return 0; }

  //report to produce
  stringstream report;
  
  //global definitions
  setStyle();
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().getStream(0).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);

  //systematics to consider for the templates
  std::vector<TString> systVars;
  systVars.push_back("");
  if(doSyst)
    {
      systVars.push_back("puup");
      systVars.push_back("pudown");
      systVars.push_back("jesup");
      systVars.push_back("jesdown");
      systVars.push_back("jerup");
      systVars.push_back("jerdown");
    }
  const size_t nsystVars=systVars.size();
  
  TString procs[]={"Di-boson","Single top", "W#rightarrow l#nu","Z#rightarrow ll","t#bar{t}","other t#bar{t}"};

  TString ch[]               ={"ee1jets",                      "ee",                    "mumu1jets",                      "mumu",                    "emu"};
  size_t runNsystVars[]      ={1,                              nsystVars,                   1,                                 nsystVars,                    nsystVars};
  TString signalRegionHisto[]={"ee_eq1jetsdilarccosine",       "ee_dilarccosine",       "mumu_eq1jetsdilarccosine",       "mumu_dilarccosine",       "emu_mtsum"};
  TString templateHisto[]    ={"ee_eq1jetslowmetdilarccosine", "ee_lowmetdilarccosine", "mumu_eq1jetslowmetdilarccosine", "mumu_lowmetdilarccosine", "emu_mtsum"};
  TString templateTitle[]    ={"Z#rightarrow ee (=1 jets)",    "Z#rightarrow ee",       "Z#rightarrow #mu#mu (=1 jets)",  "Z#rightarrow #mu#mu",     "Z#rightarrow #tau#tau"};
  TString templateName[]     ={"dytoee",                       "dytoee",                "dytomumu",                       "dytomumu",                "dytoemu"};
  //  bool rebin[]               ={true,                           true,                    true,                             true,                      false};
  bool rebin[]               ={false,                           false,                   false,                            false,                      false};

  const size_t nchs=sizeof(ch)/sizeof(TString);
  for(size_t ich=0; ich<nchs; ich++)
    {
      //DY->tautau is taken from dedicated file
      if(ch[ich]=="emu" && dyReplacementUrl.IsNull()) continue;

      //open data-driven file
      TFile *dyReplacementFile=TFile::Open( (ch[ich]=="emu") ? dyReplacementUrl : stdUrl);
      
      //get the DY template (data-driven)
      TString dataDir("data");
      if(ch[ich]=="emu") dataDir="Z#rightarrow #tau#tau (data)";
      TH1F *dyReplacementHisto=(TH1F *) dyReplacementFile->Get(dataDir+"/"+templateHisto[ich]);
      dyReplacementHisto->SetDirectory(0); 
      dyReplacementHisto->SetName(templateName[ich]+"model");
      dyReplacementHisto->SetTitle(templateTitle[ich]);
      if(rebin[ich]) dyReplacementHisto->Rebin();
      dyReplacementFile->Close();

      //open standard file
      TFile *stdFile=TFile::Open(stdUrl);
      
      //get histogram for data
      TH1F *dataHisto=(TH1F *) stdFile->Get("data/"+signalRegionHisto[ich])->Clone("data");
      dataHisto->SetDirectory(0);
      dataHisto->SetTitle("data");
      if(rebin[ich]) dataHisto->Rebin();
      
      //get histograms for MC
      std::vector<TH1F *> dyReplacementMCHisto, dyMCHisto, otherProcsHisto, mcSumHisto;
      const size_t nprocs=sizeof(procs)/sizeof(TString);
      for(size_t ivar=0; ivar<runNsystVars[ich]; ivar++)
	{
	  for(size_t iproc=0; iproc<nprocs; iproc++)
	    {
	      TH1F *histo = (TH1F *) stdFile->Get(procs[iproc]+"/"+signalRegionHisto[ich]+systVars[ivar]);
	      if(rebin[ich]) histo->Rebin();

	      //save histogram in the appropriate category
	      if(procs[iproc].Contains("Z#rightarrow"))   
		{
		  histo = (TH1F *)histo->Clone(templateName[ich]+systVars[ivar]+"mc");
		  histo->SetTitle(templateTitle[ich]+" (MC)");
		  histo->SetDirectory(0);		  
		  dyMCHisto.push_back( histo );
		  
		  TH1F *repHistoMC = (TH1F *) stdFile->Get(procs[iproc]+"/"+templateHisto[ich]+systVars[ivar]);
		  repHistoMC       = (TH1F *) repHistoMC->Clone(templateName[ich]+systVars[ivar]+"repmc");
		  if(rebin[ich]) repHistoMC->Rebin();
		  repHistoMC->SetDirectory(0);
		  dyReplacementMCHisto.push_back( repHistoMC );
		}
	      else if(otherProcsHisto.size()==ivar)
		{
		  histo = (TH1F *) histo->Clone("otherprocs"+systVars[ivar]);
		  histo->SetDirectory(0);
		  histo->SetTitle("Other processes");
 		  otherProcsHisto.push_back( histo );
		}
	      else otherProcsHisto[ivar]->Add(histo);
	      
	      //total mc
	      if(iproc==0) 
		{
		  histo = (TH1F *)histo->Clone("mcsum"+templateName[ich]+systVars[ivar]);
		  histo->SetDirectory(0);
		  histo->SetTitle("Total MC");
		  mcSumHisto.push_back( histo );
		}
	      else mcSumHisto[ivar]->Add(histo);
	    }
	}
      stdFile->Close();

      //
      // Now fit
      //
      float totalData=dataHisto->Integral(1,dataHisto->GetXaxis()->GetNbins());
      float dyExpected=dyMCHisto[0]->Integral(1,dyMCHisto[0]->GetXaxis()->GetNbins());
      float totalOthers=otherProcsHisto[0]->Integral(1,otherProcsHisto[0]->GetXaxis()->GetNbins());
      float uncOthers=max(float(0.04*totalOthers),float(sqrt(totalOthers)));
            
      //define variable
      RooRealVar x("x","x", dataHisto->GetXaxis()->GetXmin(), dataHisto->GetXaxis()->GetXmax());
      
      //data
      RooDataHist* sumData = new RooDataHist("sumData", "sumData", RooArgList(x), dataHisto);
      
      //data (or MC) driven template for DY
      if(doSmoothing) dyReplacementHisto->Smooth();
      RooDataHist* dataTemplate = new RooDataHist("dataTemplate", "dataTemplate", RooArgList(x), dyReplacementHisto );
      RooHistPdf modelDataTemplate("modelDataTemplate", "modelDataTemplate", RooArgSet(x), *dataTemplate);
      RooRealVar ndysf("SF_{DY}","dyyieldssfactor",1.0,0.5,3.0);
      RooFormulaVar ndy("N_{DY}","@0*@1",RooArgSet(RooConst(dyExpected),ndysf));

      //MC based template for other processes
      if(doSmoothing) otherProcsHisto[0]->Smooth();
      RooDataHist *otherMCTemplate = new RooDataHist("otherMCTemplate", "otherMCTemplate", RooArgList(x), otherProcsHisto[0]);
      RooHistPdf modelOtherMCTemplate("modelOtherMCTemplate", "modelOtherMCTemplate", RooArgSet(x), *otherMCTemplate);
      RooRealVar nother("N_{other}","otheryields",totalOthers,0,totalData);
      RooGaussian other_constraint("otherconstraintpdf","otherconstraintpdf", nother, RooConst(totalOthers), RooConst(uncOthers));
      
      //the model of the data (DY+other)xconstraints
      RooAddPdf shapeModelData("shapeModelData","signal+background", RooArgSet(modelDataTemplate,modelOtherMCTemplate),RooArgSet(ndy,nother));
      RooProdPdf constrShapeModelData("constrShapeModelData","shape*constrain",RooArgSet(shapeModelData,other_constraint));
      
      report << endl
	     << "FIT REPORT FOR " << signalRegionHisto[ich] << endl
	     << "---------------------------------------------------" << endl
	     << "[Observed] " << totalData << endl 
	     << "---------------------------------------------------" << endl
	     << "[DY expected] " <<  dyExpected << endl 
	     << "[Others expected]:" << totalOthers << " +/- " << uncOthers << endl
	     << "---------------------------------------------------" << endl;

      //fit data
      constrShapeModelData.fitTo(*sumData,Extended(kTRUE),Constrain(nother),Save(kTRUE));
      report << "[SF-DY] "      << ndysf.getVal() << " +/- " << ndysf.getError()  << endl
	     << "[Others fit] " << nother.getVal() << " +/- " << nother.getError() << endl
	     << "---------------------------------------------------" << endl;
      RooNLLVar *nll = (RooNLLVar*) constrShapeModelData.createNLL(*sumData,CloneData(kFALSE),Extended(kTRUE),Constrain(nother));
      showDYFitResults(x,*sumData,constrShapeModelData,modelDataTemplate,ndysf,*nll,
		       signalRegionHisto[ich],"CMS preliminary",dyReplacementHisto->GetXaxis()->GetTitle());
      
      //do the closure tests on MC for different variations 
      report << "[#Delta SF-DY in closure test]" << endl; 
      float nominalMCSF(0);
      for(size_t ivar=0; ivar<runNsystVars[ich]; ivar++)
	{
	  RooDataHist* sumMC   = new RooDataHist("sumMC",   "sumMC",   RooArgList(x), mcSumHisto[0]);

	  //float totalMC=mcSumHisto[0]->Integral();
	  float ndyExp=dyMCHisto[ivar]->Integral(); 
	  float totalOthers=otherProcsHisto[ivar]->Integral();
	  float uncOthers=max(float(sqrt(totalOthers)),float(0.04*totalOthers)); 

	  if(doSmoothing ) dyReplacementMCHisto[ivar]->Smooth();
	  RooDataHist* mcTemplate = new RooDataHist("mcTemplate", "mcTemplate", RooArgList(x), dyReplacementMCHisto[ivar] );
	  RooHistPdf modelMCTemplate("modelMCTemplate", "modelMCTemplate", RooArgSet(x), *mcTemplate);
	  RooRealVar ndysf("SF_{DY}","dyyieldssfactor",1.0,0.5,3.0);
	  RooFormulaVar ndy("N_{DY}","@0*@1",RooArgSet(RooConst(ndyExp),ndysf));

	  if(doSmoothing && ivar>0)  otherProcsHisto[ivar]->Smooth();
	  RooDataHist *otherMCTemplate = new RooDataHist("otherMCTemplate", "otherMCTemplate", RooArgList(x), otherProcsHisto[ivar]);
	  RooHistPdf modelOtherMCTemplate("modelOtherMCTemplate", "modelOtherMCTemplate", RooArgSet(x), *otherMCTemplate);
	  RooRealVar nother("N_{other}","otheryields",totalOthers,0,totalData);
	  RooGaussian other_constraint("otherconstraintpdf","otherconstraintpdf", nother, RooConst(totalOthers), RooConst(uncOthers));
       	  
	  RooAddPdf shapeModelMC("shapeModelMC","signal+background",RooArgSet(modelMCTemplate,modelOtherMCTemplate),RooArgSet(ndy,nother));
	  RooProdPdf constrShapeModelMC("constrShapeModelMC","shape*constrain",RooArgSet(shapeModelMC,other_constraint)); 

	  constrShapeModelMC.fitTo(*sumMC,Extended(kTRUE),SumW2Error(kTRUE),Constrain(nother),Save(kTRUE));
	  report << (ivar!=0 ? systVars[ivar] + " syst. variation: " : " Nominal: ") << " " 
		 << ndy.getVal()/dyMCHisto[0]->Integral()-nominalMCSF << endl;
	  if(ivar==0) nominalMCSF=ndysf.getVal();
	}      
    }
  
  //all done
  cout << report.str() << endl;
}


//
void showDYFitResults(RooRealVar &x,
		      RooDataHist &data,
		      RooAbsPdf &model,
		      RooAbsPdf &dyModel,
		      RooRealVar &ndysf,
		      RooNLLVar &nll,
		      TString tag, 
		      TString caption, 
		      TString xvar)
{
  TCanvas *cnv = new TCanvas(tag+"c",tag+"c",600,600);
  cnv->cd();
  RooPlot *frame = x.frame();
  if(tag.Contains("mc")) data.plotOn(frame,Name("data"),DataError(RooAbsData::SumW2));
  else                   data.plotOn(frame,Name("data"));
  frame->GetXaxis()->SetTitle( xvar );     frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitle( "Events" ); frame->GetYaxis()->SetTitleOffset(1.2);

  model.plotOn(frame,Name("total"),MoveToBack());
  model.plotOn(frame,Components(dyModel),FillStyle(1001),FillColor(kGray),LineColor(kGray+1),DrawOption("lf"),Name("dytemplate"),MoveToBack());



  frame->Draw();
  Double_t ymax=frame->GetMaximum();
  frame->GetYaxis()->SetRangeUser(0.01,max(float(ymax*1.4),float(1.0)));

  TLegend *leg=new TLegend(0.2,0.82,0.55,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetNColumns(3);
  leg->AddEntry("data","data","p");
  leg->AddEntry("total","Total","l");
  leg->AddEntry("dytemplate","DY template","f");
  leg->Draw("same");

  TPaveText *pave = new TPaveText(0.2,0.85,0.5,0.92,"NDC");
  pave->SetBorderSize(0);
  pave->SetTextAlign(12);
  pave->SetFillStyle(0);
  pave->SetTextFont(42);
  pave->AddText(caption);
  char buf[200]; sprintf(buf,"SF_{DY}=%3.2f #pm %3.2f",ndysf.getVal(),ndysf.getError());
  pave->AddText(buf);
  pave->Draw("same");
	  
  //display likelihood as inset
  TPad *npad = new TPad("llpad","ll", 0.6, 0.65, 0.88, 0.93);
  npad->Draw();
  npad->cd();
  frame = ndysf.frame();
  nll.plotOn(frame,ShiftToZero(),Name("ll"));
  frame->GetXaxis()->SetTitle("SF_{DY}");  
  frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetXaxis()->SetTitleSize(0.1);
  frame->GetXaxis()->SetLabelSize(0.08);
  frame->GetYaxis()->SetTitle("-log(L/L_{max})");
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitleSize(0.1);
  frame->GetYaxis()->SetLabelSize(0.08);
  frame->GetYaxis()->SetRangeUser(0,5);
  frame->Draw();
	      
  //update and save
  cnv->Modified();
  cnv->Update();
  cnv->SaveAs(tag+".png");
  cnv->SaveAs(tag+".C");
  cnv->SaveAs(tag+".pdf");
}



