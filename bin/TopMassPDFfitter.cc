#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooAddPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooSimPdfBuilder.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooNovosibirsk.h"
#include "RooCBShape.h"
#include "RooGamma.h"
#include "RooLognormal.h"
#include "RooExponential.h"
#include "RooProdPdf.h"

#include "TList.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSystem.h"
#include "TProfile.h"

#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include<vector>

typedef std::pair<float,float> Value_t;
typedef std::map<std::string,Value_t> FitResults_t;

using namespace std; 
using namespace RooFit ;

enum Channel{INCLUSIVE,SAMEFLAVOR,OPFLAVOR};
FitResults_t SignalPDFs(TString url="EventSummaries.root",int nbtags=-1,int channel=INCLUSIVE);
FitResults_t BckgPDFs(TString url="MassPDFs.root",int channel=INCLUSIVE);
FitResults_t DYBckgPDFs(TString url,int channel=INCLUSIVE);

//
void printHelp()
{
  printf("--help    --> print this\n");
  printf("--in      --> input file with the summary trees\n");
  printf("--sig     --> run signal fitter\n");
  printf("--bkg     --> run bckg fitter\n");
}


//
int main(int argc, char* argv[])
{
  //parse command line
  TString url("");
  bool fitSignal(false), fitBckg(false);
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos) { printHelp();    return 0; }
      if(arg.find("--in")!=string::npos && i+1<argc)  { url=argv[i+1];         gSystem->ExpandPathName(url);        i++;  printf("in      = %s\n", url.Data()); }
      if(arg.find("--sig")!=string::npos)             { fitSignal=true;                                                   printf("signal PDF will be fit\n"); }
      if(arg.find("--bkg")!=string::npos)             { fitBckg=true;                                                     printf("bckg   PDF will be fit\n"); }
    }
  if(url=="" || (!fitSignal && !fitBckg) ) { printHelp(); return 0; }

  //prepare to report
  stringstream report;

  //keep RooFit quiet                                                                                                                                                                              
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().getStream(0).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(0).removeTopic(Fitting);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);


  // load framework libraries                                                                                                                                                            
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
    
  //fit the templates
  int ch[]={INCLUSIVE,SAMEFLAVOR,OPFLAVOR};
  TString chName[]={"inclusive","same flavor","op. flavor"};

  //signal fit
  if(fitSignal)
    {
      for(size_t i=0; i<3;i++)
	{
	  std::map<string,FitResults_t> sigFitResults;
	  sigFitResults["[CATEGORY #1]"]       = SignalPDFs(url,1,ch[i]);
	  sigFitResults["[CATEGORY #2]"]       = SignalPDFs(url,2,ch[i]);
	  
	  //display the results
	  report << " *************** TopMassPDFfitter signal results for " << chName[i] << " channel *********************** " << endl;
	  int catCtr(0);
	  for(std::map<string,FitResults_t>::iterator it = sigFitResults.begin(); it != sigFitResults.end(); it++,catCtr++)
	    {
	      report << it->first << endl;
	      TString catPostFix("_s"); catPostFix += catCtr;
	      for(FitResults_t::iterator itt = it->second.begin(); itt != it->second.end(); itt++)
		report << itt->first << catPostFix << ":" << itt->second.first  << "\t +/-" << itt->second.second << endl;
	    }
	}
    }

  //background fits
  if(fitBckg)
    {
      for(size_t i=0; i<3;i++)
	{
	  std::map<string,FitResults_t> bckgFitResults;
	  bckgFitResults["[NON DY]"] = BckgPDFs(url,ch[i]);      
	  bckgFitResults["[DY    ]"] = DYBckgPDFs(url,ch[i]);      
	  report << " *************** TopMassPDFfitter bckg results for " << chName[i] << " channel *********************** " << endl;
	  for(std::map<string,FitResults_t>::iterator it = bckgFitResults.begin(); it != bckgFitResults.end(); it++)
	    { 
	      for(FitResults_t::iterator itt = it->second.begin(); itt!= it->second.end(); itt++)
		report << itt->first << ":" << itt->second.first << "\t +/-" << itt->second.second << endl;
	    }
	  report << endl;
	}
    }

  //all done
  cout << report.str() << endl;
  return -1;
}


//
FitResults_t BckgPDFs(TString url,int channel)
{
  //get the histograms and build the sum
  std::map<string, TH1D *> bckgDists;
  std::map<string,std::vector<string> > bckgProcs;

  std::vector<string> topSamples;
  topSamples.push_back("SingleTbar_tW");
  topSamples.push_back("SingleT_tW");
  topSamples.push_back("SingleTbar_t");
  topSamples.push_back("SingleT_t");
  topSamples.push_back("SingleTbar_s");
  topSamples.push_back("SingleT_s");
  topSamples.push_back("TTJets");
  bckgProcs["Top"]=topSamples;
  
  std::vector<string> ewkSamples;
  ewkSamples.push_back("WW");
  ewkSamples.push_back("ZZ");
  ewkSamples.push_back("WZ");
  //  ewkSamples.push_back("WJetsToLNu");
  bckgProcs["EWK"]=ewkSamples;

  //build the mass histograms and the weighted dataset for the non DY background
  int ipt(1);
  TFile *f = TFile::Open(url);

  RooRealVar mass("m","Mass",100,300,"GeV/c^{2}");
  mass.setBins(40);
  RooRealVar eventWeight("weight","weight",1,0,100);
  RooDataSet dh("dh","dh",RooArgSet(mass,eventWeight),"weight");       
  for(std::map<string,std::vector<string> >::iterator it = bckgProcs.begin(); it != bckgProcs.end(); it++,ipt++)
    {
      TString sName("m"); sName += (ipt+1);      
      TH1D *h= new TH1D(sName,sName,80,100,500);      
      h->SetDirectory(0);
      h->GetYaxis()->SetTitle("Events / (10 GeV/c^{2})");
      h->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
      h->GetXaxis()->SetTitleOffset(0.8);
      h->GetYaxis()->SetTitleOffset(0.8);
      h->SetTitle(it->first.c_str());
      bckgDists[it->first]=h;

      for(std::vector<string>::iterator pit = it->second.begin(); pit != it->second.end(); pit++)
	{
      	  TString tname=*pit + "/data";
	  TTree *t = (TTree *) f->Get(tname);
	  if(t==0) continue;

	  Int_t cat;
	  Float_t weight,xsecWeight;
	  Float_t evmeasurements[10];
	  t->GetBranch("cat")->SetAddress(&cat);
	  t->GetBranch("weight")->SetAddress(&weight);
	  t->GetBranch("xsecWeight")->SetAddress(&xsecWeight);
	  t->GetBranch("evmeasurements")->SetAddress(evmeasurements);
	  for(Int_t i=0; i<t->GetEntriesFast(); i++)
	    {
	      t->GetEntry(i);
	      if(channel==SAMEFLAVOR && cat!=EE && cat!=MUMU) continue;
	      if(channel==OPFLAVOR && cat!= EMU) continue;
	      if(evmeasurements[0]>0) 
		{
		  h->Fill(evmeasurements[0],weight*xsecWeight);
		  mass=evmeasurements[0];
		  eventWeight=weight*xsecWeight;
		  dh.add(RooArgSet(mass,eventWeight),weight*xsecWeight,0.);
		}
	    }
	}
    }    
  f->Close();
  
  //the sum
  TH1D *sumH= (TH1D* )bckgDists["Top"]->Clone("totalbckg");
  sumH->SetDirectory(0);
  sumH->Add(bckgDists["EWK"]);
  
  //prepare the pdfs for the fit in mass
  RooRealVar mpv_l("mpv_{l}","Mpv of landau",140,100,165);
  RooRealVar sigma_l("#sigma_{l}","Sigma of landau",20,10,25);
  RooLandau lan("blandau","Mass component 1",mass,mpv_l,sigma_l);
  RooRealVar mean_g("#mu_{g}","Mean of gaus",180,160,200);
  RooRealVar sigma_g("#sigma_{g}","Sigma of gaus",20,10,25);
  RooGaussian gaus("bgauss","Mass component 2",mass,mean_g,sigma_g);   
  RooRealVar massfrac("#alpha","Fraction of component 1",0.1,0.0,0.2);
  RooLandau massmodel(lan);
  //  RooAddPdf massmodel("model","Model",RooArgList(gaus,lan),massfrac);
  //  RooNumConvPdf massmodel("model","Model",mass,lan,gaus);
  massmodel.fitTo(dh,Save(kTRUE),Range(100.,300.),SumW2Error(kTRUE));
  RooDataHist binnedData("binnedData","binned version of dh",RooArgSet(mass),dh);
  RooChi2Var chi2("chi2","chi2",massmodel,binnedData,DataError(RooAbsData::SumW2)) ;
  RooMinuit m(chi2) ;
  m.migrad() ;
  m.hesse() ;

  //show the results of the fit
  setStyle();
  TString chName("");
  if(channel==OPFLAVOR) chName += "_of";
  if(channel==SAMEFLAVOR) chName += "_sf";
  TString cnvName("NonDYPDF"); cnvName += chName; 
  TCanvas *c = new TCanvas(cnvName,cnvName,1800,600);
  c->SetWindowSize(1600,800);
  c->Divide(2,1);

  TPad *p=(TPad *)c->cd(1);
  RooPlot* frame = mass.frame(Title("Combined"));
  dh.plotOn(frame,DrawOption("pz"),DataError(RooAbsData::SumW2)); 
  massmodel.plotOn(frame,Components(lan),LineColor(kGray));
  massmodel.plotOn(frame);
  frame->Draw();
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitle("Events (A.U.)");
  frame->GetXaxis()->SetTitleOffset(0.8);
  frame->GetXaxis()->SetTitle("Reconstructed Mass [GeV/c^{2}]");
  TLegend *leg = p->BuildLegend();
  formatForCmsPublic(p,leg,"CMS simulation",1);
  leg->Delete();	  
  
  p=(TPad *)c->cd(2);
  p->Divide(1,2);
  int ibckg(1);
  for(std::map<std::string, TH1D *>::iterator it = bckgDists.begin();
      it != bckgDists.end();
      it++,ibckg++)
    {
      p->cd(ibckg);
      it->second->Rebin();
      it->second->Draw("hist");  
      it->second->GetXaxis()->SetTitleOffset(1.0);
      it->second->GetYaxis()->SetTitleOffset(1.0);
      
      TPaveText *pave = new TPaveText(0.6,0.7,0.8,0.8,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      pave->SetTextFont(42);
      pave->AddText(it->first.c_str());
      pave->Draw("same");
    }

  c->Modified();
  c->Update();
  c->SaveAs(cnvName+".C");
  c->SaveAs(cnvName+".pdf");
  c->SaveAs(cnvName+".png");
  
  //return the result
  FitResults_t fitPars;
  fitPars["bckgmpv_{l}"]    = Value_t(mpv_l.getVal(),mpv_l.getError());
  fitPars["bckg#sigma_{l}"] = Value_t(sigma_l.getVal(),sigma_l.getError());
  //  fitPars["bckg#mu_{g}"]    = Value_t(mean_g.getVal(),mean_g.getError());
  // fitPars["bckg#sigma_{g}"] = Value_t(sigma_g.getVal(),sigma_g.getError());
  // fitPars["bckg#alpha"]     = Value_t(massfrac.getVal(),massfrac.getError()); 
  return fitPars;
}

//
FitResults_t DYBckgPDFs(TString url,int channel)
{
  FitResults_t fitPars; 
  if(channel==INCLUSIVE) return fitPars;

  //build the datasets for DY backgrounds
  std::map<string,std::vector<string> > bckgProcs;

  std::vector<string>dyLLSamples;
 
  //MC
  dyLLSamples.push_back("DYJetsToLL");
  //   dyLLSamples.push_back("DYJetsToMuMu_M20to50");
  //   dyLLSamples.push_back("DYJetsToEE_M20to50");
  bckgProcs["DYLLMC"]=dyLLSamples;

  //DATA-DRIVEN
  if(channel!=OPFLAVOR)
    {
      dyLLSamples.clear();
      dyLLSamples.push_back("DoubleMuMay10ReReco");
      dyLLSamples.push_back("DoubleElectronMay10ReReco");
      dyLLSamples.push_back("MuEGMay10ReReco");
      dyLLSamples.push_back("DoubleMuPromptRecov4");
      dyLLSamples.push_back("DoubleElectronPromptRecov4");
      dyLLSamples.push_back("MuEGPromptRecov4");
      dyLLSamples.push_back("DoubleMu05AugReReco");
      dyLLSamples.push_back("DoubleElectron05AugReReco");
      dyLLSamples.push_back("MuEG05AugReReco");
      dyLLSamples.push_back("DoubleMuPromptRecov6");
      dyLLSamples.push_back("DoubleElectronPromptRecov6");
      dyLLSamples.push_back("MuEGPromptRecov6");
      bckgProcs["DYLLData"]=dyLLSamples;
    }
  if(channel==OPFLAVOR)
    {
      dyLLSamples.clear();
      dyLLSamples.push_back("DoubleMu2011AAug5thTauReplacement");
      dyLLSamples.push_back("DoubleMu2011AMay10thTauReplacement");
      dyLLSamples.push_back("DoubleMu2011APRv4TauReplacement");
      dyLLSamples.push_back("DoubleMu2011APRv6TauReplacement");
      bckgProcs["DYLLData"]=dyLLSamples;
    }
  //build the mass datasets
  RooRealVar mass("m","Mass",100,300,"GeV/c^{2}");
  RooRealVar eventWeight("weight","weight",1,0,100);
  RooDataSet mch("mch","mch",RooArgSet(mass,eventWeight),"weight");
  RooDataSet ctrlmch("ctrlmch","ctrlmch",RooArgSet(mass,eventWeight),"weight");
  RooDataSet dh("dh","dh",RooArgSet(mass,eventWeight),"weight");
  
  int ipt(1);
  TFile *f = TFile::Open(url);
  for(std::map<string,std::vector<string> >::iterator it = bckgProcs.begin(); it != bckgProcs.end(); it++,ipt++)
    {
      bool isData(true);
      if(it->first.find("MC")!=std::string::npos) isData=false;
      for(std::vector<string>::iterator pit = it->second.begin(); pit != it->second.end(); pit++)
	{
      	  TString tname=*pit + "/data";
	  TTree *t = (TTree *) f->Get(tname);
	  if(t==0) continue;
	  if(t->GetEntriesFast()==0) continue;
	  
	  Int_t cat;
	  Float_t weight,xsecWeight;
	  Float_t evmeasurements[10];
	  t->GetBranch("cat")->SetAddress(&cat);
	  t->GetBranch("weight")->SetAddress(&weight);
	  t->GetBranch("xsecWeight")->SetAddress(&xsecWeight);
	  t->GetBranch("evmeasurements")->SetAddress(evmeasurements);
	  for(Int_t i=0; i<t->GetEntriesFast(); i++)
	    {
	      t->GetEntry(i);
	      bool isZcand(false);
	      isZcand=bool(evmeasurements[2]); 
	      int nbtags=evmeasurements[5];
	      if(channel==SAMEFLAVOR)
		{
		  if(cat!=EE && cat!=MUMU) continue;
		  
		  //for data look only at Z like events
		  if(isData && !isZcand) continue;

		  //for Z like events we only care about 0 b-tags
		  if(isZcand && nbtags>0) continue;
		}
	      if(channel==OPFLAVOR)
		{
		  if(cat!= EMU) continue;
		}
	      
	      if(evmeasurements[0]>0) 
		{
		  mass=evmeasurements[0];
                  if(isData)
		    {
		      eventWeight=1;
		      dh.add(RooArgSet(mass,eventWeight),1.0,0.);
		    }
		  else
		    {
		      eventWeight=weight*xsecWeight;
		      if(!isZcand) mch.add(RooArgSet(mass,eventWeight),weight*xsecWeight,0.);
		      else         ctrlmch.add(RooArgSet(mass,eventWeight),weight*xsecWeight,0.);
		    }
		}
	    }
	}
    }
  f->Close();

  mass.setBins(20);

  //fit the MC template
  RooRealVar mpv_l("mpv_{l}^{MC}","Mpv of landau",140,100,165);
  RooRealVar sigma_l("#sigma_{l}^{MC}","Sigma of landau",20,15,25);
  RooLandau mclan("mclan","Landau component",mass,mpv_l,sigma_l);

  RooRealVar mean_g("mpv_{g}^{MC}","Mean of gauss",180,160,200);
  RooRealVar sigma_g("#sigma_{g}^{MC}","Sigma of gauss",20,15,25);
  RooGaussian mcgauss("mcgaus","Gaussian component",mass,mean_g,sigma_g);

  RooRealVar frac("frac","Model fraction",0.6,0,1.);
  //RooAddPdf mcmodel("mcmodel","Mass model",RooArgList(mcgauss,mclan),frac);
  RooLandau mcmodel(mclan);
 
  if(channel==SAMEFLAVOR)
    {
      RooDataHist binnedMC("binnedMC","binned version of mch",RooArgSet(mass),ctrlmch);
      RooChi2Var chi2("chi2","chi2",mcmodel,binnedMC,DataError(RooAbsData::SumW2)) ;
      RooMinuit m(chi2) ;
      m.migrad() ;
      m.hesse() ;
      fitPars["dybckgmpv_{l}^{ctrl MC}"]    = Value_t(mpv_l.getVal(),mpv_l.getError());
      fitPars["dybckg#sigma_{l}^{ctrl MC}"] = Value_t(sigma_l.getVal(),sigma_l.getError());
      //    fitPars["dybckmean_{g}^{ctrl MC}"]    = Value_t(mean_g.getVal(),mean_g.getError());
      // fitPars["dybckg#sigma_{g}^{ctrl MC}"] = Value_t(sigma_g.getVal(),sigma_g.getError()); 
      // fitPars["dybckgfrac^{ctrl MC}"]       = Value_t(frac.getVal(),frac.getError());
    }
  
  RooDataHist binnedMC("binnedMC","binned version of mch",RooArgSet(mass),mch);
  RooChi2Var chi2("chi2","chi2",mcmodel,binnedMC,DataError(RooAbsData::SumW2)) ;
  RooMinuit m(chi2) ;
  m.migrad() ;
  m.hesse() ;
  fitPars["dybckgmpv_{l}^{MC}"]    = Value_t(mpv_l.getVal(),mpv_l.getError());
  fitPars["dybckg#sigma_{l}^{MC}"] = Value_t(sigma_l.getVal(),sigma_l.getError());
  //  fitPars["dybckmean_{g}^{MC}"]    = Value_t(mean_g.getVal(),mean_g.getError());
  //  fitPars["dybckg#sigma_{g}^{MC}"] = Value_t(sigma_g.getVal(),sigma_g.getError()); 
  // fitPars["dybckgfrac^{MC}"]       = Value_t(frac.getVal(),frac.getError());

  //fit the data-driven template
  RooRealVar dmpv_l("mpv_{l}","Mpv of landau",140,100,165);
  RooRealVar dsigma_l("#sigma_{l}","Sigma of landau",20,10,25);
  RooLandau dlan("dlan","Landau component",mass,dmpv_l,dsigma_l);  
//   RooRealVar dtail_l("tail_{l}","Tail",1.0,-20.,20.);
//   RooExponential dexp("dexp","Mass model",mass,dtail_l);
//   RooRealVar gamma("gamma","gamma",170,100,300);
//   RooRealVar beta("beta","beta",1,0,10);
//   RooRealVar mu("mu","mu",0,0,1);
//   RooGamma dgam("dgam","Gamma component",mass,gamma,beta,mu); 

  RooRealVar dmean_g("mpv_{g}","Mean of gauss",180,160,200);
  RooRealVar dsigma_g("#sigma_{g}","Sigma of gauss",20,10,25);
  RooGaussian dgauss("dgauss","Gauss component",mass,dmean_g,dsigma_g);
    
  //RooLandau datamodel("datamodel","Mass model",mass,dmpv_l,dsigma_l);
  //RooFFTConvPdf datamodel("datamodel","Mass model",mass,dlan,dgauss);

  RooRealVar dfrac("frac","Model fraction",0.1,0.0,0.2);
  //  RooAddPdf datamodel("datamodel","Mass model",RooArgList(dgauss,dlan),dfrac);
  RooLandau datamodel(dlan);
   
  datamodel.fitTo(dh,Save(kTRUE),Range(100.,300.));
  //   RooDataHist binnedData("binnedData","binned version of dh",RooArgSet(mass),dh);
  //   RooChi2Var chi2data("chi2data","chi2data",datamodel,binnedData);
  //   RooMinuit mindata(chi2data) ;
  //   mindata.migrad() ;
  //   mindata.hesse() ;
  
  fitPars["dybckgmpv_{l}"]         = Value_t(dmpv_l.getVal(),dmpv_l.getError());
  fitPars["dybckg#sigma_{l}"]      = Value_t(dsigma_l.getVal(),dsigma_l.getError());
  //fitPars["dybckmean_{g}"]         = Value_t(dmean_g.getVal(),dmean_g.getError());
  //fitPars["dybckg#sigma_{g}"]      = Value_t(dsigma_g.getVal(),dsigma_g.getError());
  // fitPars["dybckgfrac"]            = Value_t(dfrac.getVal(),dfrac.getError());

  //show the results of the fit
  setStyle();
  TString chName("");
  if(channel==OPFLAVOR) chName += "_of";
  if(channel==SAMEFLAVOR) chName += "_sf";
  TString cnvName("DYPDF"); cnvName += chName; 
  TCanvas *c = new TCanvas(cnvName,cnvName,800,800);
  c->SetWindowSize(800,800);
  c->Divide(1,2);
  
  TPad *p=(TPad *) c->cd(1);
  RooPlot* frame = mass.frame(Title("MC Result"),Bins(20));
  mch.plotOn(frame,DataError(RooAbsData::SumW2));
  mcmodel.plotOn(frame,Components(mclan),RooFit::LineColor(kGray));
  mcmodel.plotOn(frame);
  frame->Draw();
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitle("Events (A.U.)");
  frame->GetXaxis()->SetTitleOffset(0.8);
  frame->GetXaxis()->SetTitle("Reconstructed Mass [GeV/c^{2}]");
  TLegend *leg = p->BuildLegend();
  formatForCmsPublic(p,leg,"CMS simulation",1);
  leg->Delete();	  
  
  p=(TPad *) c->cd(2);
  RooPlot* frame2 = mass.frame(Title("Data Result"),Bins(20));
  dh.plotOn(frame2);
  datamodel.plotOn(frame2,Components(dlan),RooFit::LineColor(kGray));
  datamodel.plotOn(frame2);
  frame2->Draw();
  frame2->GetYaxis()->SetTitleOffset(1.2);
  frame2->GetYaxis()->SetTitle("Events");
  frame2->GetXaxis()->SetTitleOffset(0.8);
  frame2->GetXaxis()->SetTitle("Reconstructed Mass [GeV/c^{2}]");
  leg = p->BuildLegend();
  formatForCmsPublic(p,leg,"CMS preliminary",1);
  leg->Delete();	  

  c->Modified();
  c->Update();
  c->SaveAs(cnvName+".C");
  c->SaveAs(cnvName+".pdf");
  c->SaveAs(cnvName+".png");
  
  //return the result
  return fitPars;
}



//
FitResults_t SignalPDFs(TString url,int nbtags,int channel)
{
  //the mass points
  typedef std::pair<TString,Float_t> MassPoint_t;
  std::vector<MassPoint_t> MassPointCollection;
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_161v5",161.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_163v5",163.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_166v5",166.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_169v5",169.5) ); 
  MassPointCollection.push_back( MassPoint_t("TTJets_signal",    172.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_175v5",175.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_178v5",178.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass181v5", 181.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass184v5", 184.5) );

  //define the fit range and the main variable
  float rangeMin(100), rangeMax(500);
  RooRealVar mass("m","Reconstructed Mass [GeV/c^{2}]", rangeMin,rangeMax);
  RooRealVar evWeight("w","Event weight",0,20);

  //define the dataset categorized per top quark mass point
  RooCategory massCategory("cat","cat") ;
  for(size_t ipt=0; ipt<MassPointCollection.size(); ipt++)
    {
      TString cName("m"); cName += (ipt+1);
      massCategory.defineType(cName.Data());
    }
  RooDataSet combData("combData","combData",RooArgSet(mass,massCategory));
  
  //fill the dataset
  TFile *f = TFile::Open(url);
  TGraphErrors *massProfile = new TGraphErrors;
  map<string,TH1*> hmap;
  for(size_t ipt=0; ipt<MassPointCollection.size(); ipt++)
    {
      TString tname=MassPointCollection[ipt].first + "/data";
      TTree *t = (TTree *) f->Get(tname);
      if(t==0) continue;
      
      TString cName("m"); cName += (ipt+1);      
      massCategory.setLabel(cName.Data());
      TH1F *h = new TH1F(cName,cName,20,rangeMin,rangeMax);
      h -> SetDirectory(0);
      hmap[cName.Data()] = h;
     
      Int_t cat;
      Float_t normWeight,weight;
      Float_t evmeasurements[10];
      t->GetBranch("cat")->SetAddress(&cat);
      t->GetBranch("normWeight")->SetAddress(&normWeight);
      t->GetBranch("weight")->SetAddress(&weight);
      t->GetBranch("evmeasurements")->SetAddress(evmeasurements);
      for(Int_t i=0; i<t->GetEntriesFast(); i++)
	{
	  t->GetEntry(i);
	  if(channel==SAMEFLAVOR && cat!=EE && cat!=MUMU) continue;
	  if(channel==OPFLAVOR && cat!= EMU) continue;
	  if(nbtags==0 && evmeasurements[5]!=0) continue;
	  if(nbtags==1 && evmeasurements[5]!=1) continue;
	  if(nbtags==2 && evmeasurements[5]<1) continue;
	  if(normWeight==0) continue;
	  mass=evmeasurements[0];
	  if(mass.getVal()<rangeMin || mass.getVal() > rangeMax) continue;
	  evWeight=weight;
	  combData.add(RooArgSet(mass,massCategory));
	  h->Fill(mass.getVal(),evWeight.getVal());
	}
      massProfile->SetPoint(ipt,MassPointCollection[ipt].second,h->GetMean());
      massProfile->SetPointError(ipt,0,h->GetMeanError());
    }    
  f->Close();

  //RooDataHist combData("combData", "combined data",RooArgList(mass), massCategory, hmap );

  //define the combined pdf as function of the top quark mass
  RooRealVar topmass("mtop","m_{top} [GeV/c^{2}]",rangeMin,rangeMax);  
  
  //core is the resolution
  RooRealVar g_mean_slope("#mu_{G}(slope)","g_mean_slope",0.4,0.0,0.7);    
  RooRealVar g_mean_shift("#mu_{G}(intercept)","g_mean_shift",162,100,180); 
  RooFormulaVar g_mean(  "g_mean",  "(@0-172)*@1+@2",   RooArgSet(topmass,g_mean_slope,g_mean_shift));

  RooRealVar g_sigma_slope("#sigma_{G}(slope)","g_sigma_slope",0.1,0.,1.);
  RooRealVar g_sigma_shift("#sigma_{G}(intercept)","g_sigma_shift",20,0.,30);//30);
  RooFormulaVar g_sigma( "g_sigma", "(@0-172)*@1+@2",   RooArgSet(topmass,g_sigma_slope,g_sigma_shift)); 

  RooGaussian gaus("gaus", "Mass component 1", mass, g_mean, g_sigma);
  
  //   RooRealVar cb_alpha_slope("#alpha_{CB}(slope)","cb_alpha_slope",0.4,0.0,1.5);    
  //   RooRealVar cb_alpha_shift("#alpha_{CB}(intercept)","cb_alpha_shift",10,0,50); 
  //   RooFormulaVar cb_alpha(  "cb_alpha",  "(@0-172)*@1+@2",   RooArgSet(topmass,cb_alpha_slope,cb_alpha_shift));
  //   RooRealVar cb_n_slope("n_{CB}(slope)","cb_n_slope",0.4,0.0,1.5);    
  //   RooRealVar cb_n_shift("n_{CB}(intercept)","cb_n_shift",120,100,160); 
  //   RooFormulaVar cb_n(  "cb_mean",  "(@0-172)*@1+@2",   RooArgSet(topmass,cb_n_slope,cb_n_shift));
  //   RooCBShape  cball("cball", "Mass component 1", mass, g_mean, g_sigma,cb_alpha,cb_n);
 

  RooRealVar l_mean_slope("mpv_{L}(slope)","l_mean_slope",0.,0.,1.0);
  RooRealVar l_mean_shift("mpv_{L}(intercept)","l_mean_shift",140,100,250);
  RooFormulaVar l_mean(  "l_mean",  "(@0-172)*@1+@2",   RooArgSet(topmass,l_mean_slope,l_mean_shift));
  RooRealVar l_sigma_slope("#sigma_{L}(slope)","l_sigma_slope",0.,0.,1.);
  RooRealVar l_sigma_shift("#sigma_{L}(intercept)","l_sigma_shift",20,0,25);
  RooFormulaVar l_sigma( "l_sigma", "(@0-172)*@1+@2",   RooArgSet(topmass,l_sigma_slope,l_sigma_shift)); 
  RooRealVar massfrac_slope("#alpha(slope)","massfrac_slope",0.0,0,0.01);
  RooRealVar massfrac_shift("#alpha(intercept)","massfrac_shift",0.38,0.,1.0);
  RooFormulaVar massfrac( "#alpha", "(@0-172)*@1+@2", RooArgSet(topmass,massfrac_slope,massfrac_shift)); 
  RooLandau lan("lan", "Mass component 2", mass, l_mean, l_sigma);  

  RooAddPdf massmodel("model","model",RooArgList(lan,gaus),massfrac);
  //RooAddPdf massmodel("model","model",RooArgList(lan,cball),massfrac);

  
  //now split per categories with fixed top mass
  RooSimPdfBuilder builder(massmodel) ;
  RooArgSet* config = builder.createProtoBuildConfig() ;
  config->setStringValue("physModels","model");     // Name of the PDF we are going to work with
  config->setStringValue("splitCats","cat");     // Category used to differentiate sub-datasets
  config->setStringValue("model","cat : mtop");  // Prescription to taylor PDF parameters mtop for each subset in signal
  RooSimultaneous* simPdf = builder.buildPdf(*config,&combData) ;
  config = simPdf->getParameters(combData);
  for(size_t ipt=0; ipt<MassPointCollection.size(); ipt++)
    {
      TString cName("m"); cName+=(ipt+1);
      Float_t imass=MassPointCollection[ipt].second;
      (((RooRealVar &)(*config)["mtop_"+cName])).setRange(imass,imass);
      (((RooRealVar &)(*config)["mtop_"+cName])).setVal(imass);
    }
  config->writeToStream(cout,false);

  //fit to data
  simPdf->fitTo(combData,Range(rangeMin,rangeMax),SumW2Error(kTRUE));//,Hesse(),Minos());

  //display
  setStyle();
  int ny=MassPointCollection.size()/4+1;
  TString chName("");
  if(channel==OPFLAVOR) chName += "of";
  if(channel==SAMEFLAVOR) chName += "sf";
  TString cnvName("SignalPDF_"); cnvName += chName; cnvName +="_"; cnvName += nbtags; cnvName += "btags";
  TCanvas *c = getNewCanvas(cnvName,cnvName,true);
  c->SetCanvasSize(2400,ny*600);
  c->SetWindowSize(2400,ny*600);
  c->Divide(4,ny);
  for(size_t ipt=0; ipt<MassPointCollection.size(); ipt++)
    {
      TString cName("m"); cName += (ipt+1);
      
      TPad *p = (TPad *)c->cd(ipt+1);
      //p->SetGridx();
      //p->SetGridy();
      RooPlot* frame = mass.frame(Range(100,300),Bins(20));
      RooDataSet* dataslice = (RooDataSet *)combData.reduce("cat==cat::"+cName);
      dataslice->plotOn(frame,DataError(RooAbsData::SumW2));

      RooCategory theCategory(cName,cName);
      simPdf->plotOn(frame,Slice(theCategory),ProjWData(mass,*dataslice));
      frame->GetYaxis()->SetTitleOffset(1.3);
      frame->GetYaxis()->SetTitle("Events");
      frame->GetXaxis()->SetTitleOffset(1.0);
      frame->Draw();
      TPaveText *pave = new TPaveText(0.55,0.65,0.85,0.8,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      pave->SetTextFont(42);
      char titBuf[50];
      sprintf(titBuf,"m_{top}=%3.1f GeV/c^{2}",MassPointCollection[ipt].second);
      pave->AddText(titBuf);
      pave->Draw("same");

      if(ipt==0)
	{
	  TLegend *leg = p->BuildLegend();
	  formatForCmsPublic(p,leg,"CMS simulation",1);
	  leg->Delete();	  
	}

      //superimpose the average of the gaussian fit
      TLine *l=new TLine(MassPointCollection[ipt].second,0,MassPointCollection[ipt].second,frame->GetMaximum());
      l->SetLineWidth(2);
      l->SetLineColor(kGray);
      l->SetLineStyle(2);
      l->Draw();
    }

  c->cd(MassPointCollection.size()+1);
  massProfile->SetMarkerStyle(20);
  massProfile->Draw("ap");
  massProfile->GetXaxis()->SetTitle("Generated top mass [GeV/c^{2}]");
  massProfile->GetYaxis()->SetTitle("<Reconstructed top mass> [GeV/c^{2}]");
  
  c->Modified();
  c->Update();
  c->SaveAs( cnvName + TString(".C") );
  c->SaveAs( cnvName + TString(".png") );
  c->SaveAs( cnvName + TString(".pdf") );
  delete c;

  FitResults_t fitPars;
  fitPars["#mu_{G}(slope)"]=Value_t(g_mean_slope.getVal(),g_mean_slope.getError());
  fitPars["#mu_{G}(intercept)"]=Value_t(g_mean_shift.getVal(),g_mean_shift.getError());
  fitPars["#sigma_{G}(slope)"]=Value_t(g_sigma_slope.getVal(),g_sigma_slope.getError());
  fitPars["#sigma_{G}(intercept)"]=Value_t(g_sigma_shift.getVal(),g_sigma_shift.getError());
  // fitPars["#alpha_{CB}(slope)"]=Value_t(cb_alpha_slope.getVal(),cb_alpha_slope.getError());
  // fitPars["#alpha_{CB}(intercept)"]=Value_t(cb_alpha_shift.getVal(),cb_alpha_shift.getError());
  //  fitPars["n_{CB}(slope)"]=Value_t(cb_n_slope.getVal(),cb_n_slope.getError());
  // fitPars["n_{CB}(intercept)"]=Value_t(cb_n_shift.getVal(),cb_n_shift.getError());
  fitPars["mpv_{L}(slope)"]=Value_t(l_mean_slope.getVal(),l_mean_slope.getError());
  fitPars["mpv_{L}(intercept)"]=Value_t(l_mean_shift.getVal(),l_mean_shift.getError());
  fitPars["#sigma_{L}(slope)"]=Value_t(l_sigma_slope.getVal(),l_sigma_slope.getError());
  fitPars["#sigma_{L}(intercept)"]=Value_t(l_sigma_shift.getVal(),l_sigma_shift.getError());
  //  fitPars["mpv_{L2}(slope)"]=Value_t(l2_mean_slope.getVal(),l2_mean_slope.getError());
  //  fitPars["mpv_{L2}(intercept)"]=Value_t(l2_mean_shift.getVal(),l2_mean_shift.getError());
  //  fitPars["#sigma_{L2}(slope)"]=Value_t(l2_sigma_slope.getVal(),l2_sigma_slope.getError());
  //  fitPars["#sigma_{L2}(intercept)"]=Value_t(l2_sigma_shift.getVal(),l2_sigma_shift.getError());
  fitPars["#alpha(slope)"]=Value_t(massfrac_slope.getVal(),massfrac_slope.getError());
  fitPars["#alpha(intercept)"]=Value_t(massfrac_shift.getVal(),massfrac_shift.getError());
  //  fitPars["#beta(slope)"]=Value_t(massfrac2_slope.getVal(),massfrac2_slope.getError());
  //  fitPars["#beta(intercept)"]=Value_t(massfrac2_shift.getVal(),massfrac2_shift.getError());
  return fitPars;
}




//  LocalWords:  hGexpmodel
