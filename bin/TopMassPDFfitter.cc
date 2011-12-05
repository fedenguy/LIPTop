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

#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include<vector>

typedef std::pair<float,float> Value_t;
typedef std::map<std::string,Value_t> FitResults_t;

using namespace std; 
using namespace RooFit ;


FitResults_t SignalPDFs(TString url="EventSummaries.root",int nbtags=-1);
FitResults_t BckgPDFs(TString url="MassPDFs.root");

//
int main(int argc, char* argv[])
{
  // load framework libraries                                                                                                                                                            
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  TString surl=argv[1];
  TString burl=argv[2];
  
  //fit the templates
  std::map<string,FitResults_t> allFitResults;
  allFitResults["[Signal : 0 b-tags]"]       = SignalPDFs(surl,0);
  allFitResults["[Signal : 1 b-tags]"]       = SignalPDFs(surl,1);
  allFitResults["[Signal : 2 b-tags]"]       = SignalPDFs(surl,2);
  allFitResults["[Background : inclusive]"]  = BckgPDFs(burl);

  //display the results
  cout << " *************** TopMassPDFfitter  *********************** " << endl;
  for(std::map<string,FitResults_t>::iterator it = allFitResults.begin();
      it != allFitResults.end();
      it++)
    {
      cout << it->first << endl;
      for(FitResults_t::iterator itt = it->second.begin();
	  itt != it->second.end();
	  itt++)
	cout << itt->first << " : " << itt->second.first << " +/- " << itt->second.second << endl;
      cout << endl;
    }
}


//
FitResults_t BckgPDFs(TString url)
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
  ewkSamples.push_back("WJetsToLNu");
  bckgProcs["EWK"]=ewkSamples;

  std::vector<string>dySamples;
  dySamples.push_back("DYJetsToLL");
  dySamples.push_back("DYJetsToMuMu_M20to50");
  dySamples.push_back("DYJetsToEE_M20to50");
  bckgProcs["DY"]=dySamples;

  //build the mass histograms
  int ipt(1);
  TFile *f = TFile::Open(url);
  for(std::map<string,std::vector<string> >::iterator it = bckgProcs.begin(); it != bckgProcs.end(); it++,ipt++)
    {
      TString sName("m"); sName += (ipt+1);      
      TH1D *h= new TH1D(sName,sName,80,100,500);      
      h->SetDirectory(0);
      h->GetYaxis()->SetTitle("Events / (5 GeV/c^{2})");
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

	  Float_t evmeasurements[10];
	  t->GetBranch("evmeasurements")->SetAddress(evmeasurements);
	  for(Int_t i=0; i<t->GetEntriesFast(); i++)
	    {
	      t->GetEntry(i);
	      if(evmeasurements[0]>0)  h->Fill(evmeasurements[0]);
	    }
	}
    }    
  f->Close();
  
  //the sum
  TH1D *sumH= (TH1D* )bckgDists["Top"]->Clone("totalbckg");
  sumH->SetDirectory(0);
  sumH->Add(bckgDists["EWK"]);
  //  sumH->Add(bckgDists["DY"]);
  
  //prepare the pdfs for the fit in mass
  RooRealVar mass("m","Mass",100,500,"GeV/c^{2}");
  RooRealVar mpv_l("mpv_{l}","Mpv of landau",100,300);
  RooRealVar sigma_l("#sigma_{l}","Sigma of landau",10,30);
  RooLandau lan("blandau","Mass component 1",mass,mpv_l,sigma_l);
  RooRealVar mean_g("#mu_{g}","Mean of gaus",100,200);
  RooRealVar sigma_g("#sigma_{g}","Sigma of gaus",10,30);
  RooGaussian gaus("bgauss","Mass component 2",mass,mean_g,sigma_g);   
  RooRealVar massfrac("#alpha","Fraction of component 1",0.6,1.);
  RooAddPdf massmodel("model","Model",RooArgList(lan,gaus),massfrac);
  RooDataHist dh("dh","dh",mass,sumH);       
  massmodel.fitTo(dh,Save(kTRUE),Range(100.,500.),SumW2Error(kTRUE));

  //show the results of the fit
  setStyle();
  TCanvas *c = new TCanvas("bckgpdfs","Background PDFs");
  c->SetWindowSize(1600,400);
  c->Divide(4);

  c->cd(1);
  RooPlot* frame = mass.frame(Title("Combined"));
  dh.plotOn(frame,DrawOption("pz"),DataError(RooAbsData::SumW2)); 
  massmodel.plotOn(frame,DrawOption("F"),FillColor(kAzure-4),FillStyle(3001),MoveToBack(),Range(100,500)) ;
  frame->Draw();
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitle("Events (A.U.)");
  frame->GetXaxis()->SetTitleOffset(0.8);
  frame->GetXaxis()->SetTitle("Reconstructed Mass [GeV/c^{2}]");

  
  int ibckg(1);
  for(std::map<std::string, TH1D *>::iterator it = bckgDists.begin();
      it != bckgDists.end();
      it++,ibckg++)
    {
      c->cd(ibckg+1);
      it->second->Rebin();
      it->second->Draw("hist");  
      it->second->GetXaxis()->SetTitleOffset(1.0);
      it->second->GetYaxis()->SetTitleOffset(1.0);
      TPaveText *pt = new TPaveText(0.75,0.85,0.97,0.95,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetFillStyle(0);
      pt->AddText(it->first.c_str());
      pt->Draw();
    }
  c->SaveAs("bckg_masstemplate.C");
  c->SaveAs("bckg_masstemplate.png");
  
  //return the result
  FitResults_t fitPars;
  fitPars["mpv_{l}"]    = Value_t(mpv_l.getVal(),mpv_l.getError());
  fitPars["#sigma_{l}"] = Value_t(sigma_l.getVal(),sigma_l.getError());
  fitPars["#mu_{g}"]    = Value_t(mean_g.getVal(),mean_g.getError());
  fitPars["#sigma_{g}"] = Value_t(sigma_g.getVal(),sigma_g.getError());
  fitPars["#alpha"]     = Value_t(massfrac.getVal(),massfrac.getError()); 
  return fitPars;
}



//
FitResults_t SignalPDFs(TString url,int nbtags)
{

  //the mass points
  typedef std::pair<TString,Float_t> MassPoint_t;
  std::vector<MassPoint_t> MassPointCollection;
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_161v5",161.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_163v5",163.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_166v5",166.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_169v5",169.5) ); 
  MassPointCollection.push_back( MassPoint_t("TTJets",    172.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_175v5",175.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_178v5",178.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass181v5", 181.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass184v5", 184.5) );

  std::map<std::string,TH1*> hmap;
  
  //get pdfs from file
  RooCategory sample("signal","") ;
  TFile *f = TFile::Open(url);
  for(size_t ipt=0; ipt<MassPointCollection.size(); ipt++)
    {
      TString sName("m"); sName += (ipt+1);      

      TString tname=MassPointCollection[ipt].first + "/data";
      TTree *t = (TTree *) f->Get(tname);
      if(t==0) continue;

      //fill the mass histogram
      TH1D *h= new TH1D(sName,sName,80,100,500);      

      Float_t evmeasurements[10];
      t->GetBranch("evmeasurements")->SetAddress(evmeasurements);
      for(Int_t i=0; i<t->GetEntriesFast(); i++)
	{
	  t->GetEntry(i);
	  if(nbtags==0 && evmeasurements[4]!=0) continue;
	  if(nbtags==1 && evmeasurements[4]!=1) continue;
	  if(nbtags==2 && evmeasurements[4]<1) continue;
	  h->Fill(evmeasurements[0]);
	}

      h->SetDirectory(0);
      h->GetYaxis()->SetTitle("Events / (5 GeV/c^{2})");
      h->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
      h->GetXaxis()->SetTitleOffset(0.8);
      h->GetYaxis()->SetTitleOffset(0.8);
      char titbuf[20];
      sprintf(titbuf,"m=%3.1f",MassPointCollection[ipt].second);
      h->SetTitle(titbuf);

      sample.defineType(TString(sName));
      hmap[sName.Data()] = h;
    }    
  f->Close();

  // divide the binned data in categories according to the generated top quark mass
  RooRealVar mass("m","Mass", 100, 500);
  RooDataHist combData("combData", "combined data",mass, sample, hmap );

  //the parameters to fit and the variable
  RooRealVar g_mean_slope("#mu_{G}(slope)","g_mean_slope",0.01,0.,1.);    
  RooRealVar g_mean_shift("#mu_{G}(intercept)","g_mean_shift",162,100,180); 
  RooRealVar g_sigma_slope("#sigma_{G}(slope)","g_sigma_slope",0.01,0.,1.);
  RooRealVar g_sigma_shift("#sigma_{G}(intercept)","g_sigma_shift",10,0.,25);
  RooRealVar l_mean_slope("mpv_{L}(slope)","l_mean_slope",0.,0.,1.);//1,0,10);
  RooRealVar l_mean_shift("mpv_{L}(intercept)","l_mean_shift",212,150,250); 
  RooRealVar l_sigma_slope("#sigma_{L}(slope)","l_sigma_slope",0.,0.,1.);//1,0,10);
  RooRealVar l_sigma_shift("#sigma_{L}(intercept)","l_sigma_shift",10,0,25);
  RooRealVar massfrac_slope("#alpha(slope)","massfrac_slope",0,0,0.01);
  RooRealVar massfrac_shift("#alpha(intercept)","massfrac_shift",0.38,0.,1.);

  //build the prototype pdf
  RooRealVar    topmass( "mtop","mtop",100,300);
  RooFormulaVar g_mean(  "g_mean",  "(@0-172)*@1+@2",   RooArgSet(topmass,g_mean_slope,g_mean_shift));
  RooFormulaVar g_sigma( "g_sigma", "(@0-172)*@1+@2", RooArgSet(topmass,g_sigma_slope,g_sigma_shift)); 
  RooGaussian gaus("gaus", "Mass component 1", mass, g_mean, g_sigma);
  RooFormulaVar l_mean(  "l_mean",  "(@0-172)*@1+@2",   RooArgSet(topmass,l_mean_slope,l_mean_shift));
  RooFormulaVar l_sigma( "l_sigma", "(@0-172)*@1+@2", RooArgSet(topmass,l_sigma_slope,l_sigma_shift)); 
  RooLandau lan("lan", "Mass component 2", mass, l_mean, l_sigma);  
  RooFormulaVar massfrac( "#alpha", "(@0-172)*@1+@2", RooArgSet(topmass,massfrac_slope,massfrac_shift)); 
  RooAddPdf massmodel("model","Model",RooArgList(lan,gaus),massfrac);
  //RooNumConvPdf massmodel("model","Model",topmass,lan,gaus);

  //now split per categories
  RooSimPdfBuilder builder(massmodel) ;
  RooArgSet* config = builder.createProtoBuildConfig() ;
  config->setStringValue("physModels","model");     // Name of the PDF we are going to work with
  config->setStringValue("splitCats","signal");     // Category used to differentiate sub-datasets
  config->setStringValue("model","signal : mtop");  // Prescription to taylor PDF parameters mtop for each subset in signal
  RooSimultaneous* simPdf = builder.buildPdf(*config,&combData) ;
  config = simPdf->getParameters(combData);
  for(size_t ipt=0; ipt<MassPointCollection.size(); ipt++)
    {
      TString sName("m"); sName+=(ipt+1);
      Float_t imass=MassPointCollection[ipt].second;
      (((RooRealVar &)(*config)["mtop_"+sName])).setRange(imass,imass);
      (((RooRealVar &)(*config)["mtop_"+sName])).setVal(imass);
    }
  
  //fit to data
  simPdf->fitTo(combData,Range(100.,400.));
 
  //display
  setStyle();
  TCanvas *c = 0;
  for(size_t ipt=0; ipt<MassPointCollection.size(); ipt++)
    {
      if(ipt%5==0)
	{
	  if(c!=0)
	    {
	      TString cname=c->GetName(); cname += "_"; cname += nbtags; cname+="btags"; 
	      c->SaveAs( cname + TString(".C") );
	      c->SaveAs( cname + TString(".png") );
	    }
	  TString name("SignalPDFs_");  name+=ipt;
	  c = new TCanvas(name,name);
	  c->SetBorderSize(0);
	  c->SetFillStyle(0);
	  c->SetFillColor(0);
	  c->SetWindowSize(1750,350);
	  c->Clear();
	  c->Divide(5,1);	  
	}
      
      TPad *p = (TPad *)c->cd(ipt%5+1);
      p->SetGridx();
      p->SetGridy();
      TString procName("m"); procName += (ipt+1);
      char buf[100];
      sprintf(buf,"m_{t}=%3.1f GeV/c^{2}",MassPointCollection[ipt].second);
      RooPlot* frame = mass.frame(Title(buf));
      RooDataSet* dataslice = (RooDataSet *)combData.reduce("signal==signal::"+procName);
      dataslice->plotOn(frame,DataError(RooAbsData::SumW2));
      RooCategory newCat(procName,procName);
      simPdf->plotOn(frame,Slice(newCat),ProjWData(mass,*dataslice));
      frame->GetYaxis()->SetTitleOffset(1.0);
      frame->GetYaxis()->SetTitle("Events");
      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetXaxis()->SetTitle("Reconstructed Mass [GeV/c^{2}]");
      frame->Draw();
             
      TPaveText *pt = new TPaveText(0.75,0.85,0.97,0.95,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetFillStyle(0);
      char buf2[50];
      sprintf(buf2,"%3.1f GeV/c^{2}",MassPointCollection[ipt].second);
      pt->AddText(buf2);
      pt->Draw();
    }

  //save last canvas
  if(c!=0)
    {
      TString cname=c->GetName(); cname += "_"; cname += nbtags; cname+="btags"; 
      c->SaveAs( cname + TString(".C") );
      c->SaveAs( cname + TString(".png") );
      delete c;
    }


  //return the results
  FitResults_t fitPars;
  fitPars["#mu_{G}(slope)"]=Value_t(g_mean_slope.getVal(),g_mean_slope.getError());
  fitPars["#mu_{G}(intercept)"]=Value_t(g_mean_shift.getVal(),g_mean_shift.getError());
  fitPars["#sigma_{G}(slope)"]=Value_t(g_sigma_slope.getVal(),g_sigma_slope.getError());
  fitPars["#sigma_{G}(intercept)"]=Value_t(g_sigma_shift.getVal(),g_sigma_shift.getError());
  fitPars["mpv_{L}(slope)"]=Value_t(l_mean_slope.getVal(),l_mean_slope.getError());
  fitPars["#sigma_{L}(intercept)"]=Value_t(l_sigma_shift.getVal(),l_sigma_shift.getError());
  fitPars["#alpha(slope)"]=Value_t(massfrac_slope.getVal(),massfrac_slope.getError());
  fitPars["#alpha(intercept)"]=Value_t(massfrac_shift.getVal(),massfrac_shift.getError());
  return fitPars;
}



