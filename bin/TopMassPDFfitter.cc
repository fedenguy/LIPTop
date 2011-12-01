/**
   @short  run me as
   .L TopMassPDFfitter.C+ 
   
   BckgPDFs(file_withmass_trees);
   SignalPDFs("EventSummary.root"," && evmeasurements[4]==0");
   SignalPDFs("EventSummary.root"," && evmeasurements[4]==1");
   SignalPDFs("EventSummary.root"," && evmeasurements[4]>1");
*/

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

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include<vector>

using namespace RooFit ;


RooSimultaneous *SignalPDFs(TString url="EventSummaries.root",int nbtags=-1);

//
int main(int argc, char* argv[])
{
  // load framework libraries                                                                                                                                                            
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  TString url=argv[1];
  int nbtags(-1);
  if(argc>2) sscanf(argv[2],"%d",&nbtags);
  SignalPDFs(url,nbtags);
}


//
// void BckgPDFs(TString url="MassPDFs.root",TString var="Mass")
// {
//   //get the histograms and build the sum
//   typedef std::pair<string,double> proc_t;
//   typedef std::vector<proc_t> sample_t;
//   std::map<string,sample_t> bckgSamples;
//   std::map<string, TH1D *> bckgDists;

//   double lum=36.1;
//   sample_t ttbar; 
//   ttbar.push_back(proc_t("otherttbar",1.30) ); 
//   bckgSamples["t#bar{t} (other)"] = ttbar;

//   sample_t ewk; 
//   double ewk_syst(0.0);
//   //  double ewk_syst(0.3);
//   ewk.push_back(proc_t("wjets",1.23*(1+ewk_syst)));
//   ewk.push_back(proc_t("ww",1.0274*(1+ewk_syst)));
//   ewk.push_back(proc_t("zz",0.09662*(1+ewk_syst)));
//   ewk.push_back(proc_t("wz",0.156257*(1+ewk_syst)));
//     bckgSamples["EWK"]=ewk;

//   sample_t stop; 
//   double stop_syst(0.3);
//   //double stop_syst(0);
//   stop.push_back(proc_t("singletop_tW",2.944*(1+stop_syst)) ); 
//   bckgSamples["Single top"] = stop;

//   setTDRStyle();
//   TCanvas *c = new TCanvas("bckgpdfs","Background PDFs");
//   TFile *f = TFile::Open(url);
//   int totalevts=0;
//   for(std::map<string, sample_t>::iterator bckgIt = bckgSamples.begin();
//       bckgIt != bckgSamples.end();
//       bckgIt++)
//     {
//       TH1D *sampleH=0;
//       for(sample_t::iterator pIt = bckgIt->second.begin();
// 	  pIt != bckgIt->second.end();
// 	  pIt++)
// 	{
// 	  TString tname=pIt->first+"/mass";
// 	  TTree *t = (TTree *) f->Get(tname);

// 	  if(t==0) continue;
// 	  if(var=="Mass") t->Draw("Mass[0]>>hmass(20,100,500)","IsSignal==0 && Mass[0]>0");
// 	  else t->Draw("MT2[0]>>hmass(20,100,500)","IsSignal==0 && Mass[0]>0");
// 	  TH1D *h = (TH1D*)gDirectory->Get("hmass");
// 	  if(h==0) continue;
// 	  if(sampleH==0)     
// 	    {
// 	      TString newName("hmass_"); newName += pIt->first;
// 	      sampleH=(TH1D *)h->Clone(newName);
// 	      sampleH->SetDirectory(0);
// 	      sampleH->Sumw2();
// 	      sampleH->SetTitle(TString(bckgIt->first));
// 	      sampleH->SetLineColor(1);    
// 	      sampleH->SetMarkerColor(1);  
// 	      sampleH->SetFillColor(1);
// 	      sampleH->SetMarkerStyle(24); 
// 	      sampleH->SetFillStyle(3003); 
// 	      sampleH->GetYaxis()->SetTitle("Events / (20 GeV/c^{2})");
// 	      sampleH->GetXaxis()->SetTitle("Reconstructed " + var + " [GeV/c^{2}]");
// 	      sampleH->Reset("ICE");
// 	      bckgDists[bckgIt->first]=sampleH;
// 	    }	  
// 	  if(h->Integral()==0) continue;
// 	  totalevts+=h->Integral();
// 	  sampleH->Add(h,pIt->second/h->Integral());
// 	}
//     }
//   f->Close();

//   TH1D *sumH= (TH1D* )bckgDists["Single top"]->Clone("totalbckg");
//   sumH->Reset("ICE");
//   sumH->SetDirectory(0);
//   sumH->Add(bckgDists["t#bar{t} (other)"]);
//   sumH->Add(bckgDists["Single top"]);
//   sumH->Add(bckgDists["EWK"]);
//   sumH->Scale(totalevts/sumH->Integral());

//   //prepare the pdfs for the fit in mass
//   RooRealVar mass("m",var,100,var=="Mass"? 500:300,"GeV/c^{2}");
//   RooRealVar mpv_l("mpv_{l}","Mpv of landau",100,300);
//   RooRealVar sigma_l("#sigma_{l}","Sigma of landau",10,30);
//   RooLandau lan("blandau","Mass component 1",mass,mpv_l,sigma_l);
//   RooRealVar mean_g("#mu_{g}","Mean of gaus",100,200);
//   RooRealVar sigma_g("#sigma_{g}","Sigma of gaus",10,30);
//   RooGaussian gaus("bgauss","Mass component 2",mass,mean_g,sigma_g);   
//   RooRealVar massfrac("#alpha","Fraction of component 1",0.6,1.);
//   RooAddPdf massmodel("model","Model",RooArgList(lan,gaus),massfrac);
//   RooDataHist dh("dh","dh",mass,sumH);       
//   massmodel.fitTo(dh,Save(kTRUE),Range(100.,500.),SumW2Error(kTRUE));
//   RooPlot* frame = mass.frame(Title("Combined"));
//   dh.plotOn(frame,DrawOption("pz"),DataError(RooAbsData::SumW2)); 
//   massmodel.plotOn(frame,DrawOption("F"),FillColor(kAzure-4),FillStyle(3001),MoveToBack(),Range(100,500)) ;

//   c->Clear();
//   c->SetWindowSize(1600,400);
//   c->Divide(4);
//   int ibckg(0);
//   for(std::map<std::string, TH1D *>::iterator it = bckgDists.begin();
//       it != bckgDists.end();
//       it++,ibckg++)
//     {
//       TPad *p = (TPad *) c->cd(ibckg+1);  
//       p->SetGridx();  p->SetGridy();
//       it->second->Draw("hist");  
//       it->second->GetXaxis()->SetTitleOffset(1.0);
//       it->second->GetYaxis()->SetTitleOffset(1.0);
//       TPaveText *pt = new TPaveText(0.75,0.85,0.97,0.95,"brNDC");
//       pt->SetBorderSize(0);
//       pt->SetFillColor(0);
//       pt->SetFillStyle(0);
//       pt->AddText(it->first.c_str());
//       pt->Draw();
//     }
//   c->cd(4);
//   frame->Draw();
//   frame->GetYaxis()->SetTitleOffset(1.0);
//   frame->GetYaxis()->SetTitle("Events (A.U.)");
//   frame->GetXaxis()->SetTitleOffset(0.8);
//   frame->GetXaxis()->SetTitle("Reconstructed Mass [GeV/c^{2}]");
// }


//
RooSimultaneous *SignalPDFs(TString url,int nbtags)
{
  gStyle->SetOptStat(0);

  //the mass points
  typedef std::pair<TString,Float_t> MassPoint_t;
  std::vector<MassPoint_t> MassPointCollection;
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_161v5",161.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_163v5",163.5) );
  //  MassPointCollection.push_back( MassPoint_t("TTJets_mass_166v5",166.5) );
  MassPointCollection.push_back( MassPoint_t("TTJets_mass_169v5",169.5) ); 
  MassPointCollection.push_back( MassPoint_t("TTJets",           172.5) );
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
  TCanvas *c = 0;
  for(size_t ipt=0; ipt<MassPointCollection.size(); ipt++)
    {
      if(ipt%5==0)
	{
	  if(c!=0)
	    {
	      c->SaveAs( c->GetName() + TString(".C") );
	      c->SaveAs( c->GetName() + TString(".png") );
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
      c->SaveAs( c->GetName() + TString(".C") );
      c->SaveAs( c->GetName() + TString(".png") );
      delete c;
    }

  return simPdf;
}



