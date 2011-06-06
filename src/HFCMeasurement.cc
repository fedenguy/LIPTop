#include "LIP/Top/interface/HFCMeasurement.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "TFile.h"

//
void HFCMeasurement::bookMonitoringHistograms(int maxJets)
{
  biasMonH = new TH1D("bias",";bias=#varepsilon_{b}-#bar{#varepsilon_{b}};Pseudo-experiments",100,-0.99,1.01);
  formatPlot( biasMonH,1,1,1,20,0,true,true,1,1,1);
  
  pullMonH = new TH1D("pull",";pull=(#varepsilon_{b}-#bar{#varepsilon_{b}}) / #sigma_{#varepsilon_{b}};Pseudo-experiments",100,-2.97,3.03);
  formatPlot(pullMonH,1,1,1,20,0,true,true,1,1,1);      

  statMonH = new TH1D("stat",";#sigma_{#varepsilon_{b}};Pseudo-experiments",100,0.0,1.0);
  formatPlot(pullMonH,1,1,1,20,0,true,true,1,1,1);      
  
  for(int ijets=2; ijets<=maxJets; ijets++)
    {
      TString tag("_"); tag += ijets; tag+= "jets";
      TString title(""); title += ijets; title += " jets";
      TH1D *bmh=new TH1D("btag"+tag,title+";b-tag multiplicity;Events",maxJets+1,0,maxJets+1);
      formatPlot(bmh,1,1,1,20,0,true,true,1,1,1);
      bmultH.push_back(bmh);
    }
}


//
void HFCMeasurement::showMonitoringHistograms(bool debug)
{
  //open file
  TFile *fout=TFile::Open("HFCMeasurement.root","RECREATE");
  fout->cd();
  biasMonH->Fit("gaus");  biasMonH->Write(); 
  pullMonH->Fit("gaus");  pullMonH->Write();
  statMonH->Write();
  
  if(debug)
    {
      TCanvas *c = getNewCanvas("monh","monh",false);
      c->SetWindowSize(1500,500);
      c->Divide(3,1);
      c->cd(1);
      biasMonH->Draw("e1");
      c->cd(2);
      statMonH->Draw("e1");
      c->cd(3);
      pullMonH->Draw("e1");
    }
  
  //close file
  fout->Close();
  fout->Delete();
}


//
void HFCMeasurement::initHFCModel(int maxJets,TString wp,int fitType)
{
  //global inputs
  model.bmult = new RooRealVar("bmult","N_{btags}",0.,4.);
  if(fitType==0) model.r = new RooRealVar("r","R",1.0);
  else           model.r = new RooRealVar("r","R",1.0,0.0,2.0);
  model.lfacceptance = new RooRealVar("lfacceptance","A(R=0)",1.0);

  if(fitType==0) model.eb = new RooRealVar("eb","#varepsilon_{b}",btageff[wp],0.0,1.0);
  else           
    {
      model.eb = new RooRealVar("eb","#varepsilon_{b}",btageff[wp],0.0,1.0);
      model.eb_mean_constrain = new RooRealVar("meaneb","#bar{#varepsilon_{b}}",btageff[wp]);
      model.eb_sigma_constrain = new RooRealVar("unceb","#sigma_{#varepsilon_{b}}",btageffunc[wp]);
      model.eb_constrain = new RooGaussian("eb_constrain","#varepsilon_{b} constrain",*model.eb,*model.eb_mean_constrain,*model.eb_sigma_constrain);
    }
  model.eq = new RooRealVar("eq","#varepsilon_{q}",mistagrate[wp],0.0,1.0); 
  model.eq_mean_constrain = new RooRealVar("meaneq","#bar{#varepsilon_{q}}",mistagrate[wp]);
  model.eq_sigma_constrain = new RooRealVar("unceq","#sigma_{#varepsilon_{q}}",mistagrateunc[wp]);
  model.eq_constrain = new RooGaussian("eq_constrain","#varepsilon_{q} constrain",*model.eq,*model.eq_mean_constrain,*model.eq_sigma_constrain);

  //event type composition (alpha_i)
  model.alpha2=new RooRealVar("alpha2","#alpha_{2}",0.5,0.7);
  model.alpha2_mean_constrain = new RooRealVar("meanalpha2","#bar{#alpha}_{2}",0.63);
  model.alpha2_sigma_constrain = new RooRealVar("uncalpha2","#sigma_{#alpha_{2}}",sqrt(pow(0.01,2)+pow((0.03+0.01)*0.5,2)));
  model.alpha2_constrain = new RooGaussian("ctr_alpha2","#alpha_{2} constrain",*model.alpha2,*model.alpha2_mean_constrain,*model.alpha2_sigma_constrain);
  
  model.alpha0 = new RooRealVar("alpha0","#alpha_{0}",0,0.2);
  model.alpha0_mean_constrain = new RooRealVar("meanalpha0","#bar{#alpha}_{0}",0.135);
  model.alpha0_sigma_constrain = new RooRealVar("uncalpha0","#sigma_{#alpha_{0}}",sqrt(pow(0.007,2)+pow((0.006+0.003)/2,2)));
  model.alpha0_constrain = new RooGaussian("ctralpha0","#alpha_{0} constrain",*model.alpha0,*model.alpha0_mean_constrain,*model.alpha0_sigma_constrain);

  model.alpha1 = new RooFormulaVar("alpha1","1-@0-@1",RooArgSet(*model.alpha2,*model.alpha0));

  //add exclusive jet multiplicity models
  for(int jmult=2; jmult<=maxJets; jmult++)
    {
      TString tag("_"); tag+=jmult;      
      
      //create the model
      HeavyFlavorPDF *mhfc = new HeavyFlavorPDF("hfcmodel"+tag,"hfcmodel"+tag,*model.bmult,*model.r,*model.eb,*model.eq,*model.alpha2,*model.alpha1,*model.alpha0);
      mhfc->setJetMultiplicity(jmult);
      model.pdfSet.add( *mhfc );
      
      //add the constrains
      RooProdPdf *modelconstr=0;
      //if(fitType==0) modelconstr=new RooProdPdf("modelconstr"+tag,"model x product of constrains",RooArgSet(*mhfc,*model.alpha2_constrain,*model.alpha0_constrain,*model.eq_constrain));
      if(fitType==0) modelconstr=new RooProdPdf("modelconstr"+tag,"model x product of constrains",RooArgSet(*mhfc,*model.alpha2_constrain,*model.alpha0_constrain));
      else           modelconstr=new RooProdPdf("modelconstr"+tag,"model x product of constrains",RooArgSet(*mhfc,*model.alpha2_constrain,*model.alpha0_constrain,*model.eq_constrain,*model.eb_constrain));
      model.constrPDFSet.add( *modelconstr );
    }
  
}

//
void HFCMeasurement::fitHFCtoEnsemble(TTree *t, EventSummary_t *evt, bool debug)
{
  if(t==0 || evt==0) return;
  
  //reset histos
  for(std::vector<TH1D *>::iterator it = bmultH.begin();
      it != bmultH.end();
      it++) 
    (*it)->Reset("ICE");
  
  //run over the events
  for(unsigned int i=0; i<t->GetEntriesFast(); i++)
    {
      //veto if category is exclusive (emu  or ee+mumu channels)
      t->GetEntry(i);
      if(evt->cat != dilepton::EMU) continue;
      //if((evt->categ != 0xbb && evt->categ != 0xdd)) continue;

      //count the number of btags
      int nbtags(0),njets(0);
      for(int ipart=0; ipart<evt->nparticles; ipart++)
	{
	  int id = evt->id[ipart];
	  if(id !=1 ) continue;
	  njets++;
	  if(wp_=="tight")
	    {
	      if( evt->info2[ipart] > btagcut[wp_] ) nbtags++;
	    }
	  else
	    {
	      if( evt->info1[ipart] > btagcut[wp_] ) nbtags++;
	    }
	}
      if(njets<2) continue;
      if(njets>maxJets_) continue;
      bmultH[njets-2]->Fill(nbtags);
    }

  if(debug)
    { 
      setStyle();
      TCanvas *c=getNewCanvas("ensemble","ensemble",false);
      c->cd();
      c->SetWindowSize(maxJets_*500,500);
      c->SetGridx();
      c->SetGridy();     
      c->Divide(maxJets_,1);
      for(int njets=2; njets<=maxJets_; njets++)
	{
	  c->cd(njets-1);
	  bmultH[njets-2]->Draw("hist");
	}
    }

  //run the fit
  fitHFCto(bmultH,debug);
}


//
void HFCMeasurement::fitHFCto(std::vector<TH1D *> &bmultH, bool debug)
{
  using namespace RooFit;

  RooArgSet allLL,jetOccupancies;  
  RooDataHist* alldata = new RooDataHist("data","data", RooArgList(*model.bmult));
  TIterator *pdfIt = model.constrPDFSet.createIterator();
  for(int jmult=2; jmult<=maxJets_; jmult++)
    {
      TString tag("_"); tag+=jmult;      
      TH1D *h = bmultH[jmult-2];
      RooDataHist* ds = new RooDataHist("data"+tag,"data"+tag, RooArgList(*model.bmult), h);
      alldata->add(*ds);
      RooRealVar *jetocc = new RooRealVar("jetocc"+tag,"occ"+tag,h->Integral());
      jetOccupancies.add(*jetocc);
      RooProdPdf *modelconstr = (RooProdPdf *) pdfIt->Next();
      RooAbsReal *nll=0;
      if(fitType_==0) nll = modelconstr->createNLL(*ds,Constrain(RooArgSet(*model.alpha2,*model.alpha0,*model.eq_constrain)));
      //if(fitType_==0) nll = modelconstr->createNLL(*ds,Constrain(RooArgSet(*model.alpha2,*model.alpha0)));
      else            nll = modelconstr->createNLL(*ds,Constrain(RooArgSet(*model.alpha2,*model.alpha0,*model.eq_constrain,*model.eb_constrain)));
      allLL.add(*nll);      
    }
  
  //maximize combined likelihood 
  RooAddition *combll = new RooAddition("combll","combll",allLL);      
  RooMinuit minuit(*combll); 
  minuit.migrad();
  minuit.setErrorLevel(0.5);
  minuit.hesse();
  RooFitResult *r=minuit.save();

  //draw the result
  if(debug)
    {
      TCanvas *c=getNewCanvas("combination","combination",false);
      c->cd();
      c->SetWindowSize(600,600);
      c->SetGridx();
      c->SetGridy();     

      //likelihood (main frame)
      RooPlot *frame = 0;
      if(fitType_==0) 
	{
	  frame = model.eb->frame(Title("Likelihood"),Range(0.,1.));
	  frame->GetXaxis()->SetTitle("#varepsilon_{b}");
	}
      else
	{
	  frame = model.r->frame(Title("Likelihood"),Range(0,1.3)) ;    
	  frame->GetXaxis()->SetTitle("R");
	}
      combll->plotOn(frame,ShiftToZero(),Name("ll"));

      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetYaxis()->SetTitle("-log(L/L_{max})");
      frame->GetYaxis()->SetTitleOffset(1);
      frame->Draw();
      
      //fit result

      TPaveText *pt = new TPaveText(0.15,0.75,0.45,0.95,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetFillStyle(0);
      pt->SetTextAlign(10);
      pt->AddText("CMS preliminary");
      //      pt->AddText("36.1 pb^{-1} at #sqrt{s}=7 TeV");
      char buf[100];
      if(fitType_==0)  sprintf(buf,"#varepsilon_{b}=%3.3f #pm %3.2f",model.eb->getVal(),model.eb->getError());
      else sprintf(buf,"R=%3.2f #pm %3.2f",model.r->getVal(),model.r->getError());
      pt->AddText(buf);
      pt->Draw();

      //the model fit
      TPad *npad = new TPad("llpad","ll", 0.6, 0.6, 0.9, 0.9);
      npad->Draw();
      npad->cd();
      frame = model.bmult->frame();

      RooAddPdf incmodel("incmodel","Inclusive Model",model.pdfSet,jetOccupancies);
      alldata->plotOn(frame,Binning(4),DrawOption("pz"));
      incmodel.plotOn(frame,Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
      frame->Draw();
      /*
      c = new TCanvas("contour","contour");
      c->cd();
      c->SetWindowSize(600,600);
      c->SetGridx();
      c->SetGridy();     
      RooPlot *plot = minuit.contour(*model.eb,*model.eq,1,2,3) ;
      plot->SetTitle("Contour for 1s,2s,3s between eb and eq") ;
      plot->Draw();
      */
    }
  cout << "Here2 ?" << endl;
  float bias=(fitType_==0 ? model.eb->getVal()-btageff[wp_] : model.r->getVal()-1.0);
  float unc =(fitType_==0 ? model.eb->getError() :  model.r->getError());
  if(unc>0.01)
    {
      biasMonH->Fill(bias);
      statMonH->Fill(unc);
      pullMonH->Fill(bias/unc);
 
      // Access basic information
      cout << "EDM = " << r->edm() << endl
	   << "-log(L) at minimum = " << r->minNll() << endl   
	   << "final value of floating parameters" << endl ;
      r->floatParsFinal().Print("s") ;
    }
  cout << "Somewhere" << endl;
}
