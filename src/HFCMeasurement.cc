#include "LIP/Top/interface/HFCMeasurement.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

using namespace RooFit;
using namespace std;

//
void HFCMeasurement::bookMonitoringHistograms()
{
  //main histogram
  TH1D *bmh=new TH1D("btags",";b-tag multiplicity;Events",maxJets_+1,0,maxJets_+1);
  for(int ibin=1; ibin<=bmh->GetXaxis()->GetNbins(); ibin++)
    {
      TString label("="); label += (ibin-1); label+="b-tags";  
      if(ibin==bmh->GetXaxis()->GetNbins()) label.ReplaceAll("=","#geq");
      bmh->GetXaxis()->SetBinLabel(ibin,label);
    }
  controlHistos_.addHistogram( bmh );      
  controlHistos_.addHistogram( (TH1D *)bmh->Clone(TString("avg")+bmh->GetName() ) );
   
  //replicate for exclusive jet categories
  for(int ijets=2; ijets<=maxJets_; ijets++)
    {
      TString tag(""); tag += ijets; tag+= "jets";
      controlHistos_.initMonitorForStep(tag);
    }

  //control pseudo experiments
  controlHistos_.addHistogram( new TH1D("bias",";bias=#varepsilon_{b}-#bar{#varepsilon_{b}};Pseudo-experiments",100,-0.99,1.01) );
  controlHistos_.addHistogram( new TH1D("pull",";pull=(#varepsilon_{b}-#bar{#varepsilon_{b}}) / #sigma_{#varepsilon_{b}};Pseudo-experiments",100,-2.97,3.03) );
  controlHistos_.addHistogram( new TH1D("stat",";#sigma_{#varepsilon_{b}};Pseudo-experiments",100,0.0,1.0) );
}

//
void HFCMeasurement::resetHistograms()
{
  SelectionMonitor::StepMonitor_t &mons=controlHistos_.getAllMonitors();
  for(SelectionMonitor::StepMonitor_t::iterator it=mons.begin(); it!=mons.end(); it++)
    {
      for(SelectionMonitor::Monitor_t::iterator hit=it->second.begin(); hit!= it->second.end(); hit++)
	{
	  TString hname=hit->second->GetName();
	  if(hname.Contains("btags") && !hname.Contains("avg") ) hit->second->Reset("ICE");
	}
    }
}

//
void HFCMeasurement::saveMonitoringHistograms(TString tag)
{
  //open file
  TFile *fout=TFile::Open("HFCMeasurement.root","UPDATE");
  fout->cd();

  TDirectory *baseOutDir=fout->mkdir("localAnalysis/"+tag);
  SelectionMonitor::StepMonitor_t &mons=controlHistos_.getAllMonitors();
  for(SelectionMonitor::StepMonitor_t::iterator it=mons.begin(); it!=mons.end(); it++)
    {
      TDirectory *dir=baseOutDir->mkdir(it->first);
      dir->cd();
      for(SelectionMonitor::Monitor_t::iterator hit=it->second.begin(); hit!= it->second.end(); hit++)
	{
	  fixExtremities(hit->second,true,true);
	  
	  TString hname=hit->second->GetName();
	  if(hname.Contains("bias") || hname.Contains("pull") || hname.Contains("stat") )
	    {
	      if(nMeasurements_>0) hit->second->Scale(1./nMeasurements_);
	      hit->second->Fit("gaus","Q");
	    }
        }
    }
  
  //close file
  fout->Close();
  fout->Delete();
}

//
void HFCMeasurement::initHFCModel()
{
  //the observable
  model.bmult = new RooRealVar("bmult","N_{btags}",0.,float(maxJets_));
  model.bmult->setBins(maxJets_);  //if this is not set roofit will scan in very small steps
  
  //top decay modelling
  if(fitType_==FIT_R || fitType_==FIT_R_AND_EB || fitType_==FIT_R_AND_XSEC) model.r = new RooRealVar("r","R",smR_);
  else                                                                      model.r = new RooRealVar("r","R",1.0,0.0,2.0);
  model.lfacceptance = new RooRealVar("lfacceptance","A(R=0)",1.0);

  //b-tagging effiency modelling
  model.abseb = new RooRealVar("abseb","abs#varepsilon_{b}",effb[btagAlgo_]);
  if(fitType_==FIT_EB || fitType_==FIT_EB_AND_XSEC || fitType_==FIT_EB_AND_EQ) model.sfeb = new RooRealVar("sfeb","SF #varepsilon_{b}",sfb[btagAlgo_]);
  else                                                                        model.sfeb = new RooRealVar("sfeb","SF #varepsilon_{b}",sfb[btagAlgo_],0.0,2.0);
  model.sfeb_mean_constrain = new RooRealVar("meansfeb","#bar{SF #varepsilon_{b}}",sfb[btagAlgo_]);
  model.sfeb_sigma_constrain = new RooRealVar("uncsfeb","#sigma_{SF #varepsilon_{b}}",sfbUnc[btagAlgo_]);
  model.sfeb_constrain = new RooGaussian("sfeb_constrain","#varepsilon_{b} constrain",*model.sfeb,*model.sfeb_mean_constrain,*model.sfeb_sigma_constrain);
  model.eb = new RooFormulaVar("eb","@0*@1",RooArgSet(*model.abseb,*model.sfeb));

  //mistag efficiency modelling
  model.abseq = new RooRealVar("abseq","abs#varepsilon_{q}",effq[btagAlgo_]);
  if(fitType_==FIT_EB_AND_EQ) model.sfeq = new RooRealVar("sfeq","SF #varepsilon_{q}",sfq[btagAlgo_]);
  else                        model.sfeq = new RooRealVar("sfeq","SF #varepsilon_{q}",sfq[btagAlgo_],0.0,2.0);
  model.sfeq = new RooRealVar("sfeq","SF #varepsilon_{q}",sfq[btagAlgo_],0.0,2.0);
  model.sfeq_mean_constrain = new RooRealVar("meansfeq","#bar{SF #varepsilon_{q}}",sfq[btagAlgo_]);
  model.sfeq_sigma_constrain = new RooRealVar("uncsfeq","#sigma_{SF #varepsilon_{q}}",sfqUnc[btagAlgo_]);
  model.sfeq_constrain = new RooGaussian("sfeq_constrain","#varepsilon_{q} constrain",*model.sfeq,*model.sfeq_mean_constrain,*model.sfeq_sigma_constrain);
  model.eq = new RooFormulaVar("eq","@0*@1",RooArgSet(*model.abseq,*model.sfeq));

  //event type composition (alpha_i)
  model.alpha2=new RooRealVar("alpha2","#alpha_{2}",0.5,0.7);
  model.alpha2_mean_constrain = new RooRealVar("meanalpha2","#bar{#alpha}_{2}",alpha2[2]);
  model.alpha2_sigma_constrain = new RooRealVar("uncalpha2","#sigma_{#alpha_{2}}",alpha2Unc[0]);
  model.alpha2_constrain = new RooGaussian("ctr_alpha2","#alpha_{2} constrain",*model.alpha2,*model.alpha2_mean_constrain,*model.alpha2_sigma_constrain);
  
  model.alpha0 = new RooRealVar("alpha0","#alpha_{0}",0,0.2);
  model.alpha0_mean_constrain = new RooRealVar("meanalpha0","#bar{#alpha}_{0}",alpha0[2]);
  model.alpha0_sigma_constrain = new RooRealVar("uncalpha0","#sigma_{#alpha_{0}}",alpha0Unc[0]);
  model.alpha0_constrain = new RooGaussian("ctralpha0","#alpha_{0} constrain",*model.alpha0,*model.alpha0_mean_constrain,*model.alpha0_sigma_constrain);

  model.alpha1 = new RooFormulaVar("alpha1","1-@0-@1",RooArgSet(*model.alpha2,*model.alpha0));

  //add exclusive jet multiplicity models
  for(int jmult=2; jmult<=maxJets_; jmult++)
    {
      TString tag("_"); tag+=jmult;      
      
      //create the model
      HeavyFlavorPDF *mhfc = new HeavyFlavorPDF("hfcmodel"+tag,"hfcmodel"+tag,*model.bmult,*model.r,*model.eb,*model.eq,*model.alpha2,*model.alpha1,*model.alpha0);
      mhfc->setJetMultiplicity(jmult);
      model.pdfSet.add( *mhfc );
      
      //add the constrains
      RooProdPdf *modelconstr=0;
      if(fitType_==FIT_EB) 
	modelconstr=new RooProdPdf("modelconstr"+tag,"model x product of constrains",RooArgSet(*mhfc,*model.alpha2_constrain,*model.alpha0_constrain));
      else if(fitType_==FIT_R) 
	modelconstr=new RooProdPdf("modelconstr"+tag,"model x product of constrains",RooArgSet(*mhfc,*model.alpha2_constrain,*model.alpha0_constrain,*model.sfeq_constrain,*model.sfeb_constrain));

      model.constrPDFSet.add( *modelconstr );
    }
}

//
void HFCMeasurement::fitHFCtoEnsemble(EventSummaryHandler &evHandler, TString btagAlgo)
{
  if(evHandler.getEntries()==0) return;
  nMeasurements_++;
  resetHistograms();

  for(int i=0; i<evHandler.getEntries(); i++)
    {
      EventSummary_t &ev = evHandler.getEvent();
      if(eventCategory_!=0)
	{
	  if(eventCategory_==SFDileptons && ev.cat == dilepton::EMU) continue;
	  else if(eventCategory_==OFDileptons && ev.cat != dilepton::EMU) continue;
	  else if (eventCategory_ != ev.cat ) continue;
	}

      //count the number of b-tags
      PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      if(int(phys.jets.size())>maxJets_) continue;

      std::map<TString,int> nbtags;
      for(std::map<TString,float>::iterator it = algoCut.begin(); it!= algoCut.end(); it++) nbtags[it->first]=0;
      for(unsigned int ijet=0; ijet<phys.jets.size(); ijet++)
	{
	  nbtags["TCHEL"]  += (phys.jets[ijet].btag1>algoCut["TCHEL"]);
	  nbtags["TCHEM"]  += (phys.jets[ijet].btag1>algoCut["TCHEM"]);
	  nbtags["TCHPT"]  += (phys.jets[ijet].btag2>algoCut["TCHPT"]);
	  nbtags["SSVHEM"] += (phys.jets[ijet].btag3>algoCut["SSVHEM"]);
	  nbtags["JBPL"]   += (phys.jets[ijet].btag4>algoCut["JBPL"]);
	  nbtags["JBPM"]   += (phys.jets[ijet].btag4>algoCut["JBPM"]);
	  nbtags["JBPT"]   += (phys.jets[ijet].btag4>algoCut["JBPT"]);
	}

      std::vector<TString> catsToFill;
      catsToFill.push_back("all");
      TString tag(""); tag += phys.jets.size(); tag+= "jets";
      catsToFill.push_back(tag);
      for(std::vector<TString>::iterator it = catsToFill.begin(); it!= catsToFill.end(); it++)
	controlHistos_.fillHisto("btags",*it,nbtags[btagAlgo]);
    }
 

  //run the fit
  runHFCFit();
}


//
void HFCMeasurement::runHFCFit()
{

  RooArgSet allLL,jetOccupancies;  
  RooDataHist* alldata = new RooDataHist("data","data", RooArgList(*model.bmult));
  TIterator *pdfIt     = model.constrPDFSet.createIterator();
  for(int jmult=2; jmult<=maxJets_; jmult++)
    {
      TString tag(""); tag += jmult; tag+= "jets";
      TH1 *h = controlHistos_.getHisto("btags",tag);

      RooDataHist* ds = new RooDataHist("data"+tag,"data"+tag, RooArgList(*model.bmult), h);
      alldata->add(*ds);

      RooRealVar *jetocc = new RooRealVar("jetocc"+tag,"occ"+tag,h->Integral());
      jetOccupancies.add(*jetocc);

      RooProdPdf *modelconstr = (RooProdPdf *) pdfIt->Next();
      RooAbsReal *nll=0;
      if(fitType_==FIT_EB)      nll = modelconstr->createNLL(*ds,Constrain(RooArgSet(*model.alpha2,*model.alpha0,*model.sfeq_constrain)));
      else if (fitType_==FIT_R) nll = modelconstr->createNLL(*ds,Constrain(RooArgSet(*model.alpha2,*model.alpha0,*model.sfeq_constrain,*model.sfeb_constrain)));
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
  if(false)
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
	  frame = model.sfeb->frame(Title("Likelihood"),Range(0.,1.));
	  frame->GetXaxis()->SetTitle("SF #varepsilon_{b}");
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
      if(fitType_==FIT_EB)  sprintf(buf,"SF #varepsilon_{b}=%3.3f #pm %3.2f",model.sfeb->getVal(),model.sfeb->getError());
      else                  sprintf(buf,"R=%3.2f #pm %3.2f",model.r->getVal(),model.r->getError());
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
  float bias=(fitType_==0 ? model.sfeb->getVal()-1.0 : model.r->getVal()-smR_);
  float unc =(fitType_==0 ? model.sfeb->getError()   :  model.r->getError());
  if(unc>0.01)
    {
      controlHistos_.fillHisto("bias","all",bias);
      controlHistos_.fillHisto("stat","all",unc);
      controlHistos_.fillHisto("pull","all",bias/unc);
 
      // Access basic information
      cout << "EDM = " << r->edm() << endl
	   << "-log(L) at minimum = " << r->minNll() << endl   
	   << "final value of floating parameters" << endl ;
      r->floatParsFinal().Print("s") ;
    }
  cout << "Somewhere" << endl;
}
