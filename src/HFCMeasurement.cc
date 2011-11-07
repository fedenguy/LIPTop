#include "LIP/Top/interface/HFCMeasurement.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "RooNumIntConfig.h"
#include "RooNLLVar.h"

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
  controlHistos_.addHistogram( (TH1D *)bmh->Clone("avgbtags") );
   
  //replicate for exclusive jet categories
  for(int ijets=2; ijets<=maxJets_; ijets++)
    {
      TString tag(""); tag += ijets; tag+= "jets";
      controlHistos_.addHistogram( (TH1D *)bmh->Clone( "btags_" + tag ) );
      controlHistos_.addHistogram( (TH1D *)bmh->Clone( "avgbtags_" + tag ) );
    }

  //replicate for different categories in 2 jet bin
  for(int icat=0; icat<MAXCATEGORIES*(MAXCATEGORIES-1); icat++)
    {
      TString tag(""); tag += icat;
      controlHistos_.addHistogram( (TH1D *)bmh->Clone( "btags_" + tag ) );
      controlHistos_.addHistogram( (TH1D *)bmh->Clone( "avgbtags_" + tag ) );
    }

  //control pseudo experiments
  controlHistos_.addHistogram( new TH1D("bias",";bias=#varepsilon_{b}-#bar{#varepsilon_{b}};Pseudo-experiments",100,-0.99,1.01) );
  controlHistos_.addHistogram( new TH1D("pull",";pull=(#varepsilon_{b}-#bar{#varepsilon_{b}}) / #sigma_{#varepsilon_{b}};Pseudo-experiments",100,-2.97,3.03) );
  controlHistos_.addHistogram( new TH1D("stat",";#sigma_{#varepsilon_{b}};Pseudo-experiments",100,0.0,1.0) );

  //instantiate for different categories
  controlHistos_.initMonitorForStep("ee");
  controlHistos_.initMonitorForStep("emu");
  controlHistos_.initMonitorForStep("mumu");
  controlHistos_.initMonitorForStep("ll");
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
  if(isInit_) 
    {
      //cout << "[HFCMeasurement::initHFCModel] model has already been initiated - nothing done" << endl;
      return;
    }

  //the observable
  model.bmult = new RooRealVar("bmult","N_{btags}",0.,float(maxJets_));
  //model.bmult->setBins(maxJets_+1);
  RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("maxSteps",30);

  //top decay modelling
  if(fitType_==FIT_R || fitType_==FIT_R_AND_XSEC) model.r = new RooRealVar("r","R",1.0,0.85,1.25);
  else if(fitType_==FIT_R_AND_EB)                 model.r = new RooRealVar("r","R",1.0,0.95,1.05);
  else                                            model.r = new RooRealVar("r","R",smR_);
  model.lfacceptance = new RooRealVar("lfacceptance","A(R=0)",1.0);

  //
  //FIX-ME do this per jet multiplicity bin
  //
  //b-tagging effiency modelling
  model.abseb = new RooRealVar("abseb","abs#varepsilon_{b}",effb_[0]);
  if(fitType_==FIT_EB || fitType_==FIT_R_AND_EB || fitType_==FIT_EB_AND_XSEC || fitType_==FIT_EB_AND_EQ) 
    model.sfeb = new RooRealVar("sfeb","SF #varepsilon_{b}",sfb_[0],0.7,min(1./effb_[0],1.3));
  else
    model.sfeb = new RooRealVar("sfeb","SF #varepsilon_{b}",sfb_[0]);
  model.sfeb_mean_constrain = new RooRealVar("meansfeb","#bar{SF #varepsilon_{b}}",sfb_[0]);
  model.sfeb_sigma_constrain = new RooRealVar("uncsfeb","#sigma_{SF #varepsilon_{b}}",sfbUnc_[0]);
  model.sfeb_constrain = new RooGaussian("sfeb_constrain","#varepsilon_{b} constrain",*model.sfeb,*model.sfeb_mean_constrain,*model.sfeb_sigma_constrain);
  model.eb = new RooFormulaVar("eb","@0*@1",RooArgSet(*model.abseb,*model.sfeb));

  for(int icat=0; icat<MAXCATEGORIES; icat++)
    {
      TString icatLabel(""); icatLabel +=icat;
      model.diff_sfeb[icat] = new RooRealVar("sfeb_"+icatLabel,"SF #varepsilon_{b}^{"+icatLabel+"}",sfb_[0],0.7,min(1./effb_[0],1.3));
      model.diff_eb[icat]   = new RooFormulaVar("eb_"+icatLabel,"@0*@1",RooArgSet(*model.abseb,*model.diff_sfeb[icat]));
    }

  //mistag efficiency modelling
  model.abseq = new RooRealVar("abseq","abs#varepsilon_{q}",effq_[0]);
  if(fitType_==FIT_EB_AND_EQ) model.sfeq = new RooRealVar("sfeq","SF #varepsilon_{q}",sfq_[0],0.5,min(1./effq_[0],1.5));
  else                        model.sfeq = new RooRealVar("sfeq","SF #varepsilon_{q}",sfq_[0]);
  model.sfeq_mean_constrain = new RooRealVar("meansfeq","#bar{SF #varepsilon_{q}}",sfq_[0]);
  model.sfeq_sigma_constrain = new RooRealVar("uncsfeq","#sigma_{SF #varepsilon_{q}}",sfqUnc_[0]);
  model.sfeq_constrain = new RooGaussian("sfeq_constrain","#varepsilon_{q} constrain",*model.sfeq,*model.sfeq_mean_constrain,*model.sfeq_sigma_constrain);
  model.eq = new RooFormulaVar("eq","@0*@1",RooArgSet(*model.abseq,*model.sfeq));

  //add exclusive jet multiplicity models
  for(int jmult=2; jmult<=maxJets_; jmult++)
    {
      TString tag(""); tag+=jmult;      

      model.jetocc[jmult-2] = new RooRealVar("jetocc"+tag,"occ"+tag,1);
      
      Float_t val(0),valerr(0);
      if(fcorrect_.find(jmult) != fcorrect_.end() ) { val = fcorrect_[jmult];  valerr = fcorrectUnc_[jmult]; }
      else                                          { val = fcorrect_[0];      valerr = fcorrectUnc_[0]; }
      model.fcorrect[jmult-2]=new RooRealVar("fcorrect_"+tag,"f_{correct}^{"+tag+"}",val);//,0.0,0.5);
      //       model.fcorrect_mean_constrain[jmult-2] = new RooRealVar("meanfcorrect_"+tag,"#bar{f_{correct}^{"+tag+"}}",val);
      //       model.fcorrect_sigma_constrain[jmult-2] = new RooRealVar("uncfcorrect_"+tag,"#sigma_{f_{correct}^{"+tag+"}}",valerr);
      //       model.fcorrect_constrain[jmult-2] = new RooGaussian("ctr_fcorrect_"+tag,"f_{correct}^{"+tag+"} constrain",*model.fcorrect[jmult-2],*model.fcorrect_mean_constrain[jmult-2],*model.fcorrect_sigma_constrain[jmult-2]);
      
      if(fttbar_.find(jmult) != fttbar_.end() ) { val = fttbar_[jmult];  valerr = fttbarUnc_[jmult]; }
      else                                      { val = fttbar_[0];      valerr = fttbarUnc_[0]; }
      model.fttbar[jmult-2]=new RooRealVar("fttbar_"+tag,"f_{t#bar{t}}^{"+tag+"}",val,val-3*valerr,val+3*valerr);//,0.0,0.5);
      model.fttbar_mean_constrain[jmult-2] = new RooRealVar("meanfttbar_"+tag,"#bar{f_{t#bar{t}}^{"+tag+"}}",val);
      model.fttbar_sigma_constrain[jmult-2] = new RooRealVar("uncfttbar_"+tag,"#sigma_{f_{t#bar{t}}^{"+tag+"}}",valerr);
      model.fttbar_constrain[jmult-2] = new RooGaussian("ctr_fttbar_"+tag,"f_{t#bar{t}}^{"+tag+"} constrain",*model.fttbar[jmult-2],*model.fttbar_mean_constrain[jmult-2],*model.fttbar_sigma_constrain[jmult-2]);
      
      if(fsingletop_.find(jmult) != fsingletop_.end() ) { val = fsingletop_[jmult];  valerr = fsingletopUnc_[jmult]; }
      else                                              { val = fsingletop_[0];      valerr = fsingletopUnc_[0]; }
      model.fsingletop[jmult-2]=new RooRealVar("fsingletop_"+tag,"f_{t}^{"+tag+"}",val,val-3*valerr,val+3*valerr);//,0.0,0.5);
      model.fsingletop_mean_constrain[jmult-2] = new RooRealVar("meanfsingletop_"+tag,"#bar{f_{t}^{"+tag+"}}",val);
      model.fsingletop_sigma_constrain[jmult-2] = new RooRealVar("uncfsingletop_"+tag,"#sigma_{f_{t}^{"+tag+"}}",valerr);
      model.fsingletop_constrain[jmult-2] = new RooGaussian("ctr_fsingletop_"+tag,"f_{t}^{"+tag+"} constrain",*model.fsingletop[jmult-2],*model.fsingletop_mean_constrain[jmult-2],*model.fsingletop_sigma_constrain[jmult-2]);

      //create the model
      HeavyFlavorPDF *mhfc = new HeavyFlavorPDF("hfcmodel_"+tag,"hfcmodel_"+tag,*model.bmult,*model.r,*model.eb,*model.eq,*model.fcorrect[jmult-2],*model.fttbar[jmult-2],*model.fsingletop[jmult-2],*model.jetocc[jmult-2]);
      mhfc->setJetMultiplicity(jmult);
      model.pdfSet.add( *mhfc );
      
      //add the constrains
      RooProdPdf *modelconstr=0;
      if(fitType_==FIT_EB || fitType_==FIT_R || fitType_==FIT_R_AND_EB) 
	{
	  //modelconstr=new RooProdPdf("modelconstr_"+tag,"model x product of constrains",RooArgSet(*mhfc/*,*model.fttbar_constrain[jmult-2],*model.fsingletop_constrain[jmult-2]*/));
	  modelconstr=new RooProdPdf("modelconstr_"+tag,"model x product of constrains",RooArgSet(*mhfc,*model.fttbar_constrain[jmult-2],*model.fsingletop_constrain[jmult-2]));
	}
      else if(fitType_==FIT_EB_AND_EQ /*|| fitType_==FIT_R_AND_EQ*/)
	modelconstr=new RooProdPdf("modelconstr_"+tag,"model x product of constrains",RooArgSet(*mhfc/*,*model.fttbar_constrain[jmult-2],*model.fsingletop_constrain[jmult-2]*/,*model.sfeq_constrain,*model.sfeb_constrain));
      
      model.constrPDFSet.add( *modelconstr );
    }

  isInit_=true;
}


//
void HFCMeasurement::resetModelValues()
{
  if(!isInit_) return;
  model.r->setVal(smR_);
  model.abseb->setVal(effb_[0]);
  model.sfeb->setVal(sfb_[0]);
  model.abseq->setVal(effq_[0]);
  model.sfeq->setVal(sfq_[0]);
}


//
void HFCMeasurement::fitHFCtoEnsemble(top::EventSummaryHandler &evHandler, TString dilCategory)
{
  if(evHandler.getEntries()==0) return;

  //restart all over again
  initHFCModel();
  resetModelValues();
  resetHistograms();

  nMeasurements_++;
  for(int i=0; i<evHandler.getEntries(); i++)
    {
      evHandler.getEntry(i);
      top::EventSummary_t &ev = evHandler.getEvent();
      
      //filter the even type
      if(dilCategory=="ll" && ev.cat == EMU)         continue;
      else if(dilCategory=="emu" && ev.cat != EMU)   continue;
      else if(dilCategory=="ee" && ev.cat != EE)     continue;
      else if(dilCategory=="mumu" && ev.cat != MUMU) continue;
      else if(dilCategory!="all") continue;

      //count the number of b-tags
      top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      
      int njets=0;
      int nbtags=0;
      for(unsigned int ijet=0; ijet<phys.jets.size(); ijet++)
	{
	  //	  if(phys.jets[ijet].pt()<120 || phys.jets[ijet].pt()>10000) continue;
	  njets++;
	  double btag(-9999.);
	  if(btagAlgo_.Contains("TCHE") )  btag = phys.jets[ijet].btag1;
	  else if(btagAlgo_.Contains("TCHP") )  btag = phys.jets[ijet].btag2;
	  else if(btagAlgo_.Contains("SSVHE") ) btag = phys.jets[ijet].btag3;
	  else if(btagAlgo_.Contains("JBP") )   btag = phys.jets[ijet].btag4;
	  else if(btagAlgo_.Contains("JP") )    btag = phys.jets[ijet].btag5;
	  else if(btagAlgo_.Contains("SSVHP") )  btag = phys.jets[ijet].btag6;
	  else if(btagAlgo_.Contains("CSV") )  btag = phys.jets[ijet].btag7;
	  nbtags += (btag>algoCut_);
	}
      if(njets>maxJets_) continue;

      TString tag("btags_"); tag += njets; tag+= "jets";
      controlHistos_.fillHisto("btags",dilCategory,nbtags);
      controlHistos_.fillHisto(tag,dilCategory,nbtags);

      if(njets>2) continue;
      float pt1=max(phys.jets[0].pt(),phys.jets[1].pt());
      float pt2=min(phys.jets[0].pt(),phys.jets[1].pt());

      int icat(0);
      if     (pt2<40 && pt1<40)   icat=0;
      else if(pt2<40 && pt1<50)   icat=1;
      else if(pt2<40 && pt1<60)   icat=2;
      else if(pt2<40 && pt1<80)   icat=3;
      else if(pt2<40 && pt1<120)  icat=4;
      else if(pt2<40 && pt1>=120) icat=5;

      else if(pt2<50 && pt1<50)   icat=6;
      else if(pt2<50 && pt1<60)   icat=7;
      else if(pt2<50 && pt1<80)   icat=8;
      else if(pt2<50 && pt1<120)  icat=9;
      else if(pt2<50 && pt1>=120) icat=10;

      else if(pt2<60 && pt1<60)   icat=11;
      else if(pt2<60 && pt1<80)   icat=12;
      else if(pt2<60 && pt1<120)  icat=13;
      else if(pt2<60 && pt1>=120) icat=14;

      else if(pt2<80 && pt1<80)   icat=15;
      else if(pt2<80 && pt1<120)  icat=16;
      else if(pt2<80 && pt1>=120) icat=17;

      else if(pt2<120 && pt1<120)  icat=18;
      else if(pt2<120 && pt1>=120) icat=19;

      else                         icat=20;

      tag=""; tag += icat;
      controlHistos_.fillHisto( "btags_" + tag, dilCategory, nbtags );
      controlHistos_.fillHisto( "avgbtags_" + tag, dilCategory, nbtags );      
    }
  
  //run the fit
  runHFCFit(dilCategory);
  runHFCDiffFit(dilCategory);
}


//
void HFCMeasurement::runHFCFit(TString dilCategory)
{

  //all data
  RooDataHist* alldata = new RooDataHist("data","data", RooArgList(*model.bmult));

  //build the likelihoods per category
  RooArgSet allLL, jetOccupancies;  
  TIterator *pdfIt     = model.constrPDFSet.createIterator();
  TIterator *modelpdfIt     = model.pdfSet.createIterator();
  std::vector<RooNLLVar *> theLLsPerCategory;
  for(int jmult=2; jmult<=maxJets_; jmult++)
    {
      TString tag(""); tag += jmult; tag+= "jets";
      TH1 *h = controlHistos_.getHisto("btags_"+tag,dilCategory);

      //convert to a data hist
      RooDataHist* ds = new RooDataHist("data"+tag,"data"+tag, RooArgList(*model.bmult), h);
      alldata->add(*ds);

      //build the model for this category
      model.jetocc[jmult-2]->setVal(h->Integral());
      jetOccupancies.add(*(model.jetocc[jmult-2]));

      //HeavyFlavorPDF *hfcmodel = dynamic_cast<HeavyFlavorPDF *>( modelpdfIt->Next() );
      //cout << jmult << " " << h->Integral() << endl;
      //hfcmodel->setEventMultiplicity(h->Integral());

      RooNLLVar *nll=0;
      RooProdPdf *modelconstr = dynamic_cast<RooProdPdf *>( pdfIt->Next() );
      if(fitType_==FIT_EB || fitType_==FIT_R || fitType_==FIT_R_AND_EB) 
	{
	  //nll = (RooNLLVar *) modelconstr->createNLL(*ds);
	  nll = (RooNLLVar *) modelconstr->createNLL(*ds,Constrain(RooArgSet(*model.fttbar[jmult-2],*model.fsingletop[jmult-2])));
	}
      else if(fitType_==FIT_EB_AND_EQ /*|| fitType_==FIT_R_AND_EQ*/)    nll = (RooNLLVar *) modelconstr->createNLL(*ds,Constrain(RooArgSet(/**model.fttbar[jmult-2],*model.fsingletop[jmult-2],*/*model.sfeq)));
      TString itit("=");    itit+=jmult; itit += " jets";
      nll->SetTitle(itit);
      
      allLL.add(*nll); 
      theLLsPerCategory.push_back(nll);
    }

  //add up all the likelihoods
  RooAddition *combll = new RooAddition("combll","combll",allLL);      

  //maximize the likelihood (reinforce the error should be taken from +/- 0.5 contour
  RooMinuit minuit(*combll); 
  minuit.migrad();
  minuit.setErrorLevel(0.5);
  minuit.hesse();
  RooFitResult *r=minuit.save();

  //draw the result
  if(true)
    {
      setStyle();

      TCanvas *c=getNewCanvas("combination","combination",false);
      c->cd();
      c->SetWindowSize(600,600);
      c->SetGridx();
      c->SetGridy();     

      //likelihood (main frame)
      RooPlot *frame = 0;
      TString label("CMS preliminary, #sqrt{s}=7 TeV \\");
      if(fitType_==FIT_EB || fitType_==FIT_EB_AND_EQ) 
	{
	  frame = model.sfeb->frame(Title("Likelihood"),Range(0.8,1.1));
	  frame->GetXaxis()->SetTitle("SF #varepsilon_{b}");
	  char buf[100];
	  sprintf(buf,"SF #varepsilon_{b}=%3.3f #pm %3.2f\\",model.sfeb->getVal(),model.sfeb->getError());
	  label += buf;
	  char buf2[100];
	  sprintf(buf2,"#varepsilon_{b}=%3.3f #pm %3.2f",model.abseb->getVal()*model.sfeb->getVal(),model.abseb->getVal()*model.sfeb->getError());
	  label +=buf2;
	  cout << buf << endl << buf2 << endl;
	}
      else if (fitType_==FIT_R || fitType_==FIT_R_AND_EB) 
	{
	  frame = model.r->frame(Title("Likelihood"),Range(0.7,1.4)) ;    
	  frame->GetXaxis()->SetTitle("R");
	  char buf[100];
	  sprintf(buf,"R=%3.2f #pm %3.2f",model.r->getVal(),model.r->getError());
	  label += buf;
	}
      combll->plotOn(frame,ShiftToZero(),Name("ll"));
      
      for(int jmult=2; jmult<=maxJets_; jmult++)
	{
	  TString itit("=");    itit+=jmult; itit += " jets";
	  TString iname("ll_"); iname +=jmult;
	  theLLsPerCategory[jmult-2]->plotOn(frame,ShiftToZero(),
					     LineColor(kGreen+4-2*(jmult-2)),
					     LineStyle(kDashed),
					     LineWidth(1),
					     Name(iname),
					     MoveToBack());
	}
      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetYaxis()->SetTitle("-Log(L/L_{Max})");
      frame->GetYaxis()->SetTitleOffset(1);
      frame->Draw();
      
      //fit result
      formatForCmsPublic(c,0,label,1);

      //the model fit
      TPad *npad = new TPad("llpad","ll", 0.6, 0.6, 0.9, 0.9);
      npad->Draw();
      npad->cd();
      frame = model.bmult->frame();

      RooAddPdf incmodel("incmodel","Inclusive Model",model.pdfSet,jetOccupancies);
      alldata->plotOn(frame,Binning(maxJets_),DrawOption("pz"));
      incmodel.plotOn(frame,Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
      frame->Draw();

      if(fitType_ == FIT_EB_AND_EQ || fitType_==FIT_R_AND_EB)
	{
	  c = new TCanvas("contour","contour");
	  c->cd();
	  c->SetWindowSize(600,600);
	  c->SetGridx();
	  c->SetGridy();     
	  RooPlot *plot=0;
	  if(fitType_==FIT_R_AND_EB)
	    {
	      plot = minuit.contour(*model.r,*model.sfeb,1,2,3) ;
              plot->SetTitle("Contour for 1s,2s,3s between r and eb") ;
	    }
	  else
	    {
	      plot = minuit.contour(*model.sfeb,*model.sfeq,1,2,3) ;
	      plot->SetTitle("Contour for 1s,2s,3s between eb and eq") ;
	    }	 
	  plot->Draw();
	  formatForCmsPublic(c,0,label,3);
	}
    }
  //  float bias=(fitType_==0 ? model.sfeb->getVal()-1.0 : model.r->getVal()-smR_);
  //  float unc =(fitType_==0 ? model.sfeb->getError()   :  model.r->getError());
  //  if(unc>0.01)
  //    {
  //      controlHistos_.fillHisto("bias","all",bias);
  //      controlHistos_.fillHisto("stat","all",unc);
  //      controlHistos_.fillHisto("pull","all",bias/unc);
  
  // Access basic information
  cout << "EDM = " << r->edm() << endl
       << "-log(L) at minimum = " << r->minNll() << endl   
       << "final value of floating parameters" << endl ;
  r->floatParsFinal().Print("s") ;
}


//
void HFCMeasurement::runHFCDiffFit(TString dilCategory)
{

  //all data
  RooDataHist* alldata = new RooDataHist("data","data", RooArgList(*model.bmult));

  //build the likelihoods per category
  RooArgSet allLL, jetOccupancies;  
  TIterator *pdfIt       = model.constrPDFSet.createIterator();
  TIterator *modelpdfIt  = model.pdfSet.createIterator();
  std::vector<RooNLLVar *> theLLsPerCategory;
  for(int jmult=2; jmult<=maxJets_; jmult++)
    {
      TString tag(""); tag += jmult; tag+= "jets";
      TH1 *h = controlHistos_.getHisto("btags_"+tag,dilCategory);

      //convert to a data hist
      RooDataHist* ds = new RooDataHist("data"+tag,"data"+tag, RooArgList(*model.bmult), h);
      alldata->add(*ds);

      //build the model for this category
      model.jetocc[jmult-2]->setVal(h->Integral());
      jetOccupancies.add(*(model.jetocc[jmult-2]));

      //HeavyFlavorPDF *hfcmodel = dynamic_cast<HeavyFlavorPDF *>( modelpdfIt->Next() );
      //cout << jmult << " " << h->Integral() << endl;
      //hfcmodel->setEventMultiplicity(h->Integral());

      RooNLLVar *nll=0;
      RooProdPdf *modelconstr = dynamic_cast<RooProdPdf *>( pdfIt->Next() );
      if(fitType_==FIT_EB || fitType_==FIT_R || fitType_==FIT_R_AND_EB) nll = (RooNLLVar *) modelconstr->createNLL(*ds);//,Constrain(RooArgSet(*model.fttbar[jmult-2],*model.fsingletop[jmult-2])));
      else if(fitType_==FIT_EB_AND_EQ /*|| fitType_==FIT_R_AND_EQ*/)    nll = (RooNLLVar *) modelconstr->createNLL(*ds,Constrain(RooArgSet(/**model.fttbar[jmult-2],*model.fsingletop[jmult-2],*/*model.sfeq)));
      TString itit("=");    itit+=jmult; itit += " jets";
      nll->SetTitle(itit);
      
      allLL.add(*nll); 
      theLLsPerCategory.push_back(nll);
    }

  //add up all the likelihoods
  RooAddition *combll = new RooAddition("combll","combll",allLL);      

  //maximize the likelihood (reinforce the error should be taken from +/- 0.5 contour
  RooMinuit minuit(*combll); 
  minuit.migrad();
  minuit.setErrorLevel(0.5);
  minuit.hesse();
  RooFitResult *r=minuit.save();

  //draw the result
  if(true)
    {
      setStyle();

      TCanvas *c=getNewCanvas("combination","combination",false);
      c->cd();
      c->SetWindowSize(600,600);
      c->SetGridx();
      c->SetGridy();     

      //likelihood (main frame)
      RooPlot *frame = 0;
      TString label("CMS preliminary, #sqrt{s}=7 TeV \\");
      if(fitType_==FIT_EB || fitType_==FIT_EB_AND_EQ) 
	{
	  frame = model.sfeb->frame(Title("Likelihood"),Range(0.8,1.1));
	  frame->GetXaxis()->SetTitle("SF #varepsilon_{b}");
	  char buf[100];
	  sprintf(buf,"SF #varepsilon_{b}=%3.3f #pm %3.2f\\",model.sfeb->getVal(),model.sfeb->getError());
	  label += buf;
	  char buf2[100];
	  sprintf(buf2,"#varepsilon_{b}=%3.3f #pm %3.2f",model.abseb->getVal()*model.sfeb->getVal(),model.abseb->getVal()*model.sfeb->getError());
	  label +=buf2;
	}
      else if (fitType_==FIT_R || fitType_==FIT_R_AND_EB) 
	{
	  frame = model.r->frame(Title("Likelihood"),Range(0.7,1.4)) ;    
	  frame->GetXaxis()->SetTitle("R");
	  char buf[100];
	  sprintf(buf,"R=%3.2f #pm %3.2f",model.r->getVal(),model.r->getError());
	  label += buf;
	}
      combll->plotOn(frame,ShiftToZero(),Name("ll"));
      
      for(int jmult=2; jmult<=maxJets_; jmult++)
	{
	  TString itit("=");    itit+=jmult; itit += " jets";
	  TString iname("ll_"); iname +=jmult;
	  theLLsPerCategory[jmult-2]->plotOn(frame,ShiftToZero(),
					     LineColor(kGreen+4-2*(jmult-2)),
					     LineStyle(kDashed),
					     LineWidth(1),
					     Name(iname),
					     MoveToBack());
	}
      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetYaxis()->SetTitle("-Log(L/L_{Max})");
      frame->GetYaxis()->SetTitleOffset(1);
      frame->Draw();
      
      //fit result
      formatForCmsPublic(c,0,label,1);

      //the model fit
      TPad *npad = new TPad("llpad","ll", 0.6, 0.6, 0.9, 0.9);
      npad->Draw();
      npad->cd();
      frame = model.bmult->frame();

      RooAddPdf incmodel("incmodel","Inclusive Model",model.pdfSet,jetOccupancies);
      alldata->plotOn(frame,Binning(maxJets_),DrawOption("pz"));
      incmodel.plotOn(frame,Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
      frame->Draw();

      if(fitType_ == FIT_EB_AND_EQ || fitType_==FIT_R_AND_EB)
	{
	  c = new TCanvas("contour","contour");
	  c->cd();
	  c->SetWindowSize(600,600);
	  c->SetGridx();
	  c->SetGridy();     
	  RooPlot *plot=0;
	  if(fitType_==FIT_R_AND_EB)
	    {
	      plot = minuit.contour(*model.r,*model.sfeb,1,2,3) ;
              plot->SetTitle("Contour for 1s,2s,3s between r and eb") ;
	    }
	  else
	    {
	      plot = minuit.contour(*model.sfeb,*model.sfeq,1,2,3) ;
	      plot->SetTitle("Contour for 1s,2s,3s between eb and eq") ;
	    }	 
	  plot->Draw();
	  formatForCmsPublic(c,0,label,3);
	}
    }
  //  float bias=(fitType_==0 ? model.sfeb->getVal()-1.0 : model.r->getVal()-smR_);
  //  float unc =(fitType_==0 ? model.sfeb->getError()   :  model.r->getError());
  //  if(unc>0.01)
  //    {
  //      controlHistos_.fillHisto("bias","all",bias);
  //      controlHistos_.fillHisto("stat","all",unc);
  //      controlHistos_.fillHisto("pull","all",bias/unc);
  
  // Access basic information
  cout << "EDM = " << r->edm() << endl
       << "-log(L) at minimum = " << r->minNll() << endl   
       << "final value of floating parameters" << endl ;
  r->floatParsFinal().Print("s") ;
}
