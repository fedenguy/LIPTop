#include "LIP/Top/interface/HFCMeasurement.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "RooNumIntConfig.h"
#include "RooNLLVar.h"

using namespace RooFit;
using namespace std;

//
void HFCMeasurement::bookMonitoringHistograms()
{
  setStyle();

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
      TString tag(""); tag += ijets;
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
  if(isInit_) return;

  //b-tag multiplicity is the observable
  //model.bmult = new RooRealVar("bmult","N_{btags}",0.,float(maxJets_));
  model.bmult = new RooRealVar("bmult","N_{btags}",0.,float(maxJets_+1));
  RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("maxSteps",50);

  //what do we want to fit from the observable
  bool fitR(fitType_==FIT_R || fitType_==FIT_R_AND_XSEC || fitType_==FIT_R_AND_EB);
  bool fitEb(fitType_==FIT_EB || fitType_==FIT_R_AND_EB || fitType_==FIT_EB_AND_XSEC || fitType_==FIT_EB_AND_EQ);
  bool fitEq(fitType_==FIT_EB_AND_EQ);

  //top decay modelling
  if(!fitR)                        model.r = new RooRealVar("r","R",smR_);
  else if(fitType_==FIT_R_AND_EB)  model.r = new RooRealVar("r","R",1.0,0.9,1.1);
  else                             model.r = new RooRealVar("r","R",1.0,0.,2.0);



  //this is a correction for acceptance (keep constant for now)
  model.lfacceptance = new RooRealVar("lfacceptance","A(R=0)",1.0);

  //b-tagging effiency modelling
  model.abseb = new RooRealVar("abseb","abs#varepsilon_{b}",effb_);
  if(fitEb)    model.sfeb = new RooRealVar("sfeb","SF #varepsilon_{b}",sfb_,0.7,min(1./effb_,1.3));
  else         model.sfeb = new RooRealVar("sfeb","SF #varepsilon_{b}",sfb_);
  model.sfeb_mean_constrain = new RooRealVar("meansfeb","#bar{SF #varepsilon_{b}}",sfb_);
  model.sfeb_sigma_constrain = new RooRealVar("uncsfeb","#sigma_{SF #varepsilon_{b}}",sfbUnc_);
  model.sfeb_constrain = new RooGaussian("sfeb_constrain","#varepsilon_{b} constrain",*model.sfeb,*model.sfeb_mean_constrain,*model.sfeb_sigma_constrain);
  model.eb = new RooFormulaVar("eb","@0*@1",RooArgSet(*model.abseb,*model.sfeb));
  
  //mistag efficiency modelling
  model.abseq = new RooRealVar("abseq","abs#varepsilon_{q}",effq_);
  if(fitEq) model.sfeq = new RooRealVar("sfeq","SF #varepsilon_{q}",sfq_,0.5,min(1./effq_,1.5));
  else      model.sfeq = new RooRealVar("sfeq","SF #varepsilon_{q}",sfq_);
  model.sfeq_mean_constrain  = new RooRealVar("meansfeq","#bar{SF #varepsilon_{q}}",sfq_);
  model.sfeq_sigma_constrain = new RooRealVar("uncsfeq","#sigma_{SF #varepsilon_{q}}",sfqUnc_);
  model.sfeq_constrain       = new RooGaussian("sfeq_constrain","#varepsilon_{q} constrain",*model.sfeq,*model.sfeq_mean_constrain,*model.sfeq_sigma_constrain);
  model.eq                   = new RooFormulaVar("eq","@0*@1",RooArgSet(*model.abseq,*model.sfeq));

  //add exclusive jet multiplicity models
  int icat(0);
  for(std::set<TString>::iterator cIt=categoryKeys_.begin(); cIt!=categoryKeys_.end(); cIt++,icat++)
    {
      TString tag=*cIt;
      model.jetocc[icat] = new RooRealVar("jetocc"+tag,"occ"+tag,1);
      
      if(fitEq || fitType_==FIT_R_AND_EB)  model.fcorrect[icat] = new RooRealVar("fcorrect_"+tag,"f_{correct}^{"+tag+"}",fcorrect_[tag]);
      else                                 model.fcorrect[icat] = new RooRealVar("fcorrect_"+tag,"f_{correct}^{"+tag+"}",fcorrect_[tag],fcorrect_[tag]-3*fcorrectUnc_[tag],fcorrect_[tag]+3*fcorrectUnc_[tag]);
      model.fcorrect_mean_constrain[icat]  = new RooRealVar("meanfcorrect_"+tag,"#bar{f_{correct}^{"+tag+"}}",fcorrect_[tag]);
      model.fcorrect_sigma_constrain[icat] = new RooRealVar("uncfcorrect_"+tag,"#sigma_{f_{correct}^{"+tag+"}}",fcorrectUnc_[tag]);
      model.fcorrect_constrain[icat]       = new RooGaussian("ctr_fcorrect_"+tag,"f_{correct}^{"+tag+"} constrain",*model.fcorrect[icat],*model.fcorrect_mean_constrain[icat],*model.fcorrect_sigma_constrain[icat]);
	  	  
      if(fitEq ||  fitType_==FIT_R_AND_EB)  model.fttbar[icat] = new RooRealVar("fttbar_"+tag,"f_{t#bar{t}}^{"+tag+"}",fttbar_[tag]);
      else                                  model.fttbar[icat] = new RooRealVar("fttbar_"+tag,"f_{t#bar{t}}^{"+tag+"}",fttbar_[tag],fttbar_[tag]-3*fttbarUnc_[tag],fttbar_[tag]+3*fttbarUnc_[tag]);
      model.fttbar_mean_constrain[icat]    = new RooRealVar("meanfttbar_"+tag,"#bar{f_{t#bar{t}}^{"+tag+"}}",fttbar_[tag]);
      model.fttbar_sigma_constrain[icat]   = new RooRealVar("uncfttbar_"+tag,"#sigma_{f_{t#bar{t}}^{"+tag+"}}",fttbarUnc_[tag]);
      model.fttbar_constrain[icat]         = new RooGaussian("ctr_fttbar_"+tag,"f_{t#bar{t}}^{"+tag+"} constrain",*model.fttbar[icat],*model.fttbar_mean_constrain[icat],*model.fttbar_sigma_constrain[icat]);
      
      if(fitEq ||  fitType_==FIT_R_AND_EB)  model.fsingletop[icat] = new RooRealVar("fsingletop_"+tag,"f_{t}^{"+tag+"}",fsingletop_[tag]);
      else                                  model.fsingletop[icat] = new RooRealVar("fsingletop_"+tag,"f_{t}^{"+tag+"}",fsingletop_[tag],fsingletopUnc_[tag]-3*fsingletopUnc_[tag],fsingletop_[tag]+3*fsingletopUnc_[tag]);
      model.fsingletop_mean_constrain[icat]  = new RooRealVar("meanfsingletop_"+tag,"#bar{f_{t}^{"+tag+"}}",fsingletop_[tag]);
      model.fsingletop_sigma_constrain[icat] = new RooRealVar("uncfsingletop_"+tag,"#sigma_{f_{t}^{"+tag+"}}",fsingletopUnc_[tag]);
      model.fsingletop_constrain[icat]       = new RooGaussian("ctr_fsingletop_"+tag,"f_{t}^{"+tag+"} constrain",*model.fsingletop[icat],*model.fsingletop_mean_constrain[icat],*model.fsingletop_sigma_constrain[icat]);

      //create the model
      HeavyFlavorPDF *mhfc = new HeavyFlavorPDF("hfcmodel_"+tag,"hfcmodel_"+tag,*model.bmult,*model.r,*model.eb,*model.eq,*model.fcorrect[icat],*model.fttbar[icat],*model.fsingletop[icat],*model.jetocc[icat]);

      TString dilCategory(tag);
      if(tag.Contains("ee"))   dilCategory="ee";
      if(tag.Contains("emu"))  dilCategory="emu";
      if(tag.Contains("mumu")) dilCategory="mumu";
      TString jmultkey=tag;
      jmultkey=jmultkey.ReplaceAll(dilCategory,"");
      int jmult=jmultkey.Atoi();
      mhfc->setJetMultiplicity(jmult);
      model.pdfSet.add( *mhfc );
      
      //add the constraints
      RooProdPdf *modelconstr=0;
      if(fitType_==FIT_EB || fitType_==FIT_R)
	{
	  modelconstr=new RooProdPdf("modelconstr_"+tag,"model x product of constrains",RooArgSet(*mhfc,*model.fttbar_constrain[icat],*model.fsingletop_constrain[icat]));
	}
      else if (fitType_==FIT_R_AND_EB) 
	{
	  modelconstr=new RooProdPdf("modelconstr_"+tag,"model x product of constrains",RooArgSet(*mhfc));
	}
      else if(fitType_==FIT_EB_AND_EQ)
	modelconstr=new RooProdPdf("modelconstr_"+tag,"model x product of constrains",RooArgSet(*mhfc,*model.sfeq_constrain));
	  
      model.constrPDFSet.add( *modelconstr );
    }
    
  //all done here
  isInit_=true;
}


//
void HFCMeasurement::resetModelValues()
{
  if(!isInit_) return;
  model.r->setVal(smR_);
  model.abseb->setVal(effb_);
  model.sfeb->setVal(sfb_);
  model.abseq->setVal(effq_);
  model.sfeq->setVal(sfq_);
  int icat=0;
  for(std::set<TString>::iterator cIt=categoryKeys_.begin(); cIt!=categoryKeys_.end(); cIt++,icat++)
    {
      TString tag=*cIt;
      model.fcorrect[icat]->setVal(fcorrect_[tag]);
      model.fttbar[icat]->setVal(fttbar_[tag]);
      model.fsingletop[icat]->setVal(fsingletop_[tag]);
    }
}


//
void HFCMeasurement::fitHFCtoEnsemble(top::EventSummaryHandler &evHandler)
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

      //the physics
      top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      
      //check dilepton
      LorentzVector dilepton = phys.leptons[0]+phys.leptons[1];
      float dilcharge=(phys.leptons[0].id/fabs(phys.leptons[0].id)) *(phys.leptons[1].id/fabs(phys.leptons[1].id));
      if((ev.cat==EE || ev.cat==MUMU) && (fabs(dilepton.mass()-91)<15 || phys.met.pt()<30)) continue;
      if(dilcharge>0) continue;

      //check jets: kinematics + count b-tags
      int njets=0;
      int nbtags=0;
      for(unsigned int ijet=0; ijet<phys.jets.size(); ijet++)
	{
	  if(phys.jets[ijet].pt()<30 || fabs(phys.jets[ijet].eta())>2.4) continue;

	  njets++;
	  double btag(-9999.);
	  if(btagAlgo_.Contains("TCHE") )        btag = phys.jets[ijet].btag1;
	  else if(btagAlgo_.Contains("TCHP") )   btag = phys.jets[ijet].btag2;
	  else if(btagAlgo_.Contains("SSVHE") )  btag = phys.jets[ijet].btag3;
	  else if(btagAlgo_.Contains("JBP") )    btag = phys.jets[ijet].btag4;
	  else if(btagAlgo_.Contains("JP") )     btag = phys.jets[ijet].btag5;
	  else if(btagAlgo_.Contains("SSVHP") )  btag = phys.jets[ijet].btag6;
	  else if(btagAlgo_.Contains("CSV") )    btag = phys.jets[ijet].btag7;
	  nbtags += (btag>algoCut_);
	}
      if(njets>maxJets_) continue;

      std::vector<TString> dilCategories;
      dilCategories.push_back("");
      if(ev.cat==EE)   dilCategories.push_back("ee");
      if(ev.cat==MUMU) dilCategories.push_back("mumu");
      if(ev.cat==EMU)  dilCategories.push_back("emu");

      TString tag("btags_"); tag += njets;

      for(size_t icat=0; icat<dilCategories.size(); icat++)
	{
	  TString ctf= dilCategories[icat];
	  controlHistos_.fillHisto("btags",ctf,nbtags);
	  controlHistos_.fillHisto(tag,ctf,nbtags);
	  controlHistos_.fillHisto("avg"+tag, ctf, nbtags );      
	  controlHistos_.fillHisto("avgbtags", ctf, nbtags );      
	}
    }
  
  //run the fit
  runHFCFit();
}


//
void HFCMeasurement::runHFCFit()
{

  //all data
  RooDataHist* alldata = new RooDataHist("data","data", RooArgList(*model.bmult));

  //build the likelihoods per category
  RooArgSet allLL, jetOccupancies;  
  TIterator *pdfIt     = model.constrPDFSet.createIterator();
  // TIterator *modelpdfIt     = model.pdfSet.createIterator();
  std::vector<RooNLLVar *> theLLsPerCategory;
  int icat(0);
  for(std::set<TString>::iterator cIt = categoryKeys_.begin(); cIt != categoryKeys_.end(); cIt++, icat++)
    {
      //decode this in channel and jet multiplicity
      TString tag = *cIt;

      TString dilCategory(tag);
      if(tag.Contains("ee"))   dilCategory="ee";
      if(tag.Contains("emu"))  dilCategory="emu";
      if(tag.Contains("mumu")) dilCategory="mumu";
      
      TString jmultkey=tag;
      jmultkey = jmultkey.ReplaceAll(dilCategory,"");
            
      //convert to a binned dataset
      TH1 *h = controlHistos_.getHisto("btags_"+jmultkey,dilCategory);
      RooDataHist* ds = new RooDataHist("data"+tag,"data"+tag, RooArgList(*model.bmult), h);
      alldata->add(*ds);

      //set the overall normalization for this category
      model.jetocc[icat]->setVal(h->Integral());
      jetOccupancies.add(*(model.jetocc[icat]));

      //build the likelihood
      RooNLLVar *nll=0;
      RooProdPdf *modelconstr = dynamic_cast<RooProdPdf *>( pdfIt->Next() );
      if(fitType_==FIT_EB || fitType_==FIT_R || fitType_==FIT_R_AND_EB) 
	{
	  nll = (RooNLLVar *) modelconstr->createNLL(*ds,Constrain(RooArgSet(*model.fttbar[icat],*model.fsingletop[icat])),NumCPU(2));
	}
      else if(  fitType_==FIT_R_AND_EB)
	{
	  nll = (RooNLLVar *) modelconstr->createNLL(*ds,NumCPU(2));
	}
      else if(fitType_==FIT_EB_AND_EQ /*|| fitType_==FIT_R_AND_EQ*/)  
	{
	  nll = (RooNLLVar *) modelconstr->createNLL(*ds,Constrain(RooArgSet(*model.sfeq)),NumCPU(2));
	}

      TString itit="="+jmultkey + " jets";
      nll->SetTitle(itit);
      allLL.add(*nll); 
      theLLsPerCategory.push_back(nll);
    }

  // add up all the log likelihoods and fit them 
  // reinforce that the error is to be taken from  +/- 1/2
  RooAddition *combll = new RooAddition("combll","combll",allLL);      
  RooMinuit minuit(*combll); 
  minuit.migrad();
  minuit.setErrorLevel(0.5);     
  minuit.hesse();
  minuit.minos();
  minuit.setErrorLevel(0.5);     
  RooFitResult *r=minuit.save();

  printf("\n\n#varepsilon_{b}=%3.3f #pm %3.2f\n\n",model.abseb->getVal()*model.sfeb->getVal(),model.abseb->getVal()*model.sfeb->getError());

  //draw the result
  if(true)
    {
      //
      // SUMMARY OF RESULTS
      //
      TCanvas *c = new TCanvas("combination","combination",600,600);
      c->SetWindowSize(600,600);

      //plot the likelihood in the main frame
      RooPlot *frame = 0;

      TString resLabel;
      if(fitType_==FIT_EB || fitType_==FIT_EB_AND_EQ) 
	{
	  frame = model.sfeb->frame(Title("Likelihood"),Range(0.8,1.1));
	  frame->GetXaxis()->SetTitle("SF #varepsilon_{b}");
	  char buf[100];
	  sprintf(buf,"SF #varepsilon_{b}=%3.3f #pm %3.2f, ",model.sfeb->getVal(),model.sfeb->getError());
	  resLabel += buf;
	  char buf2[100];
	  sprintf(buf2,"#varepsilon_{b}=%3.3f #pm %3.2f",model.abseb->getVal()*model.sfeb->getVal(),model.abseb->getVal()*model.sfeb->getError());
	  resLabel +=buf2;
	}
      else if (fitType_==FIT_R || fitType_==FIT_R_AND_EB) 
	{
	  frame = model.r->frame(Title("Likelihood"),Range(0.7,1.4)) ;    
	  frame->GetXaxis()->SetTitle("R");
	  char buf[100];
	  sprintf(buf,"R=%3.2f #pm %3.2f",model.r->getVal(),model.r->getError());
	  resLabel += buf;
	}
      combll->plotOn(frame,ShiftToZero(),Name("ll"));

      TLegend *leg = new TLegend(0.2,0.65,0.45,0.9,NULL,"brNDC");
      formatForCmsPublic(c,leg,"",5);
      leg->AddEntry("ll","Combined","l");

      TPaveText *pave = new TPaveText(0.15,0.96,0.51,0.99,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      pave->SetTextAlign(12);
      pave->SetTextFont(42);
      pave->AddText("CMS preliminary");
      pave->Draw();
      
      pave = new TPaveText(0.5,0.96,0.94,0.99,"NDC");
      pave->SetFillStyle(0);
      pave->SetBorderSize(0);
      pave->SetTextAlign(32);
      pave->SetTextFont(42);
      pave->AddText(resLabel.Data());
      pave->Draw();
      
      //debug
      cout << "****** SUMMARY RESULT ******" << endl
	   << resLabel << endl
	   << "EDM = " << r->edm() << endl
	   << "-log(L) at minimum = " << r->minNll() << endl   
	   << "final value of floating parameters" << endl ;
      r->floatParsFinal().Print("s") ;
      cout << "****************************" << endl;

      //superimpose the other likelihoods now
      int icat(0);
      for(std::set<TString>::iterator cIt=categoryKeys_.begin(); cIt!=categoryKeys_.end(); cIt++,icat++)
	{
	  TString tag=*cIt;
	  TString lltitle("");
	  TString dilCategory("");
	  int lcolor(1), lstyle(1);
	  if(tag.Contains("ee"))   { dilCategory="ee";  lltitle="ee";      lcolor=kGreen+4;}
	  if(tag.Contains("emu"))  { dilCategory="emu"; lltitle="e#mu";    lcolor=kRed;    }
	  if(tag.Contains("mumu")) { dilCategory="mumu"; lltitle="#mu#mu"; lcolor=kOrange; }
	  TString jmult=tag;
	  jmult=jmult.ReplaceAll(dilCategory,"");
	  if(jmult=="3") lstyle=9;
	  lltitle += " (" + jmult + " jets)";
	  theLLsPerCategory[icat]->plotOn(frame, ShiftToZero(),  LineColor(lcolor), LineStyle(lstyle), LineWidth(2),  MoveToBack(), Name(tag+"ll") );
	  leg->AddEntry(tag+"ll",lltitle,"l");
	}
      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetYaxis()->SetTitle("-Log(L/L_{Max})");
      frame->GetYaxis()->SetTitleOffset(1);
      frame->Draw();
      leg->SetFillColor(0);
      leg->SetFillStyle(3001);
      leg->Draw();
      
      //draw the data and the model sum in a sub-pad
      TPad *npad = new TPad("llpad","ll", 0.6, 0.6, 0.9, 0.9);
      npad->Draw();
      npad->cd();
      frame = model.bmult->frame();

      RooAddPdf incmodel("incmodel","Inclusive Model",model.pdfSet,jetOccupancies);
      alldata->plotOn(frame,Binning(maxJets_+1),DrawOption("pz"));
      incmodel.plotOn(frame,Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
      frame->Draw();
      c->Modified();
      c->Update();
      c->SaveAs("HFCLikelihood.png");
      c->SaveAs("HFCLikelihood.pdf");
      c->SaveAs("HFCLikelihood.C");

      //
      // CONTOUR PLOT FOR COMBINED FITS
      //
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
	  formatForCmsPublic(c,0,"",3);
	  c->Modified();
	  c->Update();
	  c->SaveAs("HFCContour.C");
	  c->SaveAs("HFCContour.pdf");
	  c->SaveAs("HFCContour.png");

	}
    }
}


//
void HFCMeasurement::runHFCDiffFit(TString dilCategory)
{

  //all data
  RooDataHist* alldata = new RooDataHist("data","data", RooArgList(*model.bmult));

  //build the likelihoods per category
  RooArgSet allLL, jetOccupancies;  
  TIterator *pdfIt       = model.constrPDFSet.createIterator();
  //TIterator *modelpdfIt  = model.pdfSet.createIterator();
  std::vector<RooNLLVar *> theLLsPerCategory;
  for(int jmult=2; jmult<=maxJets_; jmult++)
    {
      TString tag(""); tag += jmult;
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
      if(fitType_==FIT_EB || fitType_==FIT_R || fitType_==FIT_R_AND_EB) nll = (RooNLLVar *) modelconstr->createNLL(*ds,Constrain(RooArgSet(*model.fttbar[jmult-2],*model.fsingletop[jmult-2])));
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
      TLegend *leg=c->BuildLegend();
      formatForCmsPublic(c,leg,label,1);
      leg->Delete();

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
void HFCMeasurement::printConfiguration(std::ostream &os)
{
  os << "[HFCMeasurement::printConfiguration]" << endl;
  os << "Fit type:" << fitTypeTitle_ << endl;
  os << "Using: 2<= N_{jets}<=" << maxJets_ << endl;
  os << "Counting b-tags with: " << btagAlgo_ << " cutting @ " << algoCut_ << endl;
  os << "\t epsilon_b:" << effb_ << " SF_b: " << sfb_ << " +/- " << sfbUnc_ << endl;
  os << "\t epsilon_q:" << effq_ << " SF_q: " << sfq_ << " +/- " << sfqUnc_ << endl;
  os << "Data will be split in: " << categoryKeys_.size() << " categories" << endl;
  for(std::set<TString>::iterator cit=categoryKeys_.begin(); cit!=categoryKeys_.end(); cit++)
    {
      TString tag(*cit);
      TString dilCategory(tag);
      if(tag.Contains("ee"))   dilCategory="ee";
      if(tag.Contains("emu"))  dilCategory="emu";
      if(tag.Contains("mumu")) dilCategory="mumu";
      TString jmult=tag;
      jmult = jmult.ReplaceAll(dilCategory,"");

      os << "Category: " << dilCategory << "(=" << jmult << " jets)" << endl; 
      os << "\t f_{correct}:" << fcorrect_[tag] << " +/- " << fcorrectUnc_[tag] << endl;
      os << "\t f_{ttbar}:" << fttbar_[tag] << " +/- " << fttbarUnc_[tag] << endl;
      os << "\t k_{single top}:" << fsingletop_[tag] << " +/- " << fsingletopUnc_[tag] << endl;
    }
}
