#include "LIP/Top/interface/HFCMeasurement.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ModelConfig.h"

#include "RooNumIntConfig.h"
#include "RooNLLVar.h"
#include "RooConstVar.h"
#include "RooExtendPdf.h"

#include  "TGraphErrors.h"

using namespace RooFit;
using namespace RooStats;
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
	  bool doAverage(hname.Contains("avg") && nMeasurements_>0);
	  if(doAverage) hit->second->Scale(1./nMeasurements_);
	  hit->second->Write();
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
  model.bmult    = new RooRealVar("bmult","N_{btags}",0.,float(maxJets_+1));
  model.bmult->setBins(maxJets_+1);
  model.bmultObs = new RooRealVar("bmultobs","events",0.,9999999999999999999.);

  //what do we want to fit from the observables
  bool fitR(fitType_==FIT_R || fitType_==FIT_R_AND_XSEC || fitType_==FIT_R_AND_EB || fitType_==FIT_R_CONSTRAINED);
  bool fitEb(fitType_==FIT_EB || fitType_==FIT_R_AND_EB || fitType_==FIT_EB_AND_XSEC || fitType_==FIT_EB_AND_EQ);
  bool fitEq(fitType_==FIT_EB_AND_EQ);

  //R=B(t->Wb)/B(t->Wq)
  if(!fitR)                              model.r = new RooRealVar("r","R",smR_);
  else if(fitType_==FIT_R_AND_EB)        model.r = new RooRealVar("r","R",1.0,0.96,1.06);
  else if(fitType_==FIT_R_CONSTRAINED)   model.r = new RooRealVar("r","R",1.0,0.,1.0);
  else                                   model.r = new RooRealVar("r","R",1.0,0.,2.0);
	

  //b-tagging effiency
  model.abseb                = new RooRealVar("abseb","abs#varepsilon_{b}",effb_);
  if(fitType_==FIT_R_AND_EB) model.sfeb    = new RooRealVar("sfeb","SF #varepsilon_{b}",sfb_,0.85,1.05); 
  else if(fitEb)             model.sfeb    = new RooRealVar("sfeb","SF #varepsilon_{b}",sfb_,0.7,min(1./effb_,1.3));
  else                     { model.sfeb    = new RooRealVar("sfeb","SF #varepsilon_{b}",sfb_);                          model.sfeb->setError( sfbUnc_ );  }
  model.sfeb_constrain       = new RooGaussian("sfeb_constrain","#varepsilon_{b} constrain",*model.sfeb,RooConst(sfb_),RooConst(sfbUnc_));
  model.eb                   = new RooFormulaVar("eb","@0*@1",RooArgSet(*model.abseb,*model.sfeb));
  
  //mistag efficiency 
  model.abseq = new RooRealVar("abseq","abs#varepsilon_{q}",effq_);
  if(fitEq)   model.sfeq       = new RooRealVar("sfeq","SF #varepsilon_{q}",sfq_,0.5,min(1./effq_,1.5));
  else      { model.sfeq       = new RooRealVar("sfeq","SF #varepsilon_{q}",sfq_);                                      model.sfeq->setError(sfqUnc_); }
  model.sfeq_constrain         = new RooGaussian("sfeq_constrain","#varepsilon_{q} constrain",*model.sfeq,RooConst(sfq_),RooConst(sfqUnc_));
  model.eq                     = new RooFormulaVar("eq","@0*@1",RooArgSet(*model.abseq,*model.sfeq));
  if(fitType_==FIT_EB_AND_EQ)    model.pdfConstrains.add(*model.sfeq_constrain);

  //this is a correction for acceptance (keep constant for now)
  model.acc1  = new RooRealVar("accr1","A(R=1)",1.0);
  model.acc05 = new RooRealVar("accr05","A(R=0.5)",1.0);
  model.acc0  = new RooRealVar("accr0","A(R=0)",1.0);

  //the fit is divided in categories (i b-tags, j jets)  : add exclusive jet multiplicity models
  int icat(0);
  model.sample = new RooCategory("sample","sample");
  for(std::set<TString>::iterator cIt=categoryKeys_.begin(); cIt!=categoryKeys_.end(); cIt++,icat++)
    {
      TString tag=*cIt;

      //define the categories for this sub-sample
      if(icat==EE_3JETS || icat==MUMU_3JETS || icat==EMU_3JETS)
	{
	  model.sample->defineType("n0btags_"+tag);
	  model.sample->defineType("n1btags_"+tag);
	  model.sample->defineType("n2btags_"+tag);
	  model.sample->defineType("n3btags_"+tag);
	}
      else
	{
	  model.sample->defineType("n0btags_"+tag);
	  model.sample->defineType("n1btags_"+tag);
	  model.sample->defineType("n2btags_"+tag);
	}

      //occupancy of this category
      model.jetocc[icat] = new RooRealVar("jetocc"+tag,"occ"+tag,1);
      
      //fraction of correct assignments
      bool addFcorrectAsConstrain(true);
      if(fitEq || fitType_==FIT_R_AND_EB)  { model.fcorrect[icat] = new RooRealVar("fcorrect_"+tag,"f_{correct}^{"+tag+"}",fcorrect_[tag]); addFcorrectAsConstrain=false; }
      else                                 model.fcorrect[icat] = new RooRealVar("fcorrect_"+tag,"f_{correct}^{"+tag+"}",fcorrect_[tag],fcorrect_[tag]-3*fcorrectUnc_[tag],fcorrect_[tag]+3*fcorrectUnc_[tag]);
      if(nuisanceType_==GAUSSIAN)
	{
	  model.fcorrect_constrain[icat]       = new RooGaussian("ctr_fcorrect_"+tag,"f_{correct}^{"+tag+"} constrain",*model.fcorrect[icat],RooConst(fcorrect_[tag]),RooConst(fcorrectUnc_[tag]));
	}
      else if(nuisanceType_==UNIFORM)
	{
	  //	  double centralVal=fcorrect_[tag];
	  //	  double sigma=fcorrectUnc_[tag];
	  //	  model.fcorrect[icat]->setRange(max(0.,centralVal-sigma),min(centralVal-sigma,1.));
	  //	  model.fcorrect_constrain[icat]       = new RooUniform("ctr_fcorrect_"+tag,"f_{correct}^{"+tag+"} constrain",*model.fcorrect[icat]);
	  addFcorrectAsConstrain=false;
	}
      else if(nuisanceType_==LOGNORMAL)
	{
	  model.fcorrect_constrain[icat]       = new RooLognormal("ctr_fcorrect_"+tag,"f_{correct}^{"+tag+"} constrain",*model.fcorrect[icat],RooConst(fcorrect_[tag]),RooConst(exp(fcorrectUnc_[tag])));
	}
      if(addFcorrectAsConstrain) model.pdfConstrains.add(*model.fcorrect_constrain[icat]);

      //ttbar fraction in the sample
      bool addFttbarAsContrain(true);
      if(fttbar_[tag]==1.0 || fitEq ||  fitType_==FIT_R_AND_EB)  { model.fttbar[icat] = new RooRealVar("fttbar_"+tag,"f_{t#bar{t}}^{"+tag+"}",fttbar_[tag]); addFttbarAsContrain=false; }
      else                                                         model.fttbar[icat] = new RooRealVar("fttbar_"+tag,"f_{t#bar{t}}^{"+tag+"}",fttbar_[tag],fttbar_[tag]-3*fttbarUnc_[tag],fttbar_[tag]+3*fttbarUnc_[tag]);
      model.fttbar_constrain[icat]         = new RooGaussian("ctr_fttbar_"+tag,"f_{t#bar{t}}^{"+tag+"} constrain",*model.fttbar[icat],RooConst(fttbar_[tag]),RooConst(fttbarUnc_[tag]));
      if(addFttbarAsContrain) model.pdfConstrains.add(*model.fttbar_constrain[icat]);

      //single top fraction in the sample
      bool addFsingleTopAsContrain(true);
      if(fsingletop_[tag]==0 || fitEq ||  fitType_==FIT_R_AND_EB)  { model.fsingletop[icat] = new RooRealVar("fsingletop_"+tag,"f_{t}^{"+tag+"}",fsingletop_[tag]); addFsingleTopAsContrain=false; }
      else                                                           model.fsingletop[icat] = new RooRealVar("fsingletop_"+tag,"f_{t}^{"+tag+"}",fsingletop_[tag],fsingletopUnc_[tag]-3*fsingletopUnc_[tag],fsingletop_[tag]+3*fsingletopUnc_[tag]);
      model.fsingletop_constrain[icat]       = new RooGaussian("ctr_fsingletop_"+tag,"f_{t}^{"+tag+"} constrain",*model.fsingletop[icat],RooConst(fsingletop_[tag]),RooConst(fsingletopUnc_[tag]));
      if(addFsingleTopAsContrain) model.pdfConstrains.add(*model.fsingletop_constrain[icat]);
      
      //sample composition variables i.e. alpha and alpha_i
      float nPairsVal=4;
      if(icat==EE_3JETS || icat==MUMU_3JETS || icat==EMU_3JETS) nPairsVal=6;
      RooRealVar *npairs = new RooRealVar("npairs_"+tag,"npairs_"+tag,nPairsVal);
      model.alpha[icat]  = new RooFormulaVar("alpha_"+tag,  "(@0*@1)/(@2*(2.+@3))",   RooArgSet(*npairs, *model.fcorrect[icat], *model.fttbar[icat], *model.fsingletop[icat]));
      model.alpha2[icat] = new RooFormulaVar("alpha2_"+tag, "pow(@0,2)*@1",           RooArgSet(*model.alpha[icat],*model.fttbar[icat]));
      model.alpha1[icat] = new RooFormulaVar("alpha1_"+tag, "(2*@0*(1-@0)+@0*@2)*@1", RooArgSet(*model.alpha[icat],*model.fttbar[icat],*model.fsingletop[icat]));
      model.alpha0[icat] = new RooFormulaVar("alpha0_"+tag, "1-@0-@1",                RooArgSet(*model.alpha2[icat],*model.alpha1[icat]));
      
      //modify locally the b-tag/mistag rates with correction factors
      model.ebcorr[icat] = new RooFormulaVar("eb_"+tag,"@0*@1", RooArgSet(*model.eb,RooConst(btagEffCorr_[tag])));
      model.eqcorr[icat] = new RooFormulaVar("eq_"+tag,"@0*@1", RooArgSet(*model.eq,RooConst(ltagEffCorr_[tag])));
      
      //define the base PDFs (kernel)
      RooArgList stdArgList;
      stdArgList.add(*model.r);                     // @0  R
      stdArgList.add(*model.ebcorr[icat]);          // @1  eff_b
      stdArgList.add(*model.eqcorr[icat]);          // @2  eff_q
      stdArgList.add(*model.alpha2[icat]);          // @3  alpha_2
      stdArgList.add(*model.alpha1[icat]);          // @4  alpha_1
      stdArgList.add(*model.alpha0[icat]);          // @5  alpha_0 
      stdArgList.add(*model.acc1);                  // @6  acc_1
      stdArgList.add(*model.acc05);                 // @7  acc05
      stdArgList.add(*model.acc0);                  // @8  acc0
      
      RooGenericPdf *basepdf0btags2jets = new RooGenericPdf("kernelpdf0btags_"+tag,
							    "@3*( pow(@0*(1-@1),2)*@6 + 2*@0*(1-@0)*(1-@1)*(1-@2)*@7 + pow((1-@0)*(1-@2),2)*@8) +"
							    "@4*( pow(@0,2)*(1-@1)*(1-@2)*@6 + @0*(1-@0)*((1-@1)+(1-@2))*(1-@2)*@7 + pow((1-@0)*(1-@2),2)*@8)+"
							    "@5*( pow(1-@2,2) )",
							    stdArgList);
      
      RooGenericPdf *basepdf1btags2jets = new RooGenericPdf("kernelpdf1btags_"+tag,
							    "@3*( 2*pow(@0,2)*(1-@1)*@1*@6 + 2*@0*(1-@0)*((1-@1)*@2+@1*(1-@2))*@7 + 2*pow((1-@0),2)*(1-@2)*@2*@8) +"
							    "@4*( pow(@0,2)*(@1*(1-@2)+(1-@1)*@2)*@6 + @0*(1-@0)*(@1*(1-@2)+(1-@1)*@2+2*@2*(1-@2))*@7 + 2*pow((1-@0),2)*@2*(1-@2)*@8 ) +"
							    "@5*( 2*@2*(1-@2) )",
							    stdArgList);

      RooGenericPdf *basepdf2btags2jets = new RooGenericPdf("kernelpdf2btags_"+tag,
							    "@3*( pow(@0*@1,2)*@6 + 2*@0*(1-@0)*@1*@2*@7 + pow((1-@0)*@2,2)*@8 ) + "
							    "@4*( pow(@0,2)*@1*@2*@6 + @0*(1-@0)*((@1+@2)*@2)*@7 + pow((1-@0)*@2,2)*@8 ) +"
							    "@5*( pow(@2,2) )",
							    stdArgList);
      
      //check if extension is needed
      if(icat==EE_3JETS || icat==MUMU_3JETS || icat==EMU_3JETS)
	{
	  //b-tags from extra jets are assumed to be mistags
	  RooGenericPdf *basepdf0btagsextrajet = new RooGenericPdf("pdf0btagsextrajet_"+tag,"(1-@0)",RooArgList(*model.eq));
	  RooGenericPdf *basepdf1btagsextrajet = new RooGenericPdf("pdf1btagsextrajet_"+tag,"@0",RooArgList(*model.eq));
	  
	  //
	  RooGenericPdf *basepdf0btags3jets = new RooGenericPdf("basepdf0btags3Jets_"+tag, "@0*@1",       RooArgList(*basepdf0btags2jets, *basepdf0btagsextrajet) );  
	  RooGenericPdf *basepdf1btags3jets = new RooGenericPdf("basepdf1btags3Jets_"+tag, "@0*@1+@2*@3", RooArgList(*basepdf0btags2jets, *basepdf1btagsextrajet, *basepdf1btags2jets, *basepdf0btagsextrajet) );
	  RooGenericPdf *basepdf2btags3jets = new RooGenericPdf("basepdf2btags3Jets_"+tag, "@0*@1+@2*@3", RooArgList(*basepdf1btags2jets, *basepdf1btagsextrajet, *basepdf2btags2jets, *basepdf0btagsextrajet) );
	  RooGenericPdf *basepdf3btags3jets = new RooGenericPdf("basepdf3btags3Jets_"+tag, "@0*@1",       RooArgList(*basepdf2btags2jets, *basepdf1btagsextrajet) );
	  
	  //add for the different categories
	  model.pdfForCategory["n0btags_"+tag] = basepdf0btags3jets;
	  model.pdfForCategory["n1btags_"+tag] = basepdf1btags3jets;
	  model.pdfForCategory["n2btags_"+tag] = basepdf2btags3jets;
	  model.pdfForCategory["n3btags_"+tag] = basepdf3btags3jets;
	}
      else
	{
	  //add for the different categories
	  model.pdfForCategory["n0btags_"+tag] = basepdf0btags2jets;
	  model.pdfForCategory["n1btags_"+tag] = basepdf1btags2jets;
	  model.pdfForCategory["n2btags_"+tag] = basepdf2btags2jets;
	}
    }


  // the model will be a simultaneous PDF
  RooSimultaneous *simPdf = new RooSimultaneous("simmodel","simmodel",*model.sample);
  for(std::map<TString, RooAbsPdf *>::iterator it = model.pdfForCategory.begin();   it != model.pdfForCategory.end(); it++)  simPdf->addPdf( *(it->second), it->first );
  model.eventYields = new RooRealVar("eventyields","Event yields",1,0,9999999.); 
  model.pdf = new RooExtendPdf("hfcmodel","hfcmodel",*simPdf,*model.eventYields);
      
  //this is the constrained model with nuisance parameters
  RooArgSet argSet(*model.pdf);  argSet.add(model.pdfConstrains);
  model.constrPdf = new RooProdPdf("hfcconstrmodel","hfconstrmodel",argSet);

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
      model.jetocc[icat]->setVal(0);
    }
}

//
void HFCMeasurement::fitHFCtoMeasurement(std::vector<TH1D *> &btagHistos, int runMode, bool debug)
{
  if(btagHistos.size()==0) return;

  //restart all over again
  initHFCModel();
  resetModelValues();
  resetHistograms();

  //add the histograms
  std::map<int,TString> cats;
  cats[EE_2JETS]="ee";
  cats[EE_3JETS]="ee";
  cats[MUMU_2JETS]="mumu";
  cats[MUMU_3JETS]="mumu";
  cats[EMU_2JETS]="emu";
  cats[EMU_3JETS]="emu";
  for(std::map<int,TString>::iterator cit = cats.begin(); cit!= cats.end(); cit++)
    {
      TString tag("");
      if(cit->first==EE_2JETS || cit->first==MUMU_2JETS ||cit->first==EMU_2JETS) tag="btags_2";
      if(cit->first==EE_3JETS || cit->first==MUMU_3JETS ||cit->first==EMU_3JETS) tag="btags_3";
      TH1 *bmh=btagHistos[cit->first];
      TString ctf=cit->second;
      controlHistos_.getHisto("btags",ctf)->Add(bmh);
      controlHistos_.getHisto(tag,ctf)->Add(bmh);
      controlHistos_.getHisto("avg"+tag, ctf)->Add(bmh);
      controlHistos_.getHisto("avgbtags", ctf)->Add(bmh);
    }
  
  //run the fit
  nMeasurements_++;
  runHFCFit(runMode,debug);
}


//
void HFCMeasurement::fitHFCtoEnsemble(top::EventSummaryHandler &evHandler, int runMode, bool debug )
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
  runHFCFit(runMode,debug);
}


//
void HFCMeasurement::runHFCFit(int runMode, bool debug)
{
  
  //build a categorized dataset
  int icat(0);
  RooDataSet *data = new RooDataSet("data","data",RooArgSet( *model.bmult, *model.bmultObs, *model.sample), RooFit::WeightVar(*model.bmultObs) );
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
      int njets=jmultkey.Atoi();
            
      //each bin in the histogram is a exclusive category
      TH1 *h = controlHistos_.getHisto("btags_"+jmultkey,dilCategory);
      for(int ibin=1; ibin<= h->GetXaxis()->GetNbins(); ibin++)
	{
	  if(ibin==njets+2) break;
	  TString catname("n"); catname += (ibin-1); catname += "btags_" + tag;
	  model.bmult->setVal(ibin-1);
	  model.bmultObs->setVal(h->GetBinContent(ibin));
	  model.sample->setLabel(catname);
	  data->add( RooArgSet(*model.bmult,*model.bmultObs,*model.sample), model.bmultObs->getVal() );
	}
      
      //set overall occupancy
      model.jetocc[icat]->setVal( h->Integral() );
    }
  

  //
  //fit individually each channel
  //
  std::vector<RooNLLVar *> likelihoods;
  std::vector<TString> lltitles; 
  TGraphAsymmErrors *exclusiveRgr=0;
  std::vector<TPaveText *> exclusiveFitNames;
  if(runMode==FITEXCLUSIVECATEGORIES)
    {
      exclusiveRgr = new TGraphAsymmErrors;
      icat=0;
      for(std::set<TString>::iterator cIt = categoryKeys_.begin(); cIt != categoryKeys_.end(); cIt++, icat++)
	{
	  TString tag=*cIt;
	  TString cut("sample==sample::n0btags_"+tag+" || sample==sample::n1btags_"+tag+" || sample==sample::n2btags_"+tag);
	  if(tag.Contains("3")) cut+=" || sample==sample::n3btags_"+tag;
	  
	  //fit a projection of the data
	  RooDataSet *dataslice = (RooDataSet *) data->reduce(cut);
	  model.constrPdf->fitTo(*dataslice,
				 Constrain(model.pdfConstrains),
				 RooFit::Save(),
				 RooFit::Hesse(kTRUE),
				 RooFit::Minos(kTRUE),
				 RooFit::SumW2Error(kTRUE),
				 RooFit::PrintLevel(-1),
				 RooFit::Verbose(kFALSE));
	  model.rFit[icat]=model.r->getVal();   model.rFitAsymmErrLo[icat]=model.r->getAsymErrorLo(); model.rFitAsymmErrHi[icat]=model.r->getAsymErrorHi();
	  
	  //add to the exclusive fit graph
	  exclusiveRgr->SetPoint(icat,model.rFit[icat],icat*2);
	  exclusiveRgr->SetPointError(icat,fabs(model.rFitAsymmErrLo[icat]),fabs(model.rFitAsymmErrHi[icat]),0.1,0.1);
	  TPaveText *pt = new TPaveText(1.0,icat*2+0.2,1.0,icat*2+0.6,"br");
	  pt->SetBorderSize(0);
	  pt->SetFillColor(0);
	  pt->SetFillStyle(0);
	  pt->SetTextFont(52);
	  pt->SetTextSize(0.03);
	  TString caption=tag;
	  caption=caption.ReplaceAll("mu","#mu");
	  caption=caption.ReplaceAll("2"," (2 jets)");
	  caption=caption.ReplaceAll("3"," (3 jets)");
	  pt->AddText(caption);
	  exclusiveFitNames.push_back(pt);
	  
	  //combined likelihood for the 2+3 jet bins (has to be done after fit to both multipliticty bins)
	  if(tag.Contains("2")) continue;
	  TString tag2=tag.ReplaceAll("3","2");
	  TString cut2("sample==sample::n0btags_"+tag2+" || sample==sample::n1btags_"+tag2+" || sample==sample::n2btags_"+tag2);
	  TString cut3=cut;
	  dataslice = (RooDataSet *) data->reduce(cut2 + " || " + cut3);
	  model.constrPdf->fitTo(*dataslice,
				 Constrain(model.pdfConstrains),
				 RooFit::Save(),
				 RooFit::Hesse(kTRUE),
				 RooFit::Minos(kTRUE),
				 RooFit::SumW2Error(kTRUE),
				 RooFit::PrintLevel(-1),
				 RooFit::Verbose(kFALSE));
	  likelihoods.push_back( (RooNLLVar *) model.constrPdf->createNLL(*dataslice,Constrain(model.pdfConstrains),NumCPU(2)) );
	  TString tit=tag2.ReplaceAll("2","");
	  tit=tit.ReplaceAll("mu","#mu");
	  lltitles.push_back(tit);
	}
    }

  //
  //fit globally
  //
  model.constrPdf->Print();
  model.constrPdf->fitTo(*data,
			 Constrain(model.pdfConstrains),
			 RooFit::Save(),
			 RooFit::Hesse(kTRUE),
			 RooFit::Minos(kTRUE),
			 RooFit::SumW2Error(kTRUE),
			 RooFit::PrintLevel(-1),
			 RooFit::Verbose(kFALSE));
  model.rFitResult = model.r->getVal();
  model.rFitResultAsymmErrLo = model.r->getAsymErrorLo();
  model.rFitResultAsymmErrHi = model.r->getAsymErrorHi();

  RooNLLVar *nll = (RooNLLVar *) model.constrPdf->createNLL(*data,Constrain(model.pdfConstrains),NumCPU(2));
  RooMinuit minuit(*nll); 
  minuit.setErrorLevel(0.5);     
  minuit.hesse();
  minuit.migrad();
  minuit.setErrorLevel(0.5);     
  minuit.hesse();
  minuit.minos();
  minuit.setErrorLevel(0.5);     
  RooFitResult *r=minuit.save();
  model.minLL = r->minNll();
 
  //get the interval if required
  if(fitType_==FIT_R_CONSTRAINED)
    {
      //
      // prepare the configuration for Feldman-Cousins
      //
      //create a workspace
      RooWorkspace* w = new RooWorkspace("w");
      ModelConfig modelConfig("hfcfit",w);
      modelConfig.SetPdf(*model.constrPdf);
      modelConfig.SetParametersOfInterest(RooArgSet(*model.r));
      modelConfig.SetObservables(RooArgSet(*model.bmult,*model.bmultObs,*model.sample));
      
      //configure FC
      RooStats::FeldmanCousins fc(*data,modelConfig);
      RooDataSet poiToTest("poitotest","poitotest",RooArgSet(*model.r));
      double rmin=max(double(0.),double(model.rFitResult+model.rFitResultAsymmErrLo*5));
      double rmax=1.0;
      for(double r=rmin; r<=rmax; r+= (rmax-rmin)/20.)
	{
	  model.r->setVal(r);
	  poiToTest.add(RooArgSet(*model.r));
	}
      fc.SetPOIPointsToTest(poiToTest);
      fc.SetTestSize(.05);               // set size of test 
      fc.FluctuateNumDataEntries(false); // number counting analysis: dataset always has 1 entry with N events observed
      fc.SetConfidenceLevel(0.95);       // 95% interval
      fc.UseAdaptiveSampling(true);      // speed it up a bit
      fc.CreateConfBelt(true);
      
      cout << "Starting Feldman-Cousins computation: you can go and take a loooooong coffee " << endl;
      PointSetInterval* interval = (PointSetInterval*)fc.GetInterval();
      ConfidenceBelt* belt = fc.GetConfidenceBelt();
      model.rFitLowerLimit=interval->LowerLimit(*model.r);
      model.rFitUpperLimit=interval->UpperLimit(*model.r);

      //get the thresholds
      if(debug)
	{

	  //save workspace to file
	  TFile *debugF = TFile::Open("HeavyFlavorWS.root","RECREATE");
	  debugF->cd();
	  w->import(modelConfig);
	  w->import(*data);
	  w->Write();
	  debugF->Close();
	  debugF->Delete();

	  //display in canvas
	  TCanvas *c = new TCanvas("HFCIntervalForFeldmanCousins","HFCIntervalForFeldmanCousins",600,600);
	  c->SetWindowSize(600,600);     
	  TGraph *frame = new TGraph;
	  RooDataSet* parameterScan = (RooDataSet*) fc.GetPointsToScan();
	  for(Int_t i=0; i<parameterScan->numEntries(); ++i)
	    {
	      RooArgSet *tmpPoint = (RooArgSet*) parameterScan->get(i)->clone("temp");
	      double arMax = belt->GetAcceptanceRegionMax(*tmpPoint);
	      double poiVal = tmpPoint->getRealValue(model.r->GetName()) ;	  
	      frame->SetPoint(i,poiVal,arMax);
	    }
	  frame->Draw("ap");
	  frame->GetXaxis()->SetTitle("R");
	  frame->GetYaxis()->SetTitle("Acceptance region max");
	  
	  char buf[100];
	  sprintf(buf,"%3.2f < R < %3.2f",model.rFitLowerLimit,model.rFitUpperLimit);
	  TPaveText *pave = new TPaveText(0.15,0.96,0.41,0.99,"NDC");
	  pave->SetBorderSize(0);
	  pave->SetFillStyle(0);
	  pave->SetTextAlign(12);
	  pave->SetTextFont(42);
	  pave->AddText("CMS preliminary");
	  pave->Draw();
	  
	  pave = new TPaveText(0.4,0.96,0.94,0.99,"NDC");
	  pave->SetFillStyle(0);
	  pave->SetBorderSize(0);
	  pave->SetTextAlign(32);
	  pave->SetTextFont(42);
	  pave->AddText(buf);
	  pave->Draw();
	  
	  c->SaveAs("HFCIntervalForFeldmanCousins.png");
	  c->SaveAs("HFCIntervalForFeldmanCousins.pdf");
	  c->SaveAs("HFCIntervalForFeldmanCousins.C");
	}
    }

  //
  // draw the results
  //
  if(debug)
    {
      RooPlot *frame = model.bmult->frame();

      //exclusive measurements
      TCanvas *c = new TCanvas("HFCMeasurementResults","HFCMeasurementResults",600,600);
      c->SetWindowSize(600,600);
      exclusiveRgr->SetFillStyle(0);
      exclusiveRgr->SetFillColor(0);
      exclusiveRgr->SetLineColor(kRed);
      exclusiveRgr->SetMarkerColor(1);
      exclusiveRgr->SetMarkerStyle(20);
      exclusiveRgr->SetMarkerSize(0.3);
      exclusiveRgr->Draw("ae2p");
      exclusiveRgr->GetXaxis()->SetTitle("R=B(t#rightarrow Wb)/B(t#rightarrow Wq)");
      exclusiveRgr->GetYaxis()->SetNdivisions(0);
      TLine *l=new TLine(1.0,exclusiveRgr->GetYaxis()->GetXmin(),1.0,exclusiveRgr->GetYaxis()->GetXmax());
      l->SetLineColor(kGray);
      l->SetLineStyle(9);
      l->Draw("same");
      for(size_t i=0; i<exclusiveFitNames.size(); i++)  exclusiveFitNames[i]->Draw();
      c->SaveAs("HFCMeasurementResults.png");
      c->SaveAs("HFCMeasurementResults.pdf");
      c->SaveAs("HFCMeasurementResults.C");

      //the likelihoods
      c = new TCanvas("HFCMeasurementLikelihood","HFCMeasurementLikelihood",600,600);
      c->SetWindowSize(600,600);

      //plot the likelihood in the main frame

      TString resLabel;
      if(fitType_==FIT_EB || fitType_==FIT_EB_AND_EQ) 
	{
	  frame = model.sfeb->frame(Title("Likelihood"),Range(0.5,1.3));
	  frame->GetXaxis()->SetTitle("SF #varepsilon_{b}");
	  char buf[100];
	  sprintf(buf,"SF #varepsilon_{b}=%3.3f #pm %3.2f, ",model.sfeb->getVal(),model.sfeb->getError());
	  resLabel += buf;
	  char buf2[100];
	  sprintf(buf2,"#varepsilon_{b}=%3.3f #pm %3.2f",model.abseb->getVal()*model.sfeb->getVal(),model.abseb->getVal()*model.sfeb->getError());
	  resLabel +=buf2;
	}
      else if (fitType_==FIT_R || fitType_==FIT_R_AND_EB || fitType_==FIT_R_CONSTRAINED) 
	{
	  frame = model.r->frame(Title("Likelihood"),Range(0.8,1.2)) ;    
	  frame->GetXaxis()->SetTitle("R");
	  char buf[100];
	  sprintf(buf,"R=%3.2f #pm %3.2f",model.r->getVal(),model.r->getError());
	  resLabel += buf;
	
	  //add the Feldman-Cousins results
	  if(fitType_==FIT_R_CONSTRAINED)
	    {
	      sprintf(buf," (FC: ]%3.2f, %3.2f[)",model.rFitLowerLimit,model.rFitUpperLimit);
	      resLabel += buf;
	    }
	}
  
      nll->plotOn(frame,ShiftToZero(),Name("ll"),FillStyle(0),LineColor(kBlue),LineWidth(2));
      Int_t llColors[]={809,590,824,831};
      for(size_t ill=0; ill<likelihoods.size(); ill++)
	{
	  if(likelihoods[ill]==0) continue;
	  TString llname("ll"); llname +=ill;
	  likelihoods[ill]->plotOn(frame, ShiftToZero(), LineColor(llColors[ill]), LineWidth(2),Name(llname));
	}

      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetYaxis()->SetTitle("-Log(L/L_{Max})");
      frame->GetYaxis()->SetTitleOffset(1);
      frame->GetYaxis()->SetRangeUser(0,5);
      frame->Draw();

      TPaveText *pave = new TPaveText(0.15,0.96,0.41,0.99,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      pave->SetTextAlign(12);
      pave->SetTextFont(42);
      pave->AddText("CMS preliminary");
      pave->Draw();
      
      pave = new TPaveText(0.4,0.96,0.94,0.99,"NDC");
      pave->SetFillStyle(0);
      pave->SetBorderSize(0);
      pave->SetTextAlign(32);
      pave->SetTextFont(42);
      pave->AddText(resLabel.Data());
      pave->Draw();

      TLegend *leg = new TLegend(0.2,0.65,0.45,0.9,NULL,"brNDC");
      leg->AddEntry("ll","Combined","l");     
      for(size_t ill=0; ill<likelihoods.size(); ill++)
	{
	  TString llname("ll"); llname +=ill;
	  leg->AddEntry(llname,lltitles[ill],"l");
	}

      formatForCmsPublic(c,leg,"",5);
      leg->SetFillColor(0);
      leg->SetFillStyle(3001);
      leg->Draw();
      
      //draw the data and the model sum in a sub-pad
      TPad *npad = new TPad("llpad","ll", 0.6, 0.6, 0.9, 0.9);
      npad->Draw();
      npad->cd();
      frame = model.bmult->frame();

      data->plotOn(frame,DrawOption("pz"));
      frame->Draw();
      TH1D *h = new TH1D("model","model",4,0,4);
      h->SetLineWidth(2);
      h->SetLineColor(kBlue);
      int icat(0);
      for(std::set<TString>::iterator cIt=categoryKeys_.begin(); cIt!=categoryKeys_.end(); cIt++,icat++)
	{
	  TString tag=*cIt;
	  for(int ibin=1;ibin<=h->GetXaxis()->GetNbins(); ibin++)
	    {
	      TString catname("n"); catname += (ibin-1); catname += "btags_" + tag;
	      if(model.pdfForCategory.find(catname)==model.pdfForCategory.end()) continue;
	      h->Fill(ibin-1,model.pdfForCategory[catname]->getVal()*model.jetocc[icat]->getVal());
	    }
	}
      h->Draw("histsame");
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
	  c = new TCanvas("HeavyFlavorContour","HeavyFlavorContour");
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
	  c->SaveAs("HeavyFlavorContour.C");
 	  c->SaveAs("HeavyFlavorContour.pdf");
 	  c->SaveAs("HeavyFlavorContour.png");
 	}


      //draw the model per category
      icat=0;
      TGraph *frameGr= new TGraph; frameGr->SetPoint(0,0,0); frameGr->SetPoint(1,1,1);  frameGr->SetMarkerStyle(1);
      Int_t modelColors[]={809,590,824,831};
      for(std::set<TString>::iterator cIt=categoryKeys_.begin(); cIt!=categoryKeys_.end(); cIt++,icat++)
	{
	  TString tag=*cIt;
	  TString jetMult("2"); if(tag.Contains("3")) jetMult="3";
	  TString dilCh=tag; dilCh=dilCh.ReplaceAll(jetMult,"");
	  
	  if(icat%2==0)
	    {
	      c = new TCanvas("HeavyFlavorModel_"+dilCh,"HeavyFlavorContour_"+dilCh);
	      c->SetCanvasSize(1200,600);
	      c->SetWindowSize(1200,600);
	      c->Divide(2,1);
	    }

	  TPad *p=(TPad *) c->cd(icat%2+1);
	  frameGr->Draw("ap");
	  frameGr->GetXaxis()->SetRangeUser(0,1);
	  frameGr->GetXaxis()->SetTitle("R=B(t#rightarrow Wb)/B(t#rightarrow Wq)");
	  frameGr->GetYaxis()->SetRangeUser(0,1);
	  frameGr->GetYaxis()->SetTitle("Probability");

	  TLegend *leg=0;
	  if(icat%2==0) leg = new TLegend(0.2,0.65,0.45,0.9,NULL,"brNDC");
	  
	  for(int ibtags=0;ibtags<=maxJets_; ibtags++)
	    {
	      TString catname("n"); catname += ibtags; catname += "btags_" + tag;
	      if(model.pdfForCategory.find(catname)==model.pdfForCategory.end()) continue;
	      
	      TGraph * modelgr = new TGraph;
	      modelgr->SetName(catname);
	      modelgr->SetLineColor(modelColors[ibtags]);
	      modelgr->SetFillStyle(0);
	      modelgr->SetLineWidth(2);
	      for(float ir=0; ir<=1.0; ir+=0.05)
		{
		  model.r->setVal(ir);
		  int ipt=modelgr->GetN();
		  modelgr->SetPoint(ipt,ir,model.pdfForCategory[catname]->getVal());
		}
	      modelgr->Draw("l");
	      if(leg)
		{
		  TString tit("="); tit+=ibtags; tit+=" b-tags";
		  leg->AddEntry(catname,tit,"l");
		}
	    }
	  	  
	  if(leg) 
	    {
	      leg->Draw();
	      formatForCmsPublic(p,leg,"CMS preliminary",3);
	    }
	  if(icat%2==1)
	    {
	      c->cd();
	      c->Modified();
	      c->Update();
	      c->SaveAs("HeavyFlavorModel_"+dilCh+".C");
	      c->SaveAs("HeavyFlavorModel_"+dilCh+".png");
	      c->SaveAs("HeavyFlavorModel_"+dilCh+".pdf");
	    }
	}
    }

  delete data;
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
