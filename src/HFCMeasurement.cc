#include "LIP/Top/interface/HFCMeasurement.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ModelConfig.h"


#include "RooNumIntConfig.h"
#include "RooNLLVar.h"

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
  model.bmultWgt = new RooRealVar("bmultcts","events",0.,9999999999999999999.);

  //what do we want to fit from the observables
  bool fitR(fitType_==FIT_R || fitType_==FIT_R_AND_XSEC || fitType_==FIT_R_AND_EB);
  bool fitEb(fitType_==FIT_EB || fitType_==FIT_R_AND_EB || fitType_==FIT_EB_AND_XSEC || fitType_==FIT_EB_AND_EQ);
  bool fitEq(fitType_==FIT_EB_AND_EQ);

  //R=B(t->Wb)/B(t->Wq)
  if(!fitR)                        model.r = new RooRealVar("r","R",smR_);
  else if(fitType_==FIT_R_AND_EB)  model.r = new RooRealVar("r","R",1.0,0.9,1.1);
  else                             model.r = new RooRealVar("r","R",1.0,0.,2.0);

  //b-tagging effiency
  model.abseb                = new RooRealVar("abseb","abs#varepsilon_{b}",effb_);
  if(fitEb)    model.sfeb    = new RooRealVar("sfeb","SF #varepsilon_{b}",sfb_,0.7,min(1./effb_,1.3));
  else         model.sfeb    = new RooRealVar("sfeb","SF #varepsilon_{b}",sfb_);
  model.sfeb_mean_constrain  = new RooRealVar("meansfeb","#bar{SF #varepsilon_{b}}",sfb_);
  model.sfeb_sigma_constrain = new RooRealVar("uncsfeb","#sigma_{SF #varepsilon_{b}}",sfbUnc_);
  model.sfeb_constrain       = new RooGaussian("sfeb_constrain","#varepsilon_{b} constrain",*model.sfeb,*model.sfeb_mean_constrain,*model.sfeb_sigma_constrain);
  model.eb                   = new RooFormulaVar("eb","@0*@1",RooArgSet(*model.abseb,*model.sfeb));
  
  //mistag efficiency 
  model.abseq = new RooRealVar("abseq","abs#varepsilon_{q}",effq_);
  if(fitEq) model.sfeq       = new RooRealVar("sfeq","SF #varepsilon_{q}",sfq_,0.5,min(1./effq_,1.5));
  else      model.sfeq       = new RooRealVar("sfeq","SF #varepsilon_{q}",sfq_);
  model.sfeq_mean_constrain  = new RooRealVar("meansfeq","#bar{SF #varepsilon_{q}}",sfq_);
  model.sfeq_sigma_constrain = new RooRealVar("uncsfeq","#sigma_{SF #varepsilon_{q}}",sfqUnc_);
  model.sfeq_constrain       = new RooGaussian("sfeq_constrain","#varepsilon_{q} constrain",*model.sfeq,*model.sfeq_mean_constrain,*model.sfeq_sigma_constrain);
  model.eq                   = new RooFormulaVar("eq","@0*@1",RooArgSet(*model.abseq,*model.sfeq));
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
      //model.fcorrect[icat] = new RooRealVar("fcorrect_"+tag,"f_{correct}^{"+tag+"}",fcorrect_[tag]);
      if(fitEq || fitType_==FIT_R_AND_EB)  model.fcorrect[icat] = new RooRealVar("fcorrect_"+tag,"f_{correct}^{"+tag+"}",fcorrect_[tag]);
      else                                 model.fcorrect[icat] = new RooRealVar("fcorrect_"+tag,"f_{correct}^{"+tag+"}",fcorrect_[tag],fcorrect_[tag]-3*fcorrectUnc_[tag],fcorrect_[tag]+3*fcorrectUnc_[tag]);
      model.fcorrect_mean_constrain[icat]  = new RooRealVar("meanfcorrect_"+tag,"#bar{f_{correct}^{"+tag+"}}",fcorrect_[tag]);
      model.fcorrect_sigma_constrain[icat] = new RooRealVar("uncfcorrect_"+tag,"#sigma_{f_{correct}^{"+tag+"}}",fcorrectUnc_[tag]);
      model.fcorrect_constrain[icat]       = new RooGaussian("ctr_fcorrect_"+tag,"f_{correct}^{"+tag+"} constrain",*model.fcorrect[icat],*model.fcorrect_mean_constrain[icat],*model.fcorrect_sigma_constrain[icat]);
      model.pdfConstrains.add(*model.fcorrect_constrain[icat]);

      //ttbar fraction in the sample
      bool addFttbarAsContrain(true);
      if(fttbar_[tag]==1.0 || fitEq ||  fitType_==FIT_R_AND_EB)  { model.fttbar[icat] = new RooRealVar("fttbar_"+tag,"f_{t#bar{t}}^{"+tag+"}",fttbar_[tag]); addFttbarAsContrain=false; }
      else                                                         model.fttbar[icat] = new RooRealVar("fttbar_"+tag,"f_{t#bar{t}}^{"+tag+"}",fttbar_[tag],fttbar_[tag]-3*fttbarUnc_[tag],fttbar_[tag]+3*fttbarUnc_[tag]);
      model.fttbar_mean_constrain[icat]    = new RooRealVar("meanfttbar_"+tag,"#bar{f_{t#bar{t}}^{"+tag+"}}",fttbar_[tag]);
      model.fttbar_sigma_constrain[icat]   = new RooRealVar("uncfttbar_"+tag,"#sigma_{f_{t#bar{t}}^{"+tag+"}}",fttbarUnc_[tag]);
      model.fttbar_constrain[icat]         = new RooGaussian("ctr_fttbar_"+tag,"f_{t#bar{t}}^{"+tag+"} constrain",*model.fttbar[icat],*model.fttbar_mean_constrain[icat],*model.fttbar_sigma_constrain[icat]);
      if(addFttbarAsContrain) model.pdfConstrains.add(*model.fttbar_constrain[icat]);

      //single top fraction in the sample
      bool addFsingleTopAsContrain(true);
      if(fsingletop_[tag]==0 || fitEq ||  fitType_==FIT_R_AND_EB)  { model.fsingletop[icat] = new RooRealVar("fsingletop_"+tag,"f_{t}^{"+tag+"}",fsingletop_[tag]); addFsingleTopAsContrain=false; }
      else                                                           model.fsingletop[icat] = new RooRealVar("fsingletop_"+tag,"f_{t}^{"+tag+"}",fsingletop_[tag],fsingletopUnc_[tag]-3*fsingletopUnc_[tag],fsingletop_[tag]+3*fsingletopUnc_[tag]);
      model.fsingletop_mean_constrain[icat]  = new RooRealVar("meanfsingletop_"+tag,"#bar{f_{t}^{"+tag+"}}",fsingletop_[tag]);
      model.fsingletop_sigma_constrain[icat] = new RooRealVar("uncfsingletop_"+tag,"#sigma_{f_{t}^{"+tag+"}}",fsingletopUnc_[tag]);
      model.fsingletop_constrain[icat]       = new RooGaussian("ctr_fsingletop_"+tag,"f_{t}^{"+tag+"} constrain",*model.fsingletop[icat],*model.fsingletop_mean_constrain[icat],*model.fsingletop_sigma_constrain[icat]);
      if(addFsingleTopAsContrain) model.pdfConstrains.add(*model.fsingletop_constrain[icat]);
      
      //sample composition variables i.e. alpha and alpha_i
      float nPairsVal=4;
      if(icat==EE_3JETS || icat==MUMU_3JETS || icat==EMU_3JETS) nPairsVal=6;
      RooRealVar *npairs = new RooRealVar("npairs_"+tag,"npairs_"+tag,nPairsVal);
      model.alpha[icat]  = new RooFormulaVar("alpha_"+tag,  "(@0*@1)/(@2*(2.+@3))", RooArgSet(*npairs, *model.fcorrect[icat], *model.fttbar[icat], *model.fsingletop[icat]));
      model.alpha2[icat] = new RooFormulaVar("alpha2_"+tag, "pow(@0,2)*@1",           RooArgSet(*model.alpha[icat],*model.fttbar[icat]));
      model.alpha1[icat] = new RooFormulaVar("alpha1_"+tag, "(2*@0*(1-@0)+@0*@2)*@1", RooArgSet(*model.alpha[icat],*model.fttbar[icat],*model.fsingletop[icat]));
      model.alpha0[icat] = new RooFormulaVar("alpha0_"+tag, "1-@0-@1",                RooArgSet(*model.alpha2[icat],*model.alpha1[icat]));

      
      //define the base PDFs (kernel)
      //  For reference in the formulas this is the correspondence to the list of variables
      //                    R  Eb  Eq  alpha_2 alpha_1 alpha_0  acc1  acc05  acc0
      //                    @0 @1  @2  @3      @4      @5       @6    @7     @8
      RooArgList stdArgList(*model.r, *model.eb, *model.eq, *model.alpha2[icat], *model.alpha1[icat], *model.alpha0[icat],  *model.acc1, *model.acc05, *model.acc0);
      
      RooGenericPdf *pdf0btags2jets = new RooGenericPdf("kernelpdf0btags_"+tag,
							"@3*( pow(@0*(1-@1),2)*@6 + 2*@0*(1-@0)*(1-@1)*(1-@2)*@7 + pow((1-@0)*(1-@2),2)*@8) +"
							"@4*( pow(@0,2)*(1-@1)*(1-@2)*@6 + @0*(1-@0)*((1-@1)+(1-@2))*(1-@2)*@7 + pow((1-@0)*(1-@2),2)*@8)+"
							"@5*( pow(1-@2,2) )",		
							stdArgList);
  
      RooGenericPdf *pdf1btags2jets = new RooGenericPdf("kernelpdf1btags_"+tag,
							"@3*( 2*pow(@0,2)*(1-@1)*@1*@6 + 2*@0*(1-@0)*((1-@1)*@2+@1*(1-@2))*@7 + 2*pow((1-@0),2)*(1-@2)*@2*@8) +"
							"@4*( pow(@0,2)*(@1*(1-@2)+(1-@1)*@2)*@6 + @0*(1-@0)*(@1*(1-@2)+(1-@1)*@2+2*@2*(1-@2))*@7 + 2*pow((1-@0),2)*@2*(1-@2)*@8 ) +"
							"@5*( 2*@2*(1-@2) )",
							stdArgList);

      RooGenericPdf *pdf2btags2jets = new RooGenericPdf("kernelpdf2btags_"+tag,
							"@3*( pow(@0*@1,2)*@6 + 2*@0*(1-@0)*@1*@2*@7 + pow((1-@0)*@2,2)*@8 ) + "
							"@4*( pow(@0,2)*@1*@2*@6 + @0*(1-@0)*((@1+@2)*@2)*@7 + pow((1-@0)*@2,2)*@8 ) +"
							"@5*( pow(@2,2) )",
							stdArgList);
      
      //check if extension is needed
      if(icat==EE_3JETS || icat==MUMU_3JETS || icat==EMU_3JETS)
	{
	  //b-tags from extra jets are assumed to be mistags
	  RooGenericPdf *pdf0btagsextrajet = new RooGenericPdf("pdf0btagsextrajet_"+tag,"(1-@0)",RooArgList(*model.eq));
	  RooGenericPdf *pdf1btagsextrajet = new RooGenericPdf("pdf1btagsextrajet_"+tag,"@0",RooArgList(*model.eq));
	  
	  //
	  RooGenericPdf *pdf0btags3jets = new RooGenericPdf("pdf0btags3Jets_"+tag, "@0*@1",       RooArgList(*pdf0btags2jets, *pdf0btagsextrajet) );  
	  RooGenericPdf *pdf1btags3jets = new RooGenericPdf("pdf1btags3Jets_"+tag, "@0*@1+@2*@3", RooArgList(*pdf0btags2jets, *pdf1btagsextrajet, *pdf1btags2jets, *pdf0btagsextrajet) );
	  RooGenericPdf *pdf2btags3jets = new RooGenericPdf("pdf2btags3Jets_"+tag, "@0*@1+@2*@3", RooArgList(*pdf1btags2jets, *pdf1btagsextrajet, *pdf2btags2jets, *pdf0btagsextrajet) );
	  RooGenericPdf *pdf3btags3jets = new RooGenericPdf("pdf3btags3Jets_"+tag, "@0*@1",       RooArgList(*pdf2btags2jets, *pdf1btagsextrajet) );

	  model.pdfForCategory["n0btags_"+tag] = pdf0btags3jets;
	  model.pdfForCategory["n1btags_"+tag] = pdf1btags3jets;
	  model.pdfForCategory["n2btags_"+tag] = pdf2btags3jets;
	  model.pdfForCategory["n3btags_"+tag] = pdf3btags3jets;
	}
      else
	{
	  model.pdfForCategory["n0btags_"+tag] = pdf0btags2jets;
	  model.pdfForCategory["n1btags_"+tag] = pdf1btags2jets;
	  model.pdfForCategory["n2btags_"+tag] = pdf2btags2jets;
	}
    }


  // the model will be a simultaneous PDF
  model.pdf = new RooSimultaneous("hfcmodel","hfcmodel",*model.sample);
  for(std::map<TString, RooAbsPdf *>::iterator it = model.pdfForCategory.begin();   it != model.pdfForCategory.end(); it++)  model.pdf->addPdf( *(it->second), it->first );
      
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
void HFCMeasurement::fitHFCtoMeasurement(std::vector<TH1D *> &btagHistos, bool debug)
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
  runHFCFit(debug);
}


//
void HFCMeasurement::fitHFCtoEnsemble(top::EventSummaryHandler &evHandler,bool debug )
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
  runHFCFit(debug);
}


//
void HFCMeasurement::runHFCFit(bool debug)
{
  
  //build a categorized dataset
  int icat(0);
  RooDataSet *data = new RooDataSet("data","data",RooArgSet( *model.bmult, *model.bmultWgt, *model.sample), RooFit::WeightVar(*model.bmultWgt) );
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
	  model.bmultWgt->setVal(h->GetBinContent(ibin));
	  model.sample->setLabel(catname);
	  data->add( RooArgSet(*model.bmult,*model.bmultWgt,*model.sample), model.bmultWgt->getVal() );
	  cout << catname << " " << h->GetBinContent(ibin) << endl;
	}

      //set overall occupancy
      model.jetocc[icat]->setVal( h->Integral() );
    }

  data->Print();

  
  model.constrPdf->Print();
  // model.pdf->fitTo(*data,
  model.constrPdf->fitTo(*data,
			 Constrain(model.pdfConstrains),
			 RooFit::Save(),
			 RooFit::Hesse(kTRUE),
			 RooFit::Minos(kTRUE),
			 RooFit::SumW2Error(kTRUE),
			 RooFit::PrintLevel(-1),
			 RooFit::Verbose(kFALSE));


  //////// prepare the configuration for Feldman-Cousins
  RooWorkspace* w = new RooWorkspace();
  ModelConfig modelConfig("hfcfit",w);
  modelConfig.SetPdf(*model.pdf);
  modelConfig.SetParametersOfInterest(RooArgSet(*model.r)); //add the other ones?
  modelConfig.SetObservables(RooArgSet(*model.bmult,*model.bmultWgt,*model.sample));
  //  w->Print();
                                      
  // run FC
  RooStats::FeldmanCousins fc(*data,modelConfig);
  fc.SetTestSize(.05); // set size of test                                                                                                                                                                                                                                                            
  fc.UseAdaptiveSampling(true);
  fc.UseAdaptiveSampling(false);
  fc.FluctuateNumDataEntries(false); // number counting analysis: dataset always has 1 entry with N events observed                                                                                                                                                                     
  //  fc.SetNBins(100); // number of points to test per parameter
  fc.SetNBins(5); // number of points to test per parameter                                                                                                                                                                                                                                                              
  //get the interval
  try{
    PointSetInterval* interval = (PointSetInterval*)fc.GetInterval();
    
    std::cout << "is this point in the interval? " << interval->IsInInterval(model.pdfConstrains) << std::endl;
    std::cout << "interval is ["<< interval->LowerLimit(*model.r)  << ", "  << interval->UpperLimit(*model.r) << "]" << std::endl;
  }catch(std::string &e){
    cout << e << endl;
  }

  RooNLLVar *nll = (RooNLLVar *) model.constrPdf->createNLL(*data,Constrain(model.pdfConstrains),NumCPU(2));
  //RooNLLVar *nll = (RooNLLVar *) model.pdf->createNLL(*data,NumCPU(2));
  RooMinuit minuit(*nll); 
  minuit.setErrorLevel(0.5);     
  minuit.hesse();
  minuit.migrad();
  minuit.setErrorLevel(0.5);     
  minuit.hesse();
  minuit.minos();
  minuit.setErrorLevel(0.5);     
  RooFitResult *r=minuit.save();
  
  //
  // draw the result
  //
  if(debug)
    {
      //the dataset
      TCanvas *c = new TCanvas("HFCMeasurementData","HFCMeasurementData",600,600);
      c->SetWindowSize(600,600);
      RooPlot *frame = model.bmult->frame();
      data->plotOn(frame);
      frame->Draw();
      c->SaveAs("HFCMeasurementData.png");

      c = new TCanvas("combination","combination",600,600);
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
      else if (fitType_==FIT_R || fitType_==FIT_R_AND_EB) 
	{
	  frame = model.r->frame(Title("Likelihood"),Range(0.8,1.2)) ;    
	  frame->GetXaxis()->SetTitle("R");
	  char buf[100];
	  sprintf(buf,"R=%3.2f #pm %3.2f",model.r->getVal(),model.r->getError());
	  resLabel += buf;
	}
      nll->plotOn(frame,ShiftToZero(),Name("ll"));
      
      TLegend *leg = new TLegend(0.2,0.65,0.45,0.9,NULL,"brNDC");
      leg->AddEntry("ll","Combined","l");
     
      //debug
      cout << "****** SUMMARY RESULT ******" << endl
	   << resLabel << endl
       	   << "EDM = " << r->edm() << endl
       	   << "-log(L) at minimum = " << r->minNll() << endl   
       	   << "final value of floating parameters" << endl ;
      r->floatParsFinal().Print("s") ;
      cout << "****************************" << endl;
      
      //superimpose the other likelihoods now
//       int icat(0);
//       for(std::set<TString>::iterator cIt=categoryKeys_.begin(); cIt!=categoryKeys_.end(); cIt++,icat++)
// 	{
// 	  TString tag=*cIt;
// 	  TString lltitle("");
// 	  TString dilCategory("");
// 	  int lcolor(1), lstyle(1);
// 	  if(tag.Contains("ee"))   { dilCategory="ee";  lltitle="ee";      lcolor=kGreen+4;}
// 	  if(tag.Contains("emu"))  { dilCategory="emu"; lltitle="e#mu";    lcolor=kRed;    }
// 	  if(tag.Contains("mumu")) { dilCategory="mumu"; lltitle="#mu#mu"; lcolor=kOrange; }
// 	  TString jmult=tag;
// 	  jmult=jmult.ReplaceAll(dilCategory,"");
// 	  if(jmult=="3") lstyle=2;
// 	  lltitle += " (" + jmult + " jets)";
// 	  theLLsPerCategory[icat]->plotOn(frame, ShiftToZero(),  LineColor(lcolor), LineStyle(lstyle), LineWidth(2),  MoveToBack(), Name(tag+"ll") );
// 	  leg->AddEntry(tag+"ll",lltitle,"l");
// 	}

      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetYaxis()->SetTitle("-Log(L/L_{Max})");
      frame->GetYaxis()->SetTitleOffset(1);
      frame->GetYaxis()->SetRangeUser(0,5);
      frame->Draw();

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
//       if(fitType_ == FIT_EB_AND_EQ || fitType_==FIT_R_AND_EB)
// 	{
// 	  c = new TCanvas("contour","contour");
// 	  c->cd();
// 	  c->SetWindowSize(600,600);
// 	  c->SetGridx();
// 	  c->SetGridy();     
// 	  RooPlot *plot=0;
// 	  if(fitType_==FIT_R_AND_EB)
// 	    {
// 	      plot = minuit.contour(*model.r,*model.sfeb,1,2,3) ;
//               plot->SetTitle("Contour for 1s,2s,3s between r and eb") ;
// 	    }
// 	  else
// 	    {
// 	      plot = minuit.contour(*model.sfeb,*model.sfeq,1,2,3) ;
// 	      plot->SetTitle("Contour for 1s,2s,3s between eb and eq") ;
// 	    }	 
// 	  plot->Draw();
// 	  formatForCmsPublic(c,0,"",3);
// 	  c->Modified();
// 	  c->Update();
// 	  c->SaveAs("HFCContour.C");
// 	  c->SaveAs("HFCContour.pdf");
// 	  c->SaveAs("HFCContour.png");
// 	}
    }



  /*
  //build the likelihoods per category
  RooArgSet allLL, jetOccupancies;  
  TIterator *pdfIt     = model.constrPDFSet.createIterator();
  //TIterator *pdfIt     = model.pdfSet.createIterator();
  std::vector<RooNLLVar *> theLLsPerCategory;
  int icat(0);
  std::vector<RooRealVar *> roovars;
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
            
      //convert to a categorized dataset
      TH1 *h = controlHistos_.getHisto("btags_"+jmultkey,dilCategory);

      RooDataSet *ds= new RooDataSet("data_"+tag,"data_"+tag,RooArgSet( *model.bmult, *model.bmultWgt, *model.sample[icat]), RooFit::WeightVar(*model.bmultWgt) );
      cout << njets << endl;
      for(int ibin=1; ibin<= h->GetXaxis()->GetNbins(); ibin++)
	{
	  if(ibin==njets+2) break;
	  TString catname("n"); catname += (ibin-1); catname += "btags_" + tag;

	  model.bmult->setVal(ibin-1);
	  RooRealVar *wgt = new RooRealVar("cts"+catname,"cts"+catname, h->GetBinContent(ibin));
	  roovars.push_back(wgt);

	  model.bmultWgt->setVal(h->GetBinContent(ibin));
	}

      //convert to a binned dataset
      //      TH1 *h = controlHistos_.getHisto("btags_"+jmultkey,dilCategory);
      // model.bmult->setBins(maxJets_+1);
      // RooDataHist* ds = new RooDataHist("data"+tag,"data"+tag, RooArgSet(*model.bmult), h);

      alldata->append(*ds);
      ds->Print("v");
      if(icat==0) ds->plotOn(frame,MarkerStyle(20+icat));

      //set the overall normalization for this category
      model.jetocc[icat]->setVal(h->Integral());
      jetOccupancies.add(*(model.jetocc[icat]));

      //build the likelihood
      RooNLLVar *nll=0;
      RooProdPdf *modelconstr = dynamic_cast<RooProdPdf *>( pdfIt->Next() );
      modelconstr->Print("v");
      if(fitType_==FIT_EB || fitType_==FIT_R)
	{
	  //nll = (RooNLLVar *) modelconstr->createNLL(*ds,Constrain(RooArgSet(*model.fttbar[icat],*model.fsingletop[icat])),NumCPU(2));
	  //modelconstr->fitTo(*ds,Constrain(RooArgSet(*model.fttbar[icat],*model.fsingletop[icat])),Minos(),Save(kTRUE),PrintLevel(-1),Verbose(false));
	  modelconstr->fitTo(*ds,
			     RooFit::Save(),
			     RooFit::Hesse(kTRUE),
			     RooFit::Minos(kTRUE),
			     //			     RooFit::Extended(kFALSE),
			     // RooFit::SumW2Error(kTRUE),
			     RooFit::PrintLevel(-1),
			     RooFit::Verbose(kFALSE));
	  nll = (RooNLLVar *) modelconstr->createNLL(*ds,NumCPU(2));
	  
	  if(icat==0)
	    {
	      modelconstr->plotOn(frame,Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
	    }
	}
      else if(  fitType_==FIT_R_AND_EB)
	{
	  nll = (RooNLLVar *) modelconstr->createNLL(*ds,NumCPU(2));
	}
      else if(fitType_==FIT_EB_AND_EQ)//|| fitType_==FIT_R_AND_EQ)
	{
	  nll = (RooNLLVar *) modelconstr->createNLL(*ds,Constrain(RooArgSet(*model.sfeq)),NumCPU(2));
	  modelconstr->fitTo(*ds,Constrain(RooArgSet(*model.sfeq)),Minos(),Save(kTRUE),PrintLevel(-1),Verbose(false));
	}
      
      //save values fit and likelihood
      model.rFit[icat]=model.r->getVal();	
      model.rFitAsymmErrHi[icat]=model.r->getAsymErrorHi();
      model.rFitAsymmErrLo[icat]=model.r->getAsymErrorLo();

      model.ebFit[icat]=model.abseb->getVal()*model.sfeb->getVal();	
      model.ebFitAsymmErrHi[icat]=model.abseb->getVal()*model.sfeb->getAsymErrorHi();
      model.ebFitAsymmErrLo[icat]=model.abseb->getVal()*model.sfeb->getAsymErrorLo();

      model.eqFit[icat]=model.abseq->getVal()*model.sfeq->getVal();	
      model.eqFitAsymmErrHi[icat]=model.abseq->getVal()*model.sfeq->getAsymErrorHi();
      model.eqFitAsymmErrLo[icat]=model.abseq->getVal()*model.sfeq->getAsymErrorLo();
      
      RooMinuit minuit(*nll);
      minuit.setVerbose(false);
      minuit.migrad();
      minuit.hesse();
      TString itit="="+jmultkey + " jets";
      nll->SetTitle(itit);
      allLL.add(*nll); 
      theLLsPerCategory.push_back(nll);
    }

  //  alldata->plotOn(frame);
  frame->Draw();
  c->SaveAs("~/www/debug.png");
  
  
  // add up all the log likelihoods and fit them 
  // reinforce that the error is to be taken from  +/- 1/2
  RooAddition *combll = new RooAddition("combll","combll",allLL);      
  RooMinuit minuit(*combll); 
  minuit.setErrorLevel(0.5);     
  minuit.hesse();
  minuit.migrad();
  minuit.setErrorLevel(0.5);     
  minuit.hesse();
  minuit.minos();
  minuit.setErrorLevel(0.5);     
  
*/

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
