#include "LIP/Top/interface/HFCMeasurement.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuPhysicsEvent.h"
#include "CMGTools/HtoZZ2l2nu/src/JSONWrapper.cc"

#include "RooNumIntConfig.h"
#include "RooNLLVar.h"
#include "RooConstVar.h"
#include "RooExtendPdf.h"

#include "RooStats/ProfileLikelihoodTestStat.h"

#include  "TGraphErrors.h"
#include "TSystem.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;


HFCMeasurement::HFCMeasurement(int fitType,TString fitConfig, TString wpConfig)
{ 
  fitType_=fitType;
  switch(fitType_)
    {
    case FIT_EB:              fitTypeTitle_="#varepsilon_{b}";                    fitTypeName_="effb";                        break;
    case FIT_R_AND_EB:        fitTypeTitle_="R vs #varepsilon_{b}";               fitTypeName_="rvseffb";                     break;
    case FIT_R_AND_XSEC:      fitTypeTitle_="R vs #sigma";                        fitTypeName_="rvssigma";                    break;
    case FIT_EB_AND_XSEC:     fitTypeTitle_="#varepsilon_{b} vs #sigma";          fitTypeName_="effbvssigma";                 break;
    case FIT_EB_AND_EQ:       fitTypeTitle_="#varepsilon_{b} vs #varepsilon_{q}"; fitTypeName_="effbvseffq";                  break;
    default:                  fitTypeTitle_="R";                                  fitTypeName_="r";           fitType_=FIT_R; break;
    }

  mc_ = 0;
  ws_ = new RooWorkspace("w");
  parseFitConfig(fitConfig);
  parseFitConfig(wpConfig);
  initHFCModel();
}


//
void HFCMeasurement::parseFitConfig(TString url)
{
  //check if path is ok
  if(gSystem->AccessPathName(url)) return;

  char expBuf[500];  //free buffer	  
  RooArgSet nuis,constr;

  //read parameters from file
  JSONWrapper::Object jsonF(url.Data(), true);
  
  std::vector<JSONWrapper::Object> params=jsonF.daughters();
  for(size_t iparam=0; iparam<params.size(); iparam++)
    {
      std::vector<JSONWrapper::Object> &dau = params[iparam].daughters();
      for(size_t icat=0; icat<dau.size(); icat++)
	{
	  JSONWrapper::Object &descript=dau[icat];
	  string param=jsonF.key[iparam]+"_"+params[iparam].key[icat];
	  sampleCats_.insert(params[iparam].key[icat]);

	  std::vector<std::string> uncs = descript.key;

	  //for each parameter instantiate a formula, per category of the type: x\prod_i(1+theta_i)
	  //where theta_i are the nuisances which will be varied with a given prior
	  RooRealVar *cenVal=ws_->var((param+"_cen").c_str());
	  if(cenVal==0)
	    {
	      string formula("@0");
	      RooArgList varsInFormula;

	      //central value
	      sprintf(expBuf,"%s_cen[%f]",param.c_str(),descript["val"].toDouble());
	      cenVal=(RooRealVar *)ws_->factory(expBuf);
	      varsInFormula.add(*cenVal);
	  
	      //modifiers
	      int uncCntr(0);
	      for(size_t iunc=0; iunc<uncs.size(); iunc++)
		{
		  if(uncs[iunc]=="val") continue;
		  
		  //add new constraint if required
		  RooRealVar *nuisVar=ws_->var(uncs[iunc].c_str());
		  RooRealVar *modVar=ws_->var((uncs[iunc]+"_sigma").c_str());
		  if(nuisVar==0)
		    {
		      sprintf(expBuf,"Gaussian::%s_constr(%s[0,-1,1],0.0,1.0)",uncs[iunc].c_str(), uncs[iunc].c_str());
		      RooGaussian *nuisConstr=(RooGaussian *)ws_->factory(expBuf);
		      constr.add( *nuisConstr );

		      nuisVar = ws_->var( uncs[iunc].c_str() );
		      nuis.add( *nuisVar );

		      sprintf(expBuf,"%s_sigma[%f]",uncs[iunc].c_str(),descript[uncs[iunc].c_str()].toDouble());
		      modVar=(RooRealVar *)ws_->factory(expBuf);
		    }
		  
		  //add to formula
		  sprintf(expBuf,"*(1+@%d*@%d)",uncCntr+1,uncCntr+2);
		  formula+=expBuf;
		  uncCntr+=2;
		  varsInFormula.add( *modVar );
		  varsInFormula.add( *nuisVar );
		}
	      
	      //final formula
	      RooFormulaVar formulaVar(param.c_str(),formula.c_str(),varsInFormula);
	      ws_->import(formulaVar);
	    }

	  else
	    {
	      //if parameter has already been parsed, update the values only
	      cenVal->setVal(descript["val"].toDouble());
	      for(size_t iunc=0; iunc<uncs.size(); iunc++)
		{
		  if(uncs[iunc]=="val") continue;
		  
		  //add new constraint if required
		  RooRealVar *nuisVar=ws_->var((uncs[iunc]+"_sigma").c_str());
		  if(nuisVar==0) cout << "[Warning] failed to update uncertainty for: " << uncs[iunc] << " ...neglecting this source" << endl; 
		  else nuisVar->setVal(descript[ uncs[iunc].c_str() ].toDouble());
		}
	    }
	}
    }
  
  if(ws_->set("nuisances")!=0) nuis.add( *(ws_->set("nuisances")) );
  ws_->defineSet("nuisances",nuis);

  if(ws_->set("constr")!=0) constr.add( *(ws_->set("constr")) );
  ws_->defineSet("constr",constr);
 
}

//
void HFCMeasurement::initHFCModel()
{
  if(ws_==0) return;
  if(ws_->var("r")!=0) return;

  char expBuf[200];

  //observables
  ws_->factory("bmultobs[0.,99999999999.]");
  sprintf(expBuf,"bmult[0,%f]",float(maxJets_+1));
  ws_->factory(expBuf);
  ws_->var("bmult")->setBins(maxJets_+1);
  RooCategory sample("sample","sample");

  //
  //parameters of interest
  //
  RooArgSet poi;

  float minR(0.0), maxR(2.0);
  if(fitType_==FIT_R_AND_EB)                                                    { minR=0.96; maxR=1.06; }
  if(fitType_== FIT_EB || fitType_==FIT_EB_AND_XSEC || fitType_==FIT_EB_AND_EQ) { minR=1.0; maxR=1.0;   }
  sprintf(expBuf,"r[1,%f,%f]",minR,maxR);
  RooRealVar *r = (RooRealVar *)ws_->factory(expBuf);
  if(minR!=maxR) poi.add( *ws_->var("r") );
  
  float minSFb(1.0),maxSFb(1.0);
  if(fitType_==FIT_R_AND_EB)                                                    { minSFb=0.9; maxSFb=1.1; }
  if(fitType_== FIT_EB || fitType_==FIT_EB_AND_XSEC || fitType_==FIT_EB_AND_EQ) { minSFb=0.7; maxSFb=1.3; }
  sprintf(expBuf,"sfeb[1,%f,%f]",minSFb,maxSFb);
  RooRealVar *sfeb=(RooRealVar *)ws_->factory(expBuf);
  if(minSFb!=maxSFb) poi.add( *ws_->var("sfeb") );

  float minSFq(1.0),maxSFq(1.0);
  if(fitType_==FIT_EB_AND_EQ) { minSFq=0.7; maxSFq=1.8; }
  sprintf(expBuf,"sfeq[1,%f,%f]",minSFq,maxSFq);
  RooRealVar *sfeq=(RooRealVar *)ws_->factory(expBuf);
  if(minSFq!=maxSFq) poi.add( *ws_->var("sfeq") );

  float minXsec(1.0),maxXsec(1.0);
  if(fitType_== FIT_R_AND_XSEC || fitType_==FIT_EB_AND_XSEC) { minXsec=0.8; maxXsec=1.2; }
  sprintf(expBuf,"xsec[1,%f,%f]",minXsec,maxXsec);
  RooRealVar *xsec =  (RooRealVar *)ws_->factory(expBuf);
  if(minXsec!=maxXsec) poi.add( *ws_->var("xsec") );


  //
  // model
  //
  //acceptance corrections
  RooRealVar *acc1  = (RooRealVar *) ws_->factory("a1[1.0]");
  RooRealVar *acc05 = (RooRealVar *) ws_->factory("a05[1.0]");
  RooRealVar *acc0  = (RooRealVar *) ws_->factory("a0[1.0]");
  
  std::map<string,string> basePdfsPerCat; //auxiliary
  for(std::set<string>::iterator cIt=sampleCats_.begin(); cIt!=sampleCats_.end(); cIt++)
    {
      string tag=*cIt;
      int njets(2);
      if(tag.find("eq3")!=string::npos) njets=3;
      if(tag.find("eq4")!=string::npos) njets=4;

      //sample composition base parameters
      RooAbsReal *fcorrect = ws_->function(("fcorrect_"+tag).c_str());
      RooAbsReal *ktt      = ws_->function(("ktt_"+tag).c_str());
      RooAbsReal *kt       = ws_->function(("kt_"+tag).c_str());
      RooAbsReal *mu       = ws_->function(("mu_"+tag).c_str());
      if(fcorrect==0 || ktt==0 || kt==0 || mu==0)
	{
	  cout << "[Warning] unable to instantiate model for " << tag << " check the configuration file for sample composition parameters" << endl;
	  continue;
	}
      
      sprintf(expBuf,"FormulaVar::fttbar_%s('max(min(@0*@1*@2,1.0),0.0)',{%s,%s,%s})",tag.c_str(),ktt->GetName(),mu->GetName(),xsec->GetName());
      RooAbsReal * fttbar = (RooAbsReal *)ws_->factory(expBuf);

      //probability to select/reconstruct N(t->Wq)
      sprintf(expBuf,"FormulaVar::alpha_%s('max(min((2*%d*@0)/(@1*(2.+@2)),1.0),0.0)',{%s,%s,%s})",tag.c_str(),njets,fcorrect->GetName(),fttbar->GetName(),kt->GetName());
      RooAbsReal *alpha = (RooAbsReal *)ws_->factory(expBuf);

      //sample composition in terms of N(t->Wq)
      sprintf(expBuf,"FormulaVar::alpha2_%s('max(min(pow(@0,2)*@1,1.0),0.0)',{%s,%s})",tag.c_str(),alpha->GetName(),fttbar->GetName());
      RooAbsReal *alpha2=(RooAbsReal *)ws_->factory(expBuf);
      sprintf(expBuf,"FormulaVar::alpha1_%s('max(min((2*@0*(1-@0)+@0*@2)*@1,1.0),0.0)',{%s,%s,%s})",tag.c_str(),alpha->GetName(),fttbar->GetName(),kt->GetName());
      RooAbsReal *alpha1=(RooAbsReal *)ws_->factory(expBuf);
      sprintf(expBuf,"FormulaVar::alpha0_%s('max(min(1-@0-@1,1.0),0.0)',{%s,%s})",tag.c_str(),alpha2->GetName(),alpha1->GetName());
      RooAbsReal *alpha0=(RooAbsReal *)ws_->factory(expBuf);

      //b-tagging base parameters
      RooAbsReal *abseb       = ws_->function(("abseb_"+tag).c_str());
      RooAbsReal *abseq       = ws_->function(("abseq_"+tag).c_str());
      if(abseb==0 || abseq==0)
	{
	  cout << "[Warning] unable to instantiate model for " << tag << " check the configuration file for b-tagging parameters" << endl;
	  continue;
	}
      
      //modify locally the b-tag/mistag rates with correction factors
      sprintf(expBuf,"FormulaVar::eb_%s('max(min(@0*@1,1.0),0.0)',{%s,%s})",tag.c_str(),abseb->GetName(),sfeb->GetName());
      RooAbsReal *eb=(RooAbsReal *)ws_->factory(expBuf);

      sprintf(expBuf,"FormulaVar::eq_%s('max(min(@0*@1,1.0),0.0)',{%s,%s})",tag.c_str(),abseq->GetName(),sfeq->GetName());
      RooAbsReal *eq=(RooAbsReal *)ws_->factory(expBuf);

      //define the base PDFs (kernel)
      RooArgList stdArg;
      stdArg.add(*r);             // @0  R
      stdArg.add(*eb);            // @1  eff_b
      stdArg.add(*eq);            // @2  eff_q
      stdArg.add(*alpha2);        // @3  alpha_2
      stdArg.add(*alpha1);        // @4  alpha_1
      stdArg.add(*alpha0);        // @5  alpha_0 
      stdArg.add(*acc1);          // @6  acc_1
      stdArg.add(*acc05);         // @7  acc05
      stdArg.add(*acc0);          // @8  acc0
      
      RooGenericPdf *pdf2j0t = new RooGenericPdf(("pdf2j0t_"+tag).c_str(),
						 "@3*( pow(@0*(1-@1),2)*@6 + 2*@0*(1-@0)*(1-@1)*(1-@2)*@7 + pow((1-@0)*(1-@2),2)*@8) +"
						 "@4*( pow(@0,2)*(1-@1)*(1-@2)*@6 + @0*(1-@0)*((1-@1)+(1-@2))*(1-@2)*@7 + pow((1-@0)*(1-@2),2)*@8)+"
						 "@5*( pow(1-@2,2) )",
						 stdArg);
      
      RooGenericPdf *pdf2j1t = new RooGenericPdf(("pdf2j1t_"+tag).c_str(),
						 "@3*( 2*pow(@0,2)*(1-@1)*@1*@6 + 2*@0*(1-@0)*((1-@1)*@2+@1*(1-@2))*@7 + 2*pow((1-@0),2)*(1-@2)*@2*@8) +"
						 "@4*( pow(@0,2)*(@1*(1-@2)+(1-@1)*@2)*@6 + @0*(1-@0)*(@1*(1-@2)+(1-@1)*@2+2*@2*(1-@2))*@7 + 2*pow((1-@0),2)*@2*(1-@2)*@8 ) +"
						 "@5*( 2*@2*(1-@2) )",
						 stdArg);

      RooGenericPdf *pdf2j2t = new RooGenericPdf(("pdf2j2t_"+tag).c_str(),
						 "@3*( pow(@0*@1,2)*@6 + 2*@0*(1-@0)*@1*@2*@7 + pow((1-@0)*@2,2)*@8 ) + "
						 "@4*( pow(@0,2)*@1*@2*@6 + @0*(1-@0)*((@1+@2)*@2)*@7 + pow((1-@0)*@2,2)*@8 ) +"
						 "@5*( pow(@2,2) )",
						 stdArg);
      
      //b-tags from extra jets are assumed to be effective mistags (i.e. may include b-jets from non-top decays)
      if(njets>4)
	{
	  cout << "[Warning] failed to instantiate pdf for " << tag << endl;
	  break;
	} 
      else if(njets==2)
	{
	  ws_->import(*pdf2j0t);
	  ws_->import(*pdf2j1t);
	  ws_->import(*pdf2j2t);
	}
      else
	{
	  RooGenericPdf *pdf1ej0t = new RooGenericPdf(("pdf1ej0t_"+tag).c_str(), "(1-@0)",      RooArgList(*eq));
	  RooGenericPdf *pdf1ej1t = new RooGenericPdf(("pdf1ej1t_"+tag).c_str(), "@0",          RooArgList(*eq));
	  RooGenericPdf *pdf3j0t  = new RooGenericPdf(("pdf3j0t_"+tag).c_str(),  "@0*@1",       RooArgList(*pdf2j0t, *pdf1ej0t) );  
	  RooGenericPdf *pdf3j1t  = new RooGenericPdf(("pdf3j1t_"+tag).c_str(),  "@0*@1+@2*@3", RooArgList(*pdf2j0t, *pdf1ej1t, *pdf2j1t, *pdf1ej0t) );
	  RooGenericPdf *pdf3j2t  = new RooGenericPdf(("pdf3j2t_"+tag).c_str(),  "@0*@1+@2*@3", RooArgList(*pdf2j1t, *pdf1ej1t, *pdf2j2t, *pdf1ej0t) );
	  RooGenericPdf *pdf3j3t  = new RooGenericPdf(("pdf3j3t_"+tag).c_str(),  "@0*@1",       RooArgList(*pdf2j2t, *pdf1ej1t) );
	  if(njets==3)
	    {
	      ws_->import(*pdf3j0t);
	      ws_->import(*pdf3j1t);
	      ws_->import(*pdf3j2t);
	      ws_->import(*pdf3j3t);
	    }
	  else
	    {
	      RooGenericPdf *pdf4j0t  = new RooGenericPdf(("pdf4j0t_"+tag).c_str(),  "@0*@1",       RooArgList(*pdf3j0t, *pdf1ej0t) );  
	      RooGenericPdf *pdf4j1t  = new RooGenericPdf(("pdf4j1t_"+tag).c_str(),  "@0*@1+@2*@3", RooArgList(*pdf3j0t, *pdf1ej1t, *pdf3j1t, *pdf1ej0t) );
	      RooGenericPdf *pdf4j2t  = new RooGenericPdf(("pdf4j2t_"+tag).c_str(),  "@0*@1+@2*@3", RooArgList(*pdf3j1t, *pdf1ej1t, *pdf3j2t, *pdf1ej0t) );
	      RooGenericPdf *pdf4j3t  = new RooGenericPdf(("pdf4j3t_"+tag).c_str(),  "@0*@1",       RooArgList(*pdf3j2t, *pdf1ej1t) );
	      RooGenericPdf *pdf4j4t  = new RooGenericPdf(("pdf4j4t_"+tag).c_str(),  "@0*@1",       RooArgList(*pdf3j2t, *pdf1ej1t) );
	      ws_->import(*pdf4j0t);
	      ws_->import(*pdf4j1t);
	      ws_->import(*pdf4j2t);
	      ws_->import(*pdf4j3t);
	      ws_->import(*pdf4j4t);
	    }
	}
      
      //each n-jets n-b-tags is a category
      for(int ijet=0; ijet<=njets; ijet++)
	{
	  sprintf(expBuf,"%dj%dt_%s",njets,ijet,tag.c_str());
	  sample.defineType((string("n")+expBuf).c_str());
	  basePdfsPerCat[string("n")+expBuf]=string("pdf")+expBuf;
	}
    }


  //the model will be a (constrained +) extended + simultaneous PDF
  RooSimultaneous *basemodel = new RooSimultaneous("basemodel","basemodel",sample);
  for(std::map<string, string>::iterator it=basePdfsPerCat.begin(); it!=basePdfsPerCat.end(); it++)
    {
      RooAbsPdf *pdf=ws_->pdf(it->second.c_str());
      if(pdf==0) continue;
      basemodel->addPdf( *pdf, it->first.c_str() );
    }
  RooRealVar *yields = (RooRealVar *)ws_->factory("yields[1,0,999999999.]");
  RooAbsPdf *model=0;
  if(ws_->set("constr")!=0)
    {
      RooExtendPdf *extmodel = new RooExtendPdf("extmodel","extmodel",*basemodel,*yields);
      RooArgSet modelFactors(*extmodel); 
      modelFactors.add( *(ws_->set("constr")) );
      model = (RooAbsPdf *)(new RooProdPdf("model","model",modelFactors));
    }
  else
    {
      model = (RooAbsPdf *)(new RooExtendPdf("model","model",*basemodel,*yields));
    }
  ws_->import(*model);

  //finalize workspace
  ws_->defineSet("observables","bmultobs,bmult,sample");
  ws_->defineSet("poi",poi);  
  RooArgSet *nullParams = (RooArgSet *)ws_->allVars().snapshot();
  ws_->saveSnapshot("default",*nullParams,kTRUE);
  ws_->Print("v");

  //instantiate a model configurator

  mc_ = new ModelConfig("mc",ws_);
  mc_->SetPdf(*model);
  mc_->SetParametersOfInterest(*(ws_->set("poi")));
  mc_->SetObservables(*(ws_->set("observables")));
  mc_->SetNuisanceParameters(*(ws_->set("nuisances")));
}


//
void HFCMeasurement::resetModelValues()
{
  if(ws_==0) return;
  ws_->loadSnapshot("default");
}

//
void HFCMeasurement::fitHFCtoMeasurement(std::vector<TH1D *> &btagHistos, int runMode, bool debug)
{
  /*
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

  */
}


//
void HFCMeasurement::fitHFCtoEnsemble(ZZ2l2nuSummaryHandler  &evHandler, int runMode, bool debug )
{
  /*
  if(evHandler.getEntries()==0) return;

  //restart all over again
  initHFCModel();
  resetModelValues();
  resetHistograms();

  nMeasurements_++;
  for(int i=0; i<evHandler.getEntries(); i++)
    {
      evHandler.getEntry(i);
      ZZ2l2nuSummary_t &ev = evHandler.getEvent();
      PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      
      //check dilepton
      LorentzVector dilepton = phys.leptons[0]+phys.leptons[1];
      float dilcharge=(phys.leptons[0].id/fabs(phys.leptons[0].id)) *(phys.leptons[1].id/fabs(phys.leptons[1].id));
      if((ev.cat==EE || ev.cat==MUMU) && (fabs(dilepton.mass()-91)<15 || phys.met[0].pt()<30)) continue;
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
	  else if(btagAlgo_.Contains("CSV") )    btag = phys.jets[ijet].btag2;
	  else if(btagAlgo_.Contains("JP") )     btag = phys.jets[ijet].btag3;
	  else if(btagAlgo_.Contains("TCHP") )   btag = phys.jets[ijet].btag4;
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

  */
}


//
void HFCMeasurement::runHFCFit(int runMode, bool debug)
{
  
  /*

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
  std::vector<RooAbsReal *> profLikelihoods;
  std::vector<TString> lltitles; 
  TGraphAsymmErrors *exclusiveRgr=0;
  TGraphAsymmErrors *incRgr=0;
  std::vector<TPaveText *> exclusiveFitNames;
  if(runMode==FITEXCLUSIVECATEGORIES)
    {
      exclusiveRgr = new TGraphAsymmErrors;
      incRgr=new TGraphAsymmErrors;
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
	  RooNLLVar *nll = (RooNLLVar *) model.constrPdf->createNLL(*dataslice,Constrain(model.pdfConstrains),NumCPU(2));
	  likelihoods.push_back( nll );
	  profLikelihoods.push_back( nll->createProfile( (fitType_==FIT_EB || fitType_==FIT_EB_AND_EQ) ? *model.sfeb : *model.r ) ); 
	  TString tit=tag2.ReplaceAll("2","");
	  tit=tit.ReplaceAll("mu","#mu");
	  lltitles.push_back(tit);

	  model.rCombFit[icat-1]=model.r->getVal();   model.rCombFitAsymmErrLo[icat-1]=model.r->getAsymErrorLo(); model.rCombFitAsymmErrHi[icat-1]=model.r->getAsymErrorHi();
	  model.rCombFit[icat]=model.r->getVal();     model.rCombFitAsymmErrLo[icat]=model.r->getAsymErrorLo();   model.rCombFitAsymmErrHi[icat]=model.r->getAsymErrorHi();
	}
  
      //inclusive fit to two jets bin
      TString cut2inc("sample==sample::n0btags_ee2 || sample==sample::n1btags_ee2 || sample==sample::n2btags_ee2 ||"
		      "sample==sample::n0btags_mumu2 || sample==sample::n1btags_mumu2 || sample==sample::n2btags_mumu2 ||"
		      "sample==sample::n0btags_emu2 || sample==sample::n1btags_emu2 || sample==sample::n2btags_emu2");
      RooDataSet *dataslice = (RooDataSet *) data->reduce(cut2inc);
      model.constrPdf->fitTo(*dataslice,
			     Constrain(model.pdfConstrains),
			     RooFit::Save(),
			     RooFit::Hesse(kTRUE),
			     RooFit::Minos(kTRUE),
			     RooFit::SumW2Error(kTRUE),
			     RooFit::PrintLevel(-1),
			     RooFit::Verbose(kFALSE));
      incRgr->SetPoint(0,1,model.r->getVal());
      incRgr->SetPointError(0,0,0,fabs(model.r->getAsymErrorLo()),fabs(model.r->getAsymErrorHi()));

      //inclusive fit to two jets bin
      TString cut3inc("sample==sample::n0btags_ee3   || sample==sample::n1btags_ee3   || sample==sample::n2btags_ee3   || sample==sample::n3btags_ee3 ||"
		      "sample==sample::n0btags_mumu3 || sample==sample::n1btags_mumu3 || sample==sample::n2btags_mumu3 || sample==sample::n3btags_mumu3 ||"
		      "sample==sample::n0btags_emu3  || sample==sample::n1btags_emu3  || sample==sample::n2btags_emu3  || sample==sample::n3btags_mumu3");
      dataslice = (RooDataSet *) data->reduce(cut3inc);
      model.constrPdf->fitTo(*dataslice,
			     Constrain(model.pdfConstrains),
			     RooFit::Save(),
			     RooFit::Hesse(kTRUE),
			     RooFit::Minos(kTRUE),
			     RooFit::SumW2Error(kTRUE),
			     RooFit::PrintLevel(-1),
			     RooFit::Verbose(kFALSE));
      incRgr->SetPoint(1,4,model.r->getVal());
      incRgr->SetPointError(1,0,0,fabs(model.r->getAsymErrorLo()),fabs(model.r->getAsymErrorHi()));
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
  RooAbsReal *pll= nll->createProfile( (fitType_==FIT_EB || fitType_==FIT_EB_AND_EQ) ? *model.sfeb : *model.r ); 

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
      
      //configure FC
      if(fc_==0) fc_ = new RooStats::FeldmanCousins(*data,*model.modelConfig);
      RooDataSet poiToTest("poitotest","poitotest",RooArgSet(*model.r));
      double rmin=0.;//max(double(0.),double(model.rFitResult+model.rFitResultAsymmErrLo*5));
      double rmax=1.0;
      for(double r=rmin; r<=rmax; r+= (rmax-rmin)/20.)
	{
	  model.r->setVal(r);
	  poiToTest.add(RooArgSet(*model.r));
	}
      fc_->SetPOIPointsToTest(poiToTest);
      fc_->SetTestSize(.05);               // set size of test 
      fc_->FluctuateNumDataEntries(false); // number counting analysis: dataset always has 1 entry with N events observed
      fc_->SetConfidenceLevel(0.95);       // 95% interval
      fc_->UseAdaptiveSampling(true);      // speed it up a bit
      fc_->CreateConfBelt(true);
      fc_->SaveBeltToFile(true);
      
      //set as one sided
      //ToyMCSampler*  toymcsampler = (ToyMCSampler*) fc_->GetTestStatSampler(); 
      //ProfileLikelihoodTestStat* testStat = dynamic_cast<ProfileLikelihoodTestStat*>(toymcsampler->GetTestStatistic());
      //testStat->SetOneSided(true);

// 	cout << "Starting Feldman-Cousins computation: you can go and take a loooooong coffee " << endl;
// 	PointSetInterval* interval = (PointSetInterval*)fc_->GetInterval();
// 	ConfidenceBelt* belt = fc_->GetConfidenceBelt();
// 	model.rFitLowerLimit=interval->LowerLimit(*model.r);
// 	model.rFitUpperLimit=interval->UpperLimit(*model.r);

      
      //get the thresholds
      if(debug)
	{
	  //save workspace to file
	  model.ws->import(*model.modelConfig);
	  model.ws->import(*data);
	  model.ws->writeToFile("HeavyFlavorWS.root");

	  //display in canvas
// 	    TCanvas *c = new TCanvas("HFCIntervalForFeldmanCousins","HFCIntervalForFeldmanCousins",600,600);
// 	    c->SetWindowSize(600,600);     
// 	    TGraph *frame = new TGraph;
// 	    RooDataSet* parameterScan = (RooDataSet*) fc_->GetPointsToScan();
// 	    for(Int_t i=0; i<parameterScan->numEntries(); ++i)
// 	    {
// 	    RooArgSet *tmpPoint = (RooArgSet*) parameterScan->get(i)->clone("temp");
// 	    double arMax = belt->GetAcceptanceRegionMax(*tmpPoint);
// 	    double poiVal = tmpPoint->getRealValue(model.r->GetName()) ;	  
// 	    frame->SetPoint(i,poiVal,arMax);
// 	    }
// 	    frame->Draw("ap");
// 	    frame->GetXaxis()->SetTitle("R");
// 	    frame->GetYaxis()->SetTitle("Acceptance region max");
	  
// 	    char buf[100];
// 	    sprintf(buf,"%3.2f < R < %3.2f",model.rFitLowerLimit,model.rFitUpperLimit);
// 	    TPaveText *pave = new TPaveText(0.15,0.96,0.41,0.99,"NDC");
// 	    pave->SetBorderSize(0);
// 	    pave->SetFillStyle(0);
// 	    pave->SetTextAlign(12);
// 	    pave->SetTextFont(42);
// 	    pave->AddText("CMS preliminary");
// 	    pave->Draw();
	    
// 	    pave = new TPaveText(0.4,0.96,0.94,0.99,"NDC");
// 	    pave->SetFillStyle(0);
// 	    pave->SetBorderSize(0);
// 	    pave->SetTextAlign(32);
// 	    pave->SetTextFont(42);
// 	    pave->AddText(buf);
// 	    pave->Draw();
	    
// 	    c->SaveAs("HFCIntervalForFeldmanCousins.png");
// 	    c->SaveAs("HFCIntervalForFeldmanCousins.pdf");
// 	    c->SaveAs("HFCIntervalForFeldmanCousins.C");
//	    }
	
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
      incRgr->SetFillStyle(0);
      incRgr->SetLineColor(kBlue);
      incRgr->SetMarkerStyle(24);
      incRgr->SetMarkerColor(kBlue);
      incRgr->Draw("p");
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
      // PROFILE LIKELIHOODS
      //
      c = new TCanvas("HFCMeasurementProfLikelihood","HFCMeasurementProfLikelihood",600,600);
      c->SetWindowSize(600,600);
      
      //plot the likelihood in the main frame
      if(fitType_==FIT_EB || fitType_==FIT_EB_AND_EQ) 
	{
	  frame = model.sfeb->frame(Title("Likelihood"),Range(0.5,1.3));
	  frame->GetXaxis()->SetTitle("SF #varepsilon_{b}");
	}
      else if (fitType_==FIT_R || fitType_==FIT_R_AND_EB || fitType_==FIT_R_CONSTRAINED) 
	{
	  frame = model.r->frame(Title("Likelihood"),Range(0.8,1.2)) ;    
	  frame->GetXaxis()->SetTitle("R");
	}
      
      pll->plotOn(frame,ShiftToZero(),Name("pll"),FillStyle(0),LineColor(kBlue),LineWidth(2));
      for(size_t ill=0; ill<profLikelihoods.size(); ill++)
	{
	  if(profLikelihoods[ill]==0) continue;
	  TString pllname("pll"); pllname +=ill;
	  profLikelihoods[ill]->plotOn(frame, ShiftToZero(), LineColor(llColors[ill]), LineWidth(2),Name(pllname));
	}
      
      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetYaxis()->SetTitle("-Log(L/L_{0})");
      frame->GetYaxis()->SetTitleOffset(1);
      frame->GetYaxis()->SetRangeUser(0,5);
      frame->Draw();
      
      pave = new TPaveText(0.15,0.96,0.41,0.99,"NDC");
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

      leg = new TLegend(0.2,0.65,0.45,0.9,NULL,"brNDC");
      leg->AddEntry("pll","Combined","l");     
      for(size_t ill=0; ill<profLikelihoods.size(); ill++)
	{
	  TString pllname("pll"); pllname +=ill;
	  leg->AddEntry(pllname,lltitles[ill],"l");
	}

      formatForCmsPublic(c,leg,"",5);
      leg->SetFillColor(0);
      leg->SetFillStyle(3001);
      leg->Draw();
      
      //draw the data and the model sum in a sub-pad
      npad = new TPad("pllpad","pll", 0.6, 0.6, 0.9, 0.9);
      npad->Draw();
      npad->cd();
      frame = model.bmult->frame();

      data->plotOn(frame,DrawOption("pz"));
      frame->Draw();
      h->DrawClone("histsame");
      c->Modified();
      c->Update();
      c->SaveAs("HFCProfLikelihood.png");
      c->SaveAs("HFCProfLikelihood.pdf");
      c->SaveAs("HFCProfLikelihood.C");



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
      std::vector<TGraph *> incPDF;
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
	  for(size_t ibtags=0;ibtags<=size_t(maxJets_); ibtags++)
	    {
	      TString catname("n"); catname += ibtags; catname += "btags_" + tag;
	      if(model.pdfForCategory.find(catname)==model.pdfForCategory.end()) continue;
	      
	      TGraph * modelgr = new TGraph;
	      modelgr->SetName(catname);
	      modelgr->SetLineColor(modelColors[ibtags]);
	      modelgr->SetFillStyle(0);
	      modelgr->SetLineWidth(2);

	      bool newIncPDF(false);
	      if(incPDF.size()<=ibtags) 
		{
		  newIncPDF=true;
		  incPDF.push_back( (TGraph *) modelgr->Clone() );
		  TString ibtagsStr("btags_"); ibtagsStr+=ibtags;
		  incPDF[ibtags]->SetName("inc"+ibtagsStr);
		}
	      
	      //loop over R
	      for(float ir=0; ir<=1.01; ir+=0.01)
		{
		  model.r->setVal(ir);
		  int ipt=modelgr->GetN();
		  Double_t ival=model.pdfForCategory[catname]->getVal();
		  Double_t iocc=model.jetocc[icat]->getVal();
		  modelgr->SetPoint(ipt,ir,ival);
		  if(newIncPDF) incPDF[ibtags]->SetPoint(ipt,ir,ival*iocc);
		  else
		    {
		      Double_t ix(0),iy(0);
		      incPDF[ibtags]->GetPoint(ipt,ix,iy);
		      incPDF[ibtags]->SetPoint(ipt,ir,iy+ival*iocc);
		    }
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

      //inclusive model
      c = new TCanvas("HeavyFlavorModelInclusive","HeavyFlavorContourInclusive");
      c->SetCanvasSize(600,600);
      c->SetWindowSize(600,600);
      leg=new TLegend(0.2,0.65,0.45,0.9,NULL,"brNDC");      
      for(size_t ibtags=0; ibtags<incPDF.size(); ibtags++) 
	{
	  incPDF[ibtags]->Draw(ibtags==0 ? "al" :"l" );
	  
	  TString tit("="); tit+=ibtags; tit+=" b-tags";
	  leg->AddEntry(incPDF[ibtags]->GetName(),tit,"l");
	  
          incPDF[ibtags]->GetXaxis()->SetTitle("R=B(t#rightarrow Wb)/B(t#rightarrow Wq)");
	  incPDF[ibtags]->GetYaxis()->SetTitle("Probability");
	}
      leg->Draw();
      formatForCmsPublic(c,leg,"CMS preliminary",3);
      c->cd();
      c->Modified();
      c->Update();
      c->SaveAs("HeavyFlavorModelInclusive.C");
      c->SaveAs("HeavyFlavorModelInclusive.pdf");
      c->SaveAs("HeavyFlavorModelInclusive.png");
    }

  delete data;
*/

}



//
void HFCMeasurement::printConfiguration(std::ostream &os)
{

  /*

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


  */
}
