#include "LIP/Top/interface/MassMeasurement.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include <memory>
#include "TCut.h"
#include "RooHist.h"

using namespace std;
using namespace RooFit;
using namespace top;

//
MassMeasurement::MassMeasurement(TString parfileURL,TString tag)
{
  fitPars_ = ParseParametersFrom(parfileURL);
  tag_=tag;
  tag_.ReplaceAll(".","");
  InitModel();
}


//
EnsembleMeasurement_t MassMeasurement::DoMassFit(top::EventSummaryHandler &evHandler, bool debug)
{
  EnsembleMeasurement_t em;
  em.status=false;
  em.nEvents = 0;
  TTree *evTree=evHandler.getTree();
  if(evTree==0) return em;
  if(evTree->GetEntriesFast()==0) return em;

  //configure the pre-selection
  float minTopMass=fitPars_["min"];
  float maxTopMass=fitPars_["max"];

  //read original tree into an ensemble measurement
  for(unsigned int i=0; i<evTree->GetEntriesFast(); i++)
    {
      evTree->GetEntry(i);

      top::EventSummary_t &ev = evHandler.getEvent();
      top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      double topMass=ev.evmeasurements[0];
      if(topMass<minTopMass || topMass>maxTopMass) continue;
      int btagMult = ev.evmeasurements[4]; 
      if(btagMult==0) continue;

      em.evMasses[em.nEvents] = topMass;
      if(btagMult==1) em.evCategories[em.nEvents]=(ev.cat==EMU ? OF_EQ1BTAGS  : SF_EQ1BTAGS);
      else            em.evCategories[em.nEvents]=(ev.cat==EMU ? OF_GEQ2BTAGS : SF_GEQ2BTAGS);
      em.nEvents++;
    }  

  
  //fit mass
  MassFitResults_t res = DoMassFit(em,debug);  
  em.status=true;
  em.mass = res.tMass;
  em.err = res.tMassErr;
  return em;
}

//
MassFitResults_t MassMeasurement::DoMassFit(EnsembleMeasurement_t &em, bool debug)
{
  MassFitResults_t result;
  result.status=false;
      
  //re-define the dataset
  inclusiveData->reset();
  for(int i=0; i<em.nEvents; i++)
    {
      recoMass->setVal( em.evMasses[i] );
      category->setVal( em.evCategories[i] );
      inclusiveData->add(RooArgSet(*recoMass,*category));
    }

  //fit
  result = CombinedMassFitter(debug);
  result.status=true;
      
  //free mem and return the result
  return result;
}

//
void MassMeasurement::InitModel()
{

  int ncategs = (int) fitPars_["ncategs"];
  double tmin(fitPars_["min"]);
  double tmax(fitPars_["max"]);

  //base variables
  recoMass = new RooRealVar( "RecoMass"+tag_, "Reconstructed mMass", tmin,tmax,"GeV/c^{2}");   
  category = new RooRealVar( "Category"+tag_, "Category", 0,MAXFITCATEGORIES);
  topMass  = new RooRealVar( "TopMass"+tag_,"Top quark mass",172,tmin,tmax,"GeV/c^{2}");

  inclusiveData = new RooDataSet("inclusivedata"+tag_,"Inclusive data"+tag_,RooArgSet(*category,*recoMass));

 
  //model prototype (per category)
  for(int icat=0; icat<ncategs; icat++)
    {
      TString sName("s"); sName += icat;

      //calibrated top mass
      char buf[50];
      sprintf(buf,"%f+@0*%f",fitPars_["resbias_"+sName],fitPars_["resslope_"+sName]);
      TString calibformula(buf);
      TString calibname("CalibratedTopMass"+tag_+sName);
      RooFormulaVar *calibtopmass = new RooFormulaVar(calibname,calibformula,*topMass);

      //signal component  
      RooRealVar *g_mean_shift   ( new RooRealVar("#mu_{G}(intercept)_"+tag_+sName,    "Sig Gaus mean (int)"+sName,    fitPars_["#mu_{G}(intercept)_"+sName]) );
      RooRealVar *g_mean_slope   ( new RooRealVar("#mu_{G}(slope)_"+tag_+sName,        "Sig Gaus mean(slope)"+sName,   fitPars_["#mu_{G}(slope)_"+sName]) );
      RooRealVar *g_sigma_shift  ( new RooRealVar("#sigma_{G}(intercept)_"+tag_+sName, "Sig Gaus width (int)"+sName,   fitPars_["#sigma_{G}(intercept)_"+sName]) );
      RooRealVar *g_sigma_slope  ( new RooRealVar("#sigma_{G}(slope)_"+tag_+sName,     "Sig Gaus width (slope)"+sName, fitPars_["#sigma_{G}(slope)_"+sName]) );
      RooFormulaVar *g_mean      ( new RooFormulaVar("g_mean_"+tag_+sName,             "(@0-172)*@1+@2",   RooArgSet(*calibtopmass,*g_mean_slope,*g_mean_shift)));
      RooFormulaVar *g_sigma     ( new RooFormulaVar("g_sigma_"+tag_+sName,            "(@0-172)*@1+@2",   RooArgSet(*calibtopmass,*g_sigma_slope,*g_sigma_shift))); 
      RooGaussian *gaus          ( new RooGaussian("gaus_"+tag_+sName,                 "Mass component #1 " + sName, *recoMass, *g_mean, *g_sigma));
	  
      RooRealVar *l_sigma_shift  ( new RooRealVar("#sigma_{L}(intercept)_"+tag_+sName, "Sig Lan width (int)"+sName,    fitPars_["#sigma_{L}(intercept)_"+sName]));
      RooRealVar *l_sigma_slope  ( new RooRealVar("#sigma_{L}(slope)_"+tag_+sName,     "Sig Lan width (slope)"+sName,  fitPars_["#sigma_{L}(slope)_"+sName]));
      RooRealVar *l_mean_shift   ( new RooRealVar("mpv_{L}(intercept)_"+tag_+sName,    "Sig Lan mpv (int)"+sName,      fitPars_["mpv_{L}(intercept)_"+sName]));
      RooRealVar *l_mean_slope   ( new RooRealVar("mpv_{L}(slope)_"+tag_+sName,        "Sig Lan mpv (slope)"+sName,    fitPars_["mpv_{L}(slope)_"+sName]));
      RooFormulaVar *l_mean      ( new RooFormulaVar("l_mean_"+tag_+sName,             "(@0-172)*@1+@2",   RooArgSet(*calibtopmass,*l_mean_slope,*l_mean_shift)));
      RooFormulaVar *l_sigma     ( new RooFormulaVar( "l_sigma_"+tag_+sName,           "(@0-172)*@1+@2",   RooArgSet(*calibtopmass,*l_sigma_slope,*l_sigma_shift))); 
      RooLandau *lan             ( new RooLandau("lan_"+tag_+sName,                    "Mass component #2" + sName, *recoMass, *l_mean, *l_sigma));  
      RooRealVar *massfrac_shift ( new RooRealVar("#alpha(intercept)_"+tag_+sName,     "Sig Lan frac (int)"+sName,     fitPars_["#alpha(intercept)_"+sName]) );
      RooRealVar *massfrac_slope ( new RooRealVar("#alpha(slope)_"+tag_+sName,         "Sig Lan frac (slope)"+sName,   fitPars_["#alpha(slope)_"+sName]) );
      RooFormulaVar *massfrac    ( new RooFormulaVar("#alpha_"+tag_+sName,             "(@0-172)*@1+@2",   RooArgSet(*calibtopmass,*massfrac_slope,*massfrac_shift)));

      //       RooRealVar *l2_sigma_shift ( new RooRealVar("#sigma_{L2}(intercept)_"+tag_+sName, "Sig Lan 2 width (int)"+sName,    fitPars_["#sigma_{L2}(intercept)_"+sName]));
      //       RooRealVar *l2_sigma_slope ( new RooRealVar("#sigma_{L2}(slope)_"+tag_+sName,     "Sig Lan 2 width (slope)"+sName,  fitPars_["#sigma_{L2}(slope)_"+sName]));
      //       RooRealVar *l2_mean_shift  ( new RooRealVar("mpv_{L2}(intercept)_"+tag_+sName,    "Sig Lan 2 mpv (int)"+sName,      fitPars_["mpv_{L2}(intercept)_"+sName]));
      //       RooRealVar *l2_mean_slope  ( new RooRealVar("mpv_{L2}(slope)_"+tag_+sName,        "Sig Lan 2 mpv (slope)"+sName,    fitPars_["mpv_{L2}(slope)_"+sName]));
      //       RooFormulaVar *l2_mean     ( new RooFormulaVar("l2_mean_"+tag_+sName,             "(@0-172)*@1+@2",   RooArgSet(*calibtopmass,*l2_mean_slope,*l2_mean_shift)));
      //       RooFormulaVar *l2_sigma    ( new RooFormulaVar( "l2_sigma_"+tag_+sName,           "(@0-172)*@1+@2",   RooArgSet(*calibtopmass,*l2_sigma_slope,*l2_sigma_shift))); 
      //       RooLandau *lan2            ( new RooLandau("lan2_"+tag_+sName,                    "Mass component #3" + sName, *recoMass, *l2_mean, *l2_sigma));  
      //       RooRealVar *massfrac2_shift ( new RooRealVar("#beta(intercept)_"+tag_+sName,     "Sig Lan 2 frac (int)"+sName,     fitPars_["#beta(intercept)_"+sName]) );
      //       RooRealVar *massfrac2_slope ( new RooRealVar("#beta(slope)_"+tag_+sName,         "Sig Lan 2 frac (slope)"+sName,   fitPars_["#beta(slope)_"+sName]) );
      //       RooFormulaVar *massfrac2    ( new RooFormulaVar("#beta_"+tag_+sName,             "(@0-172)*@1+@2",   RooArgSet(*calibtopmass,*massfrac2_slope,*massfrac2_shift)));
      
      RooAddPdf *signalmassmodel ( new RooAddPdf("signalmodel_"+tag_+sName,            "Signal Model " + sName, RooArgList(*lan,*gaus),*massfrac));     
      //RooAddPdf *signalmassmodel ( new RooAddPdf("signalmodel_"+tag_+sName,            "Signal Model " + sName, RooArgList(*lan,*lan2,*gaus),RooArgList(*massfrac,*massfrac2)) );     
      signalpdfset.add(*signalmassmodel);

      //background components
      //single top, VV, other ttbar
      RooRealVar *nondybckg_sigma_l   ( new RooRealVar("bckg#sigma_{l}_"+tag_+sName, "Non DY Bckg Lan width" + sName,     fitPars_["nondybckg#sigma_{l}_"+sName]));
      RooRealVar *nondybckg_mpv_l     ( new RooRealVar("bckgmpv_{l}_"+tag_+sName,    "Non DY Bckg Lan mpv",               fitPars_["nondybckgmpv_{l}_"+sName]));
      RooLandau *nondybckg_lan        ( new RooLandau("bckglandau_"+tag_+sName,      "Non DY Bckg Mass component #1" + sName, *recoMass, *nondybckg_mpv_l, *nondybckg_sigma_l));

      //DY->ee/mumu
      RooRealVar *sfdybckg_frac    ( new RooRealVar("sfdybckg#alpha_"+tag_+sName,     "SF DY Bckg fraction"+sName,  fitPars_["sfdybckg#alpha_"+sName]));
      RooRealVar *sfdybckg_sigma_l ( new RooRealVar("sfdybckg#sigma_{l}_"+tag_+sName, "SF DY Bckg Lan width"+sName, fitPars_["sfdybckg#sigma_{l}_"+sName]));
      RooRealVar *sfdybckg_mpv_l   ( new RooRealVar("sfdybckgmpv_{l}_"+tag_+sName,    "SF DY Bckg Lan mpv"+sName,   fitPars_["sfdybckgmpv_{l}_"+sName]));
      RooLandau *sfdybckg_lan      ( new RooLandau("sfdybckgmodel_"+tag_+sName,       "SF DY Bckg Model"+sName, *recoMass, *sfdybckg_mpv_l, *sfdybckg_sigma_l));

      //DY->tautau
      RooRealVar *ofdybckg_frac    ( new RooRealVar("ofdybckg#alpha_"+tag_+sName,     "OF DY Bckg fraction"+sName,  fitPars_["ofdybckg#alpha_"+sName]));
      RooRealVar *ofdybckg_sigma_l ( new RooRealVar("ofdybckg#sigma_{l}_"+tag_+sName, "OF DY Bckg Lan width"+sName, fitPars_["ofdybckg#sigma_{l}_"+sName]));
      RooRealVar *ofdybckg_mpv_l   ( new RooRealVar("ofdybckgmpv_{l}_"+tag_+sName,    "OF DY Bckg Lan mpv"+sName,   fitPars_["ofdybckgmpv_{l}_"+sName]));
      RooLandau  *ofdybckg_lan     ( new RooLandau("ofdybckgmodel_"+tag_+sName,       "OF DY Bckg Model"+sName, *recoMass, *ofdybckg_mpv_l, *ofdybckg_sigma_l));

      //combined background
      RooAddPdf *bckgmassmodel     ( new RooAddPdf("bckgmodel_"+tag_+sName,"Background Model",
						   RooArgList(*sfdybckg_lan,*ofdybckg_lan,*nondybckg_lan),
						   RooArgList(*sfdybckg_frac,*ofdybckg_frac)));
      bckgpdfset.add(*bckgmassmodel);

      //the data model
      RooRealVar *nsigvar        ( new RooRealVar("SignalYields"+tag_+sName, "Signal yield"+sName,1,0,10) );
      sigYieldParams.add(*nsigvar);
      RooRealVar *nbkgvar        ( new RooRealVar("BackgroundYields"+tag_+sName,"Background yield"+sName,0,0,10) );
      constrParams.add(*nbkgvar);
      RooAddPdf *shapeModel     ( new RooAddPdf("shapemodel_"+tag_+sName,"signal+background",RooArgList(*bckgmassmodel,*signalmassmodel),RooArgList(*nbkgvar,*nsigvar)));
      RooRealVar *bckg_mean_constraint  ( new RooRealVar("bckgcnstr#mu_{g}_"+tag_+sName,"Mean of bckg cnstr",fitPars_["#mu_{background}_"+sName]));
      RooRealVar *bckg_sigma_constraint ( new RooRealVar("bckgcnstr#sigma_{g}_"+tag_+sName,"Sigma of bckg cnstr",fitPars_["#sigma_{background}_"+sName]));
      RooGaussian *bckgEventConstraint  ( new RooGaussian("bckgconstraintpdf_"+tag_+sName,"Bckg constraint",*nbkgvar,*bckg_mean_constraint,*bckg_sigma_constraint)); 

      //signal+background fit
      RooAbsPdf *model=0;
      if(fitPars_["#mu_{background}_"+sName]>0)
	model = new RooProdPdf("model_"+tag_+sName,"(signal+background)*evconstraint*bkgconstraint",RooArgSet(*bckgEventConstraint,*shapeModel));
      //signal only fit
      else
	model = new RooAddPdf("model_"+sName,"Signal Only Model",RooArgList(*lan,*gaus),*massfrac);

      //      model->printCompactTree(cout);
      allPdfs.push_back( model );  
    }
}


//
MassFitResults_t MassMeasurement::CombinedMassFitter(bool debug)
{
  MassFitResults_t result;

  //
  //configure the fit parameters
  //
  int ncategs = (int) fitPars_["ncategs"];
  double tmin(fitPars_["min"]);
  double tmax(fitPars_["max"]);
 
  //  cout << "[MassMeasurement::CombinedMassFitter] with " << ncategs << " categories and calibration" << flush;
  //   calibtopmass->dumpFormula();
  //   cout << endl;

  //
  //import data and define likelihoods per category
  //
  int totalEventsUsed(0);

  //reset previous 
  RooArgSet llSet;  
  for(size_t ill=0; ill<allLL.size(); ill++)    delete allLL[ill];       allLL.clear();
  for(size_t id=0; id<allData.size(); id++)     delete allData[id];      allData.clear();
  for(size_t ic=0; ic<allCaptions.size(); ic++) delete allCaptions[ic];  allCaptions.clear();

  for(int icat=0; icat<ncategs; icat++)
    {
      TString sName("s"); sName += icat;
     
      TString cutFormula("Category"+tag_+"=="); 
      cutFormula += icat;
      RooDataSet *sData = (RooDataSet *) inclusiveData->reduce(cutFormula);
      if(getCategorizationMode()==INCLUSIVE)
	{
	  cutFormula = "Category"+tag_+"==";
	  cutFormula += (icat==0 ? OF_EQ1BTAGS : OF_GEQ2BTAGS);
	  sData->append( *((RooDataSet *) inclusiveData->reduce(cutFormula)) );
	}
      allData.push_back(sData);
      totalEventsUsed += sData->sumEntries();
      //sData->Print("v");
      
      //update values
      topMass->setVal(172);
      RooRealVar *nsigvar= (RooRealVar *)sigYieldParams.find("SignalYields"+tag_+sName);
      nsigvar->setVal(0.85*sData->sumEntries());
      nsigvar->setRange(0,sData->sumEntries());
      RooRealVar *nbkgvar= (RooRealVar *)constrParams.find("BackgroundYields"+tag_+sName);
      nbkgvar->setVal(fitPars_["#mu_{background}_"+sName]);
      nbkgvar->setRange(0,sData->sumEntries());

      //signal+background fit
      RooAbsPdf *model=allPdfs[icat];
      if(fitPars_["#mu_{background}_"+sName]>0)
	{
	  RooNLLVar *nll = (RooNLLVar*) model->createNLL(*sData,RooFit::CloneData(kFALSE),Extended(kTRUE),Constrain(*nbkgvar));  
	  //RooNLLVar *nll = (RooNLLVar*) model->createNLL(*sData,RooFit::CloneData(kFALSE),Extended(kTRUE));  
	  RooMinuit minuit(*nll);
	  minuit.setVerbose(false);
	  minuit.migrad();
	  minuit.hesse();
	  allLL.push_back(nll);
	  model->fitTo(*sData,Extended(kTRUE),Constrain(*nbkgvar),Minos(),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(tmin,tmax));	  
	  //model->fitTo(*sData,Extended(kTRUE),Minos(),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(tmin,tmax));	  
	  result.iTmass[icat]     = topMass->getVal();
	  result.iTmassErr[icat] = topMass->getError();
	  
	  //prepare label
	  TPaveText *pave = new TPaveText(0.65,0.75,0.9,0.92,"NDC");
	  pave->SetTextFont(42);
	  pave->SetFillStyle(0);
	  pave->SetBorderSize(0);
	  TString catKey("cat"); catKey += icat;
	  pave->AddText( catTitles_[icat] )->SetTextAlign(11);

	  char buf[100];
	  sprintf(buf,"m_{top}=%3.1f^{+%3.1f}_{%3.1f}",topMass->getVal(),topMass->getAsymErrorHi(),topMass->getAsymErrorLo());
	  pave->AddText(buf)->SetTextAlign(11);
	  sprintf(buf,"N_{signal}=%3.1f^{+%3.1f}_{%3.1f}",nsigvar->getVal(),nsigvar->getAsymErrorHi(),nsigvar->getAsymErrorLo());
	  pave->AddText(buf)->SetTextAlign(11);
	  sprintf(buf,"N_{background}=%3.1f^{+%3.1f}_{%3.1f}",nbkgvar->getVal(),nbkgvar->getAsymErrorHi(),nbkgvar->getAsymErrorLo());
	  pave->AddText(buf)->SetTextAlign(11);
	  allCaptions.push_back(pave);
	}
      
      //signal only fit
      else 
	{
	  RooNLLVar *nll = (RooNLLVar*) model->createNLL(*sData,RooFit::CloneData(kFALSE));
	  RooMinuit minuit(*nll);
	  minuit.setVerbose(false);
	  minuit.migrad();
	  minuit.hesse();
	  allLL.push_back(nll);
	  model->fitTo(*sData,Minos(),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(tmin,tmax));
	  result.iTmass[icat]     = topMass->getVal();
	  result.iTmassErr[icat] = topMass->getError();

	  //prepare label     
	  TPaveText *pave = new TPaveText(0.7,0.75,0.95,0.93,"NDC");
	  pave->SetFillStyle(0);
	  pave->SetBorderSize(0);
	  TString catKey("cat"); catKey += icat;
	  TString catTitle=catTitles_[icat];
	  pave->AddText( catTitle.Data() )->SetTextAlign(11);
	  char buf[100];
	  sprintf(buf,"m_{top}=%3.1f^{+%3.1f}_{%3.1f}",topMass->getVal(),topMass->getAsymErrorHi(),topMass->getAsymErrorLo());
	  pave->AddText(buf)->SetTextAlign(11);
	  allCaptions.push_back(pave);
	}    
	  
      llSet.add((RooNLLVar &)(*allLL[icat]));
	  
      // signal and background are set constant for the combined fit
      //nsigvar->setConstant();
      //nbkgvar->setConstant();
    }

       
  //minimize the combined likelihood
  RooAbsReal *combll;
  if(ncategs==1) combll = allLL[0];
  else 
    {
      combll = new RooAddition("combll"+tag_,"combll"+tag_,llSet);	  
      RooMinuit minuit(*combll); 
      minuit.setVerbose(false);
      minuit.setErrorLevel(0.5); //otherwise RooFit seems to assume the chi^2 default
      minuit.hesse();
      minuit.migrad();
      minuit.setErrorLevel(0.5);
      minuit.hesse();
      minuit.minos();
      minuit.save();
    }

  //save the result
  //cout << " ***** Total events used: " << totalEventsUsed << endl;
  result.tMass        = topMass->getVal();  
  result.tMassErr     = topMass->getError(); 
  result.tMassErrHigh = topMass->getAsymErrorHi(); 
  result.tMassErrLo   = topMass->getAsymErrorLo();

  //
  //plot the result
  //
  if(debug)
    {
      setStyle();

      //mass distributions
      TCanvas *c = new TCanvas("massfitter","Fit Result",1200,600*ncategs/2);
      c->SetWindowSize(1200,600*ncategs/2);
      c->Divide(2,ncategs/2);
      for(int icat=0; icat<ncategs; icat++)
	{
	  TString sName("s"); sName += icat;
	  c->cd(icat+1);

	  float uncalibTopMass=(result.iTmass[icat]-fitPars_["resbias_"+sName])/fitPars_["resslope_"+sName];
	  //topMass->setVal(result.iTmass[icat]);
	  topMass->setVal(uncalibTopMass);
	  RooPlot* frame = recoMass->frame(Title(sName));
	  allData[icat]->plotOn(frame,Binning(20),DrawOption("pz"));
	  allPdfs[icat]->plotOn(frame,Components("bckgmodel*"),DrawOption("lf"),FillStyle(1001),FillColor(kGray),LineColor(kGray),Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
	  //allPdfs[icat]->plotOn(frame,Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
	  allPdfs[icat]->plotOn(frame,NumEvents(allData[icat]->sumEntries()),MoveToBack());
	  frame->GetXaxis()->SetTitle("Reconstructed Mass [GeV/c^{2}]");
	  frame->GetXaxis()->SetTitleOffset(1.2);
	  frame->GetYaxis()->SetTitle("Events");
	  frame->GetYaxis()->SetTitleOffset(1.2);
	  frame->Draw();
	  allCaptions[icat]->Draw();
	}
      c->SaveAs("MassFitResults.C");
      c->SaveAs("MassFitResults.pdf");
      c->SaveAs("MassFitResults.png");
      delete c;

      //likelihoods
      
      topMass->setVal(result.tMass);
      float minForLL=topMass->getVal()-5*topMass->getError();
      float maxForLL=topMass->getVal()+5*topMass->getError();
      c = new TCanvas("massfitterll","Fit likelihood",600,600);
      c->SetWindowSize(600,600);
      RooPlot *frame = topMass->frame(Bins(100),Range(minForLL,maxForLL),Title("Likelihood")) ;    
      Int_t catColors[]={831,809,590,824};
      for(int icat=0; icat<ncategs; icat++) 
	{
	  TString catTag("cat"); catTag += icat; catTag+=tag_;
	  allLL[icat]->plotOn(frame,Name(catTag),ShiftToZero(),FillStyle(0), LineColor(catColors[icat]), LineStyle(icat%2==0 ? 9 : 1), LineWidth(2),Name(catTag));
	}
      combll->plotOn(frame,ShiftToZero(),FillStyle(0),LineWidth(3),Name("comb"+tag_));
      frame->GetXaxis()->SetTitle("Top Mass [GeV/c^{2}]");
      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetYaxis()->SetTitle("-log(L/L_{max})");
      frame->GetYaxis()->SetTitleOffset(1);
      frame->GetYaxis()->SetRangeUser(0,5);
      frame->Draw();

      TLegend *leg = new TLegend(0.6,0.65,0.9,0.9,NULL,"brNDC");
      formatForCmsPublic(c,leg,"",3);
      leg->AddEntry("comb"+tag_,"Combined","l");
      for(int icat=0; icat<ncategs; icat++) 
	{
	  TString catTag("cat"); catTag += icat; catTag += tag_;
	  leg->AddEntry(catTag,catTitles_[icat],"l");
	}
      leg->SetFillColor(0);
      leg->SetFillStyle(3001);
      leg->Draw();

      //prepare label     
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
      char buf[100];
      sprintf(buf,"m_{top}=%3.2f^{+%3.2f}_{%3.2f}",topMass->getVal(),topMass->getAsymErrorHi(),topMass->getAsymErrorLo());
      pave->AddText(buf);
      pave->Draw();

      c->Modified();
      c->Update();
      c->SaveAs("MassFitLikelihood.C");
      c->SaveAs("MassFitLikelihood.png");
      c->SaveAs("MassFitLikelihood.pdf");
      delete c;

      //the summary : inclusive distribution + inclusive likelihood
      c = new TCanvas("incmassfitter","Inclusive Fit Result",600,600);
      c->SetWindowSize(600,600);
      frame = recoMass->frame(Title("inclusivesample"));
      inclusiveData->plotOn(frame,Binning(20),DrawOption("pz"));

      RooAddPdf incbckgmassmodel("incbckg"+tag_,"Inclusive Background Model",bckgpdfset,constrParams);
      incbckgmassmodel.plotOn(frame,DrawOption("lf"),FillStyle(1001),FillColor(kGray),LineColor(kGray),Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
      signalpdfset.add(bckgpdfset);
      sigYieldParams.add(constrParams);
      RooAddPdf incmassmodel("incmodel"+tag_,"Inclusive Model",signalpdfset,sigYieldParams);
      incmassmodel.plotOn(frame,NumEvents(inclusiveData->sumEntries()),MoveToBack());
      frame->GetXaxis()->SetTitle("Reconstructed Mass [GeV/c^{2}]");
      frame->GetXaxis()->SetTitleOffset(1.2);
      frame->GetYaxis()->SetTitle("Events");
      frame->GetYaxis()->SetTitleOffset(1.2);
      frame->Draw();

      leg = c->BuildLegend();
      sprintf(buf,"m_{top}=%3.1f #pm %3.1f GeV/c^{2}",topMass->getVal(),topMass->getError());
      formatForCmsPublic(c,leg,buf,1);
      leg->Delete();

      pave = new TPaveText(0.15,0.96,0.51,0.99,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      pave->SetTextAlign(12);
      pave->SetTextFont(42);
      pave->AddText("CMS preliminary");
      pave->Draw();
	  
      TPad *npad = new TPad("llpad","ll", 0.56, 0.60, 0.9, 0.94);
      npad->Draw();
      npad->cd();
      frame = topMass->frame(Bins(100),Range(minForLL,maxForLL),Title("Likelihood")) ;    
      combll->plotOn(frame,ShiftToZero(),Name("inccomb"+tag_));
      frame->GetXaxis()->SetTitle("Top Quark Mass [GeV/c^{2}]");
      frame->GetXaxis()->SetTitleOffset(1.2);
      frame->GetYaxis()->SetTitle("-log(L/L_{max})");
      frame->GetYaxis()->SetTitleOffset(1);
      frame->GetYaxis()->SetRangeUser(0,5);
      frame->Draw();

      c->SaveAs("MassFitSummary.C");
      c->SaveAs("MassFitSummary.png");
      c->SaveAs("MassFitSummary.pdf");

      delete c;
    }
      
  delete combll;

  //all done!
  return result;
}
    
   
//
std::map<TString,Double_t> MassMeasurement::ParseParametersFrom(TString parfileURL)
{
  fitPars_.clear();
  ifstream in;
  in.open(parfileURL);
  TString line;
  while (1) {
    in >> line;
    if (!in.good()) break;
    TObjArray *tokens = line.Tokenize(":");
    if(tokens->GetEntriesFast()<2) continue;
    TString key = ((TObjString *)tokens->At(0))->GetString();
    TString val = ((TObjString *)tokens->At(1))->GetString();
    fitPars_[key]=val.Atof();
  }
  in.close();

  //translate the categories
  int cattype=(int) fitPars_["cattype"];
  if(cattype==INCLUSIVE)
    {
      catTitles_[EQ1BTAGS]  = TString("1 b-tags");
      catTitles_[GEQ2BTAGS] = TString("#geq 2 b-tags");
    }
  else
    {
      catTitles_[SF_EQ1BTAGS]  = TString("1 b-tags (ee+#mu#mu)");
      catTitles_[SF_GEQ2BTAGS] = TString("#geq 2 b-tags (ee+#mu#mu)");
      catTitles_[OF_EQ1BTAGS]  = TString("1 b-tags (e#mu)");
      catTitles_[OF_GEQ2BTAGS] = TString("#geq 2 b-tags (e#mu)");
    }

  return fitPars_;
}



