#include  <string.h>
#include "LIP/Top/interface/MassMeasurement.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "RooAddition.h"
#include "RooNLLVar.h"
#include "TPaveText.h"
#include <memory>

using namespace std;
using namespace RooFit;
using namespace top;

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
  int cattype=(int) fitPars_["cattype"];
  float minTopMass=fitPars_["min"];
  float maxTopMass=fitPars_["max"];

  //read original tree into an ensemble measurement
  for(unsigned int i=0; i<evTree->GetEntriesFast(); i++)
    {
      evTree->GetEntry(i);

      top::EventSummary_t &ev = evHandler.getEvent();
      top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      //      if(phys.met.pt()<20) continue;
      double topMass=ev.evmeasurements[0];
      if(topMass<minTopMass || topMass>maxTopMass) continue;
      int btagMult = ev.evmeasurements[4]; 
      
      em.evMasses[em.nEvents] = topMass;
      if(cattype==0) em.evSample[em.nEvents]=0;
      if(cattype==1) em.evSample[em.nEvents] = ( ev.cat==EMU );
      if(cattype==2) em.evSample[em.nEvents] = ( btagMult>=2 ? 2 : btagMult );
      em.nEvents++;
    }  

  
  //fit mass
  MassFitResults_t res = DoMassFit(em,debug);  
  em.status=true;
  em.mass = res.tMass;
  em.err = res.tMassErr;
  em.errLo = res.tMassErrLo;
  em.errHigh = res.tMassErrHigh;
  
  return em;
}

//
MassFitResults_t MassMeasurement::DoMassFit(EnsembleMeasurement_t &em, bool debug)
{
  MassFitResults_t result;
  result.status=false;
      
  //define a simple tree with the masses and the sample
  TTree *masspts = new TTree("masspts","masspts");
  Float_t imass;
  masspts->Branch("Mass",&imass,"Mass/F");
  Int_t isample;
  masspts->Branch("Sample",&isample,"Sample/I");
  for(int i=0; i<em.nEvents; i++)
    {
      imass=em.evMasses[i];
      isample=em.evSample[i];
      masspts->Fill();
    }

  //fit
  result = CombinedMassFitter(masspts,debug);
  result.status=true;
      
  //free mem
  masspts->Delete();
  return result;
}


//
MassFitResults_t MassMeasurement::CombinedMassFitter(TTree *data, bool debug)
{
  MassFitResults_t result;

  //the fit parameters
  int ncategs = (int) fitPars_["ncategs"];
  double tmin(fitPars_["min"]),tmax(fitPars_["max"]);
  RooRealVar Mass( "Mass", "Mass", tmin,tmax,"GeV/c^{2}");   
  RooRealVar Sample( "Sample", "Sample", 0,fitPars_["ncategs"]);
  RooRealVar topmass( "m_{t}","Measured Mass",172,tmin,tmax,"GeV/c^{2}");
  char buf[50];
  sprintf(buf,"%f+@%f",fitPars_["resbias"],fitPars_["resslope"]);
  TString calibformula(buf);
  RooFormulaVar calibtopmass( "m_{t}^{calib}",calibformula,topmass);
  RooArgSet constrParams,sigYieldParams;
  cout << "[MassMeasurement::CombinedMassFitter] with " << ncategs << " categories and calibration=" << buf << endl;
  int totalEventsUsed(0);

  //
  //IMPORT DATA and DEFINE THE MODEL PROTOTYPE (per category)
  //
  std::vector<RooDataSet * > allData;
  std::vector<RooAbsPdf * > allPdfs;
  RooArgSet bckgpdfset,signalpdfset;
  std::vector<TPaveText * > allCaptions;
  std::vector<TString> allCategTitles;
  std::vector<RooNLLVar *> allLL;
  RooDataSet inclusiveData("inclusivedata","Inclusive data",RooArgSet(Sample,Mass),Import(*data));
  RooArgSet llSet;
  for(int icat=0; icat<ncategs; icat++)
    {
      TString sName("s"); sName += icat;
      TString cutFormula("Sample=="); cutFormula += icat;
      RooDataSet *sData = new RooDataSet("data_"+sName,"Data in " + cutFormula,RooArgSet(Sample,Mass),Import(*data),Cut(cutFormula));
      allData.push_back(sData);
      totalEventsUsed += sData->sumEntries();
      sData->Print("v");

      //signal component  
      RooRealVar *massfrac_shift ( new RooRealVar("#alpha(intercept)_"+sName,     "Sig Lan frac (int)"+sName,     fitPars_["#alpha(intercept)_"+sName]) );
      RooRealVar *massfrac_slope ( new RooRealVar("#alpha(slope)_"+sName,         "Sig Lan frac (slope)"+sName,   fitPars_["#alpha(slope)_"+sName]) );
      RooRealVar *g_mean_shift   ( new RooRealVar("#mu_{G}(intercept)_"+sName,    "Sig Gaus mean (int)"+sName,    fitPars_["#mu_{G}(intercept)_"+sName]) );
      RooRealVar *g_mean_slope   ( new RooRealVar("#mu_{G}(slope)_"+sName,        "Sig Gaus mean(slope)"+sName,   fitPars_["#mu_{G}(slope)_"+sName]) );
      RooRealVar *g_sigma_shift  ( new RooRealVar("#sigma_{G}(intercept)_"+sName, "Sig Gaus width (int)"+sName,   fitPars_["#sigma_{G}(intercept)_"+sName]) );
      RooRealVar *g_sigma_slope  ( new RooRealVar("#sigma_{G}(slope)_"+sName,     "Sig Gaus width (slope)"+sName, fitPars_["#sigma_{G}(slope)_"+sName]) );
      RooFormulaVar *massfrac    ( new RooFormulaVar("#alpha_"+sName,             "(@0-172)*@1+@2",   RooArgSet(calibtopmass,*massfrac_slope,*massfrac_shift)));
      RooFormulaVar *g_mean      ( new RooFormulaVar("g_mean_"+sName,             "(@0-172)*@1+@2",   RooArgSet(calibtopmass,*g_mean_slope,*g_mean_shift)));
      RooFormulaVar *g_sigma     ( new RooFormulaVar("g_sigma_"+sName,            "(@0-172)*@1+@2",   RooArgSet(calibtopmass,*g_sigma_slope,*g_sigma_shift))); 
      RooGaussian *gaus          ( new RooGaussian("gaus_"+sName,                 "Mass component #1 " + sName, Mass, *g_mean, *g_sigma));
	  
      RooRealVar *l_sigma_shift  ( new RooRealVar("#sigma_{L}(intercept)_"+sName, "Sig Lan width (int)"+sName,    fitPars_["#sigma_{L}(intercept)_"+sName]));
      RooRealVar *l_sigma_slope  ( new RooRealVar("#sigma_{L}(slope)_"+sName,     "Sig Lan width (slope)"+sName,  fitPars_["#sigma_{L}(slope)_"+sName]));
      RooRealVar *l_mean_shift   ( new RooRealVar("mpv_{L}(intercept)_"+sName,    "Sig Lan mpv (int)"+sName,      fitPars_["mpv_{L}(intercept)_"+sName]));
      RooRealVar *l_mean_slope   ( new RooRealVar("mpv_{L}(slope)_"+sName,        "Sig Lan mpv (slope)"+sName,    fitPars_["mpv_{L}(slope)_"+sName]));
      RooFormulaVar *l_mean      ( new RooFormulaVar("l_mean_"+sName,             "(@0-172)*@1+@2",   RooArgSet(calibtopmass,*l_mean_slope,*l_mean_shift)));
      RooFormulaVar *l_sigma     ( new RooFormulaVar( "l_sigma_"+sName,           "(@0-172)*@1+@2",   RooArgSet(calibtopmass,*l_sigma_slope,*l_sigma_shift))); 
      RooLandau *lan             ( new RooLandau("lan_"+sName,                    "Mass component #2" + sName, Mass, *l_mean, *l_sigma));  

      RooAddPdf *signalmassmodel ( new RooAddPdf("signalmodel_"+sName,            "Signal Model " + sName, RooArgList(*lan,*gaus),*massfrac));     
      signalpdfset.add(*signalmassmodel);

      //background component
      RooRealVar *nondybckg_massfrac  ( new RooRealVar("bckg#alpha_"+sName,     "Non DY Bckg Lan frac" + sName,      fitPars_["bckg#alpha_"+sName]));
      RooRealVar *nondybckg_mean_g    ( new RooRealVar("bckg#mu_{g}_"+sName,    "Non DY Bckg Gaus mean" + sName,     fitPars_["bckg#mu_{g}_"+sName]));
      RooRealVar *nondybckg_sigma_g   ( new RooRealVar("bckg#sigma_{g}_"+sName, "Non DY Bckg Gaus width" + sName,    fitPars_["bckg#sigma_{g}_"+sName]));
      RooRealVar *nondybckg_sigma_l   ( new RooRealVar("bckg#sigma_{l}_"+sName, "Non DY Bckg Lan width" + sName,     fitPars_["bckg#sigma_{l}_"+sName]));
      RooRealVar *nondybckg_mpv_l     ( new RooRealVar("bckgmpv_{l}_"+sName,    "Non DY Bckg Lan mpv",               fitPars_["bckgmpv_{l}_"+sName]));
      RooLandau *nondybckg_lan        ( new RooLandau("bckglandau_"+sName,      "Non DY Bckg Mass component #1" + sName, Mass, *nondybckg_mpv_l, *nondybckg_sigma_l));
      RooGaussian *nondybckg_gaus     ( new RooGaussian("bckggauss_"+sName,     "Non DY Bckg Mass component #2" + sName, Mass, *nondybckg_mean_g, *nondybckg_sigma_g));   
      RooAddPdf *nondybckg_massmodel  ( new RooAddPdf("otherbckgmodel_"+sName,  "Non DY Bckg Model", RooArgList(*nondybckg_lan,*nondybckg_gaus), *nondybckg_massfrac));

      RooRealVar *dybckg_frac    ( new RooRealVar("dybckg#alpha_"+sName,     "DY Bckg fraction"+sName,  fitPars_["dybckg#alpha_"+sName]));
      RooRealVar *dybckg_sigma_l ( new RooRealVar("dybckg#sigma_{l}_"+sName, "DY Bckg Lan width"+sName, fitPars_["dybckg#sigma_{l}_"+sName]));
      RooRealVar *dybckg_mpv_l   ( new RooRealVar("dybckgmpv_{l}_"+sName,    "DY Bckg Lan mpv"+sName,   fitPars_["dybckgmpv_{l}_"+sName]));
      RooLandau *dybckg_massmodel ( new RooLandau("dybckgmodel_"+sName,      "DY Bckg Model"+sName, Mass, *dybckg_mpv_l, *dybckg_sigma_l));

      RooAddPdf *bckgmassmodel      ( new RooAddPdf("bckgmodel_"+sName,"Background Model",RooArgList(*dybckg_massmodel,*nondybckg_massmodel),*dybckg_frac));
      bckgpdfset.add(*bckgmassmodel);

      //the data model
      RooRealVar *nsigvar        ( new RooRealVar("N_{signal}_"+sName, "Signal yield"+sName,0.85*sData->sumEntries(),0,sData->sumEntries()));
      RooRealVar *nbkgvar        ( new RooRealVar("N_{background}_"+sName,"Background yield"+sName,fitPars_["#mu_{background}_"+sName],0,sData->sumEntries()));
      RooAddPdf *shapeModel     ( new RooAddPdf("shapemodel_"+sName,"signal+background",RooArgList(*bckgmassmodel,*signalmassmodel),RooArgList(*nbkgvar,*nsigvar)));
      RooRealVar *bckg_mean_constraint  ( new RooRealVar("bckgcnstr#mu_{g}_"+sName,"Mean of bckg cnstr",fitPars_["#mu_{background}_"+sName]));
      RooRealVar *bckg_sigma_constraint ( new RooRealVar("bckgcnstr#sigma_{g}_"+sName,"Sigma of bckg cnstr",fitPars_["#sigma_{background}_"+sName]));
      RooGaussian *bckgEventConstraint  ( new RooGaussian("bckgconstraintpdf_"+sName,"Bckg constraint",*nbkgvar,*bckg_mean_constraint,*bckg_sigma_constraint)); 
      constrParams.add(*nbkgvar);
      sigYieldParams.add(*nsigvar);
      if(fitPars_["#mu_{background}_"+sName]>0)
	{
	  RooAbsPdf *model = new RooProdPdf("model_"+sName,"(signal+background)*evconstraint*bkgconstraint",RooArgSet(*bckgEventConstraint,*shapeModel));	      
	  //RooAbsPdf *model = shapeModel;
	  allPdfs.push_back(model);

	  RooNLLVar *nll = (RooNLLVar*) model->createNLL(*sData,RooFit::CloneData(kFALSE),Extended(kTRUE),Constrain(*nbkgvar));  
	  //RooNLLVar *nll = (RooNLLVar*) model->createNLL(*sData,RooFit::CloneData(kFALSE),Extended(kTRUE));  
	  RooMinuit minuit(*nll); 
	  minuit.migrad();
	  minuit.hesse();
	  allLL.push_back(nll);
	  model->fitTo(*sData,Extended(kTRUE),Constrain(*nbkgvar),Minos(),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(tmin,tmax));	  
	  //model->fitTo(*sData,Extended(kTRUE),Minos(),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(tmin,tmax));	  

	  //prepare label
	  TPaveText *pave = new TPaveText(0.7,0.75,0.95,0.93,"NDC");
	  pave->SetFillStyle(0);
	  pave->SetBorderSize(0);
	  TString catKey("cat"); catKey += icat;
	  TString catTitle("");
	  for(std::map<TString,Double_t>::iterator fpIt = fitPars_.begin();
	      fpIt != fitPars_.end();
	      fpIt++)
	    {
	      if( !fpIt->first.Contains( catKey ) || !fpIt->first.Contains(")") ) continue;
	      TObjArray *tokens = fpIt->first.Tokenize(")");
	      catTitle=((TObjString *)tokens->At(1))->GetString();
	      pave->AddText( catTitle.Data() )->SetTextAlign(11);
	    }
	  char buf[100];
	  sprintf(buf,"m_{top}=%3.1f^{+%3.1f}_{%3.1f}",topmass.getVal(),topmass.getAsymErrorHi(),topmass.getAsymErrorLo());
	  pave->AddText(buf)->SetTextAlign(11);
	  sprintf(buf,"N_{signal}=%3.1f^{+%3.1f}_{%3.1f}",nsigvar->getVal(),nsigvar->getAsymErrorHi(),nsigvar->getAsymErrorLo());
	  pave->AddText(buf)->SetTextAlign(11);
	  sprintf(buf,"N_{background}=%3.1f^{+%3.1f}_{%3.1f}",nbkgvar->getVal(),nbkgvar->getAsymErrorHi(),nbkgvar->getAsymErrorLo());
	  pave->AddText(buf)->SetTextAlign(11);
	  allCaptions.push_back(pave);
	  allCategTitles.push_back(catTitle);
	}
      else 
	{
	  RooAbsPdf *model = new RooAddPdf("model_"+sName,"Signal Only Model",RooArgList(*lan,*gaus),*massfrac);
	  allPdfs.push_back( model );
	  RooNLLVar *nll = (RooNLLVar*) model->createNLL(*sData,RooFit::CloneData(kFALSE));
	  RooMinuit minuit(*nll);
	  minuit.migrad();
	  minuit.hesse();
	  allLL.push_back(nll);
	  model->fitTo(*sData,Minos(),Save(kTRUE),PrintLevel(-1),Verbose(false),Range(tmin,tmax));

	  //prepare label     
	  TPaveText *pave = new TPaveText(0.7,0.75,0.95,0.93,"NDC");
	  pave->SetFillStyle(0);
	  pave->SetBorderSize(0);
	  TString catKey("cat"); catKey += icat;
	  TString catTitle("");
	  for(std::map<TString,Double_t>::iterator fpIt = fitPars_.begin();
	      fpIt != fitPars_.end();
	      fpIt++)
	    {
	      if( !fpIt->first.Contains( catKey ) || !fpIt->first.Contains(")") ) continue;
	      TObjArray *tokens = fpIt->first.Tokenize(")");
	      catTitle=((TObjString *)tokens->At(1))->GetString();
	      pave->AddText( catTitle.Data() )->SetTextAlign(11);
	    }
	  char buf[100];
	  sprintf(buf,"m_{top}=%3.1f^{+%3.1f}_{%3.1f}",topmass.getVal(),topmass.getAsymErrorHi(),topmass.getAsymErrorLo());
	  pave->AddText(buf)->SetTextAlign(11);
	  allCaptions.push_back(pave);
	  allCategTitles.push_back(catTitle);	      
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
      combll = new RooAddition("combll","combll",llSet);	  
      RooMinuit minuit(*combll); 
      minuit.setErrorLevel(0.5); //otherwise it seems to assume the chi^2 default
      minuit.hesse();
      minuit.migrad();
      minuit.setErrorLevel(0.5);
      minuit.hesse();
    }

  //save the result
  cout << " ***** Total events used: " << totalEventsUsed << endl;
  result.tMass=topmass.getVal();                    // top mass measured from data
  result.tMassErr=topmass.getError(); 
  result.tMassErrHigh = topmass.getAsymErrorHi(); 
  result.tMassErrLo = topmass.getAsymErrorLo();
  //cout << topmass.getVal() << " " << topmass.getError() << " " << topmass.getAsymErrorHi() << " " << topmass.getAsymErrorLo() << " " << result.tMassErr << endl;
  //result.nSig=nsigvar.getVal();                      //number of signal events
  //      result.nSigErr=nsigvar.getError();  
  //      result.nSigErrHigh = nsigvar.getAsymErrorHi();  
  //      result.nSigErrLo = nsigvar.getAsymErrorLo();
  //      result.nBckg=nbkgvar.getVal();                     //number of background events
  //      result.nBckgErr=nbkgvar.getError(); 
  //      result.nBckgErrHigh = nbkgvar.getAsymErrorHi(); 
  //      result.nBckgErrLo = nbkgvar.getAsymErrorLo();
  //      RooPlot *frame2 = topmass.frame(Bins(50),Range(tmin,tmax),Title("Likelihood")) ;    
  //      RooCurve *nllplot= nll->plotOn(frame2,ShiftToZero())->getCurve();
  //      result.ll = nllplot->Eval(result.tMass);           //the max. likelihood value
  //      result.lltmass = nllplot->Eval(172);              // likelihood at the expected top mass
  //      result.lltmassHigh = nllplot->Eval(172+1.63);    
  //      result.lltmassLo = nllplot->Eval(172-1.63);
  //      result.npoints=combData.sumEntries();

  if(debug)
    {
      setStyle();
      TCanvas *c = getNewCanvas("massfitter","Fit Result",false);
      c->cd();
      c->Clear();
      c->SetWindowSize(500*(ncategs+1),500);
      c->Divide(ncategs+1,1);
	  
      //mass distributions
      for(int icat=0; icat<ncategs; icat++)
	{
	  TString sName("s"); sName += icat;
	  TPad *p = (TPad *)c->cd(icat+1);
	  p->SetGridx();
	  p->SetGridy();
	  RooPlot* frame = Mass.frame(Title(sName));
	  allData[icat]->plotOn(frame,Binning(10),DrawOption("pz"));
	  allPdfs[icat]->plotOn(frame,Components("bckgmodel*"),DrawOption("f"),FillStyle(3001),Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
	  allPdfs[icat]->plotOn(frame,Components("bckgmodel*,signalmodel*"),LineStyle(kDotted),Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
	  frame->GetXaxis()->SetTitle("Reconstructed Mass [GeV/c^{2}]");
	  frame->GetXaxis()->SetTitleOffset(0.8);
	  frame->GetYaxis()->SetTitle("Events");
	  frame->GetYaxis()->SetTitleOffset(0.8);
	  frame->Draw();
	  allCaptions[icat]->Draw();
	}

      //likelihoods
      TPad *p = (TPad *)c->cd(ncategs+1);
      p->SetGridx();
      p->SetGridy();
      RooPlot *frame = topmass.frame(Bins(100),Range(tmin,tmax),Title("Likelihood")) ;    
      combll->plotOn(frame,ShiftToZero(),Name("comb"));
      for(int icat=0; icat<ncategs; icat++) 
	{
	  TString catTag("cat"); catTag += icat;
	  allLL[icat]->plotOn(frame,ShiftToZero(),LineColor(kGreen+4-2*icat),LineStyle(kDotted),LineWidth(1),Name(catTag));
	}
      frame->GetXaxis()->SetTitle("Top Mass [GeV/c^{2}]");
      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetYaxis()->SetTitle("-log(L/L_{max})");
      frame->GetYaxis()->SetTitleOffset(1);
      frame->Draw();

      TLegend *leg =  new TLegend(0.7,0.75,0.95,0.93,NULL,"brNDC");
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry("comb","Combined","l");
      for(int icat=0; icat<ncategs; icat++) 
	{
	  TString catTag("cat"); catTag += icat;
	  leg->AddEntry(catTag,allCategTitles[icat],"l");
	}	
      leg->Draw();
      formatForCmsPublic(p,leg,"CMS preliminary",3);
	  

      c = getNewCanvas("incmassfitter","Inclusive Fit Result",false);
      c->cd();
      c->Clear();
      c->SetWindowSize(500,500);
      c->SetGridy();
      c->SetGridx();
      frame = Mass.frame(Title("inclusivesample"));
      inclusiveData.plotOn(frame,Binning(15),DrawOption("pz"));

      RooAddPdf incbckgmassmodel("incbckg","Inclusive Background Model",bckgpdfset,constrParams);
      incbckgmassmodel.plotOn(frame,DrawOption("f"),FillStyle(3001),Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
      signalpdfset.add(bckgpdfset);
      sigYieldParams.add(constrParams);
      RooAddPdf incmassmodel("incmodel","Inclusive Model",signalpdfset,sigYieldParams);
      incmassmodel.plotOn(frame,LineStyle(kDotted),Normalization(1.0,RooAbsReal::RelativeExpected),MoveToBack());
      frame->GetXaxis()->SetTitle("Reconstructed Mass [GeV/c^{2}]");
      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetYaxis()->SetTitle("Events");
      frame->GetYaxis()->SetTitleOffset(0.8);
      frame->Draw();

      char buf[100];
      sprintf(buf,"CMS preliminary\\m_{top}=%3.1f #pm %3.1f GeV/c^{2}",topmass.getVal(),topmass.getError());
      formatForCmsPublic(p,leg,buf,3);
	  
      TPad *npad = new TPad("llpad","ll", 0.6, 0.6, 0.9, 0.9);
      npad->Draw();
      npad->cd();
      npad->SetGridx();
      npad->SetGridy();     
      frame = topmass.frame(Bins(100),Range(tmin,tmax),Title("Likelihood")) ;    
      combll->plotOn(frame,ShiftToZero(),Name("inccomb"));
      frame->GetXaxis()->SetTitle("Top Quark Mass [GeV/c^{2}]");
      frame->GetXaxis()->SetTitleOffset(0.8);
      frame->GetYaxis()->SetTitle("-log(L/L_{max})");
      frame->GetYaxis()->SetTitleOffset(1);
      frame->Draw();
    }
      
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
  if(cattype==0) fitPars_["cat0)Inclusive"]=0;
  if(cattype==1) 
    {
      fitPars_["cat0)Same flavor"]=0;
      fitPars_["cat1)Op. flavor"]=1;
    }
  if(cattype==2) 
    {
      fitPars_["cat0)= 0 b-tags"]=0;
      fitPars_["cat1)= 1 b-tags"]=1;
      fitPars_["cat2)#geq 2 b-tags"]=2;
    }

  return fitPars_;
}



