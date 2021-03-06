#include "LIP/Top/interface/MisassignmentMeasurement.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "CMGTools/HtoZZ2l2nu/interface/METUtils.h"

#include "TSystem.h"

using namespace std;

//
void MisassignmentMeasurement::bookMonitoringHistograms()
{
  //m_lj plots
  //  double massAxis[]={10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,400,500,1000,2000};
  double massAxis[]={10,20,30,40,50,60,70,80,90,100,115,130,145,160,185,200,250,300,400,500,1000,2000};

  int nMassBins=sizeof(massAxis)/sizeof(double)-1;
  TString spec[]={"","avg","data"};
  for(size_t is=0; is<3; is++)
    {
      controlHistos.addHistogram( new TH1D(spec[is]+"inclusivemlj","Lepton-jet spectrum;Invariant Mass [GeV/c^{2}];Lepton-jet pairs/#Delta M",nMassBins,massAxis) );
      controlHistos.addHistogram( new TH1D(spec[is]+"correctmlj","Correct assignments;Invariant Mass [GeV/c^{2}];Lepton-jet pairs/#Delta M",nMassBins,massAxis) );     
      controlHistos.addHistogram( new TH1D(spec[is]+"wrongmlj","Misassignments;Invariant Mass [GeV/c^{2}];Lepton-jet pairs/#Delta M",nMassBins,massAxis) );          
      controlHistos.addHistogram( new TH1D(spec[is]+"rotmlj","Random Rotation model;Invariant Mass [GeV/c^{2}];Lepton-jet pairs/#Delta M",nMassBins,massAxis) );        
      controlHistos.addHistogram( new TH1D(spec[is]+"swapmlj","Swap model;Invariant Mass [GeV/c^{2}];Lepton-jet pairs/#Delta M",nMassBins,massAxis) );                 
      controlHistos.addHistogram( new TH1D(spec[is]+"wrongmodelmlj","Missassignment model;Invariant Mass [GeV/c^{2}];Lepton-jet pairs/#Delta M",nMassBins,massAxis) ); 
    }
 
  //averaged plots only
  TH1D *h=new TH1D("jetflavor","Jet flavor;Jet flavor;Jets",3,0,3);
  h->GetXaxis()->SetBinLabel(1,"udsg");
  h->GetXaxis()->SetBinLabel(2,"c");
  h->GetXaxis()->SetBinLabel(3,"b");
  controlHistos.addHistogram( h );
  controlHistos.addHistogram( new TH1D("bias",";bias=#Delta f_{correct assignments};Pseudo-experiments",25,-0.52,0.48) );
  controlHistos.addHistogram( new TH1D("pull",";pull=#Delta f_{correct assignments} / #sigma_{stat};Pseudo-experiments",25,-5.2,4.8) );
  controlHistos.addHistogram( new TH1D("staterr",";#sigma_{stat}/f_{correct assignments};Pseudo-experiments",50,0,1) );
  controlHistos.addHistogram( new TH1D("fcorr",";f_{correct};Pseudo-experiments",200,0,0.5) );
  controlHistos.addHistogram( new TH1D("truefcorr",";f_{correct} (MC truth);Pseudo-experiments",200,0,0.5) );
  controlHistos.addHistogram( new TH1D("sigonlytruefcorr",";f_{correct} (MC truth);Pseudo-experiments",200,0,0.5) );
  controlHistos.addHistogram( new TH1D("fcorrsf",";SF=f_{correct}/f_{correct}^{M>y} (MC truth);Pseudo-experiments",200,0,1.2) );
  controlHistos.addHistogram( new TH1D("knorm",";Model k-factor;Pseudo-experiments",200,0,1.5) );

  //instantiate for different categories
//   controlHistos.initMonitorForStep("ee");
//   controlHistos.initMonitorForStep("emu");
//   controlHistos.initMonitorForStep("mumu");
//   controlHistos.initMonitorForStep("ll");

}

//
void MisassignmentMeasurement::resetHistograms(bool fullReset)
{
  /*
    //FIXME
  if(fullReset) nMeasurements=0;
  SmartSelectionMonitor::Monitor_t &mons=controlHistos.getAllMonitors();
  for(SmartSelectionMonitor::Monitor_t::iterator it=mons.begin(); it!=mons.end(); it++)
    {
      TString hname=it->second->GetName();
      if(!fullReset)
	{
	  if(hname.Contains("avg")) continue; 
	  if(hname.Contains("pull")) continue; 
	  if(hname.Contains("bias")) continue; 
	  if(hname.Contains("staterr")) continue; 
	  if(hname.Contains("fcorr")) continue; 
	  if(hname.Contains("truefcorr")) continue; 
	  if(hname.Contains("sigonlytruefcorr")) continue; 
	  if(hname.Contains("jetflavor")) continue;
	  if(hname.Contains("knorm")) continue;
	}
      it->second->Reset("ICE");
    }
  */
}

//
void MisassignmentMeasurement::finishMonitoringHistograms()
{
/*FIXME
 SelectionMonitor::StepMonitor_t &mons=controlHistos.getAllMonitors();
  for(SelectionMonitor::StepMonitor_t::iterator it=mons.begin(); it!=mons.end(); it++)
    {
      for(SelectionMonitor::Monitor_t::iterator hit=it->second.begin(); hit!= it->second.end(); hit++)
	{
	  fixExtremities(hit->second,true,true);
	  TString hname=hit->second->GetName();
	  bool doAverage( (hname.Contains("jetflavor")  || hname.Contains("knorm") 
			   || hname.Contains("fcorr")   || hname.Contains("avg") 
			   || hname.Contains("staterr") || hname.Contains("bias") 
			   || hname.Contains("pull")) && nMeasurements>0);  
	  bool doGausFit(hname.Contains("bias") || hname.Contains("pull"));
	  if(doAverage)  hit->second->Scale(1./nMeasurements); 
	  if(doGausFit)  hit->second->Fit("gaus","Q");
        }
    }
*/
}

//
void MisassignmentMeasurement::saveMonitoringHistograms()
{
  //open file
  TFile *fout=TFile::Open("MisassignmentMeasurement.root","RECREATE");
  TDirectory *baseOutDir=fout->mkdir("localAnalysis");
  baseOutDir->cd();
  controlHistos.Write();
  fout->Close();
  fout->Delete();
}


//
PhysicsObjectLeptonCollection randomlyRotate( PhysicsObjectLeptonCollection &leptons, PhysicsObjectJetCollection &jets, TRandom2 &rndGen)
{
  PhysicsObjectLeptonCollection rotLeptons;

  //create a rotated copy of all leptons
  int ileptons(0);
  for(PhysicsObjectLeptonCollection::iterator lit = leptons.begin(); lit != leptons.end(); lit++, ileptons++)
    {
      if(ileptons>=2) break;
      int itry=0;
      do
	{
	  itry++;
	  if(itry>1000) { cout << "Failed to rotate lepton:" << itry << endl;  rotLeptons.push_back( *lit ); break;  }

	  //rotate lepton
	  double en    = lit->E();
	  double pabs  = lit->P();
	  double phi   = rndGen.Uniform(0,2*TMath::Pi());
	  double theta = TMath::ACos( rndGen.Uniform(-1,1) );
	  LorentzVector rotLep =  LorentzVector( pabs*TMath::Cos(phi)*TMath::Sin(theta),
						 pabs*TMath::Sin(phi)*TMath::Sin(theta),
						 pabs*TMath::Cos(theta),
						 en);
	  
	  //require selectable kinematics
	  if( TMath::Abs(rotLep.Eta())>2.4 || rotLep.Pt()<20 ) continue;	  

	  //require object separation
	  double minDR(1000);
	  for(PhysicsObjectJetCollection::iterator jit = jets.begin(); jit != jets.end(); jit++)
	    {
	      double dR = deltaR(jit->eta(),jit->phi(),rotLep.eta(),rotLep.phi());
	      if(dR>minDR) continue;
	      minDR=dR;
	    }
	  
	  if(minDR>0.4)
	    {
	      //save lepton
	      rotLeptons.push_back(PhysicsObject_Lepton(rotLep,lit->id));
	      break;
	    }
	} while( 1 );
    }

  //all done
  return rotLeptons;
}

//
void MisassignmentMeasurement::measureMisassignments(ZZ2l2nuSummaryHandler &evHandler, double mcut, double minMlj, bool isData, int jetbin, TString syst)
{
  if(evHandler.getEntries()==0) return;
 
  if(!isData) nMeasurements++;
  resetHistograms();    
  std::map<TString,int> nCorrectAssignments, nCorrectAssignmentsFullSpectrum, totalEventsUsed;
  
  //run over the entries
  TTree *evTree=evHandler.getTree();
  unsigned int ntrials(10);
  for( unsigned int itrial=0; itrial<ntrials; itrial++)
    {
      for(unsigned int i=0; i<evTree->GetEntriesFast(); i++)
	{
	  evTree->GetEntry(i);

	  ZZ2l2nuSummary_t &ev = evHandler.getEvent();
	  int evcat=ev.cat;
	  //	  bool hasSignal( ev.evmeasurements[0]>0 );

	  //save local event
	  PhysicsEvent_t phys = getPhysicsEventFrom(ev);
	  if((syst=="metcut") && ((ev.cat==EE || ev.cat==MUMU) && phys.met[0].pt()<40)) continue;

	  PhysicsObjectLeptonCollection ileptons = phys.leptons;

	  PhysicsObjectJetCollection jetColl;
	  std::vector<size_t> origJetsIdx;
	  for(size_t ijet=0; ijet<phys.jets.size(); ijet++)
	    {
	      if(phys.jets[ijet].pt()<30 || fabs(phys.jets[ijet].eta())>2.5) continue;
	      jetColl.push_back(phys.jets[ijet]);
	      origJetsIdx.push_back(ijet);
	    }
	  
	  //vary and reset if required a systematic variation
	  if(syst!="" && syst !="metcut")
	    {
	      std::vector<PhysicsObjectJetCollection> jetsVar;
	      LorentzVectorCollection metsVar;
	      METUtils::computeVariation(phys.jets, phys.leptons, phys.met[0], jetsVar, metsVar, jecUnc_);
	      if(syst.Contains("jesup"))    jetColl = jetsVar[METUtils::JES_UP];
	      if(syst.Contains("jesdown"))  jetColl = jetsVar[METUtils::JES_DOWN];
	      if(syst.Contains("jer"))      jetColl = jetsVar[METUtils::JER];
	    }
	     	     
	  if(jetbin!=0 && int(jetColl.size())!=jetbin) continue;

	  std::vector<TString> categs;
	  categs.push_back("all");
	  if(evcat==MUMU) { categs.push_back("mumu"); categs.push_back("ll"); }
	  if(evcat==EE)   { categs.push_back("ee");   categs.push_back("ll"); }
	  if(evcat==EMU)  { categs.push_back("emu"); }
	  if(itrial==0)
	    {
	      for(size_t icat=0; icat<categs.size(); icat++)
		{
		  if(totalEventsUsed.find(categs[icat]) == totalEventsUsed.end())  totalEventsUsed[categs[icat]]=0;
		  totalEventsUsed[categs[icat]]++;
		}
	    }
      	  
	  //
	  // MODEL 1 get event mixed jets in equal number to current event's jet multiplicity
	  //
	  LorentzVectorCollection mixjets;
	  do{
	    int imixtry(0);
	    unsigned int j=rndGen.Uniform(0,evTree->GetEntriesFast());
	    if(j==i) continue;
	    imixtry++;
	    
	    if(imixtry>50) 
	      { 
		cout << "Failed to mix 1 event" << endl;
		cout << imixtry << " " << j << " " << i << " " << itrial << endl;
		continue; 
	      }

	    evTree->GetEntry(j);
	    ZZ2l2nuSummary_t &mixev = evHandler.getEvent();
	    //if(evcat!= mixev.cat) continue;
	    
	    //get the jets
	    PhysicsEvent_t mixphys = getPhysicsEventFrom(mixev);
	    PhysicsObjectJetCollection mixjetColl=mixphys.jets;
	    //vary and reset if required a systematic variation
	    if(syst!="" && syst !="metcut")
	      {
		std::vector<PhysicsObjectJetCollection> mixjetsVar;
		LorentzVectorCollection mixmetsVar;
		METUtils::computeVariation(mixphys.jets, mixphys.leptons, mixphys.met[0], mixjetsVar, mixmetsVar, jecUnc_);
		if(syst.Contains("jesup"))    mixjetColl = mixjetsVar[METUtils::JES_UP];
		if(syst.Contains("jesdown"))  mixjetColl = mixjetsVar[METUtils::JES_DOWN];
		if(syst.Contains("jer"))      mixjetColl = mixjetsVar[METUtils::JER];
	      }


	    //now check kinematics and the separation with respect to the leptons in the event
	    for(size_t ijet=0; ijet< mixjetColl.size(); ijet++)
	      {
		if(mixjetColl[ijet].pt()<30 || fabs(mixjetColl[ijet].eta())>2.5) continue;
		
		double minDR(1000);
		for(PhysicsObjectLeptonCollection::iterator lit = ileptons.begin(); lit!=ileptons.end(); lit++)
		  {
		    double dR = deltaR(lit->eta(),lit->phi(),mixjetColl[ijet].eta(),mixjetColl[ijet].phi());
		    if(dR>minDR) continue;
		    minDR=dR;
		  }
		if(minDR<0.4) continue;
		
		//save jets
		if(mixjets.size()<jetColl.size()) mixjets.push_back( mixjetColl[ijet] );
	      }
	    
	    //continue until jet multiplicity is filled
	    if(mixjets.size()<jetColl.size()) continue;	
	    
	    break;
	  }while(1);
	  
	  //
	  // MODEL 2 get rotated leptons
	  //
	  PhysicsObjectLeptonCollection rotLeptons = randomlyRotate(ileptons,phys.jets,rndGen);
	  
	  //
	  // Fill the control histograms
	  //
	  //control jet flavor
	  if(itrial==0)
	    {
	      for(size_t ijet=0; ijet<jetColl.size(); ijet++)
		{
		  int flavbin=phys.jets[ origJetsIdx[ijet] ].flavid;
		  if(fabs(flavbin)==5)      flavbin=2;
		  else if(fabs(flavbin)==4) flavbin=1;
		  else                      flavbin=0;
		  for(size_t icateg=0; icateg<categs.size(); icateg++)
		    {
		      TString ctf=categs[icateg];
		      controlHistos.fillHisto("jetflavor",ctf,flavbin,1);
		    }
		}
	    }
      
	  //invariant mass spectrum
	  for(PhysicsObjectLeptonCollection::iterator lit = ileptons.begin(); lit != ileptons.end(); lit++)
	    {

	      //std pairs
	      if(itrial==0)
		{
		  for(size_t ijet=0; ijet<jetColl.size(); ijet++)
		    {
		      //the lepton-jet system
		      LorentzVector sum = *lit + jetColl[ijet];
		      double mlj=sum.M();
		  
		      //control if assignment is correct
		      int partonMatch=phys.jets[ origJetsIdx[ijet] ].genid;
		      int flavorMatch=phys.jets[ origJetsIdx[ijet] ].flavid;
		      int assignCode=(lit->genid*partonMatch);

		      //FIXME
		      bool isCorrect(assignCode<0  && fabs(flavorMatch)==5 );
		      //		      bool isCorrect( (assignCode<0) && hasSignal && fabs(flavorMatch)==5 );
		      for(size_t icateg=0; icateg<categs.size(); icateg++)
			{
			  TString ctf=categs[icateg];
			  if(nCorrectAssignments.find(ctf)==nCorrectAssignments.end())  nCorrectAssignments[ctf]=0;
			  if(mlj>minMlj)  			                        nCorrectAssignments[ctf] += isCorrect;
			  nCorrectAssignmentsFullSpectrum[ctf] += isCorrect;
		      
			  controlHistos.fillHisto("inclusivemlj",ctf,mlj,1,true);
			  controlHistos.fillHisto(isCorrect ? "correctmlj" : "wrongmlj",ctf,mlj,1,true);
			}
		    }
		}
	  
	      //event mixing spectrum
	      for(LorentzVectorCollection::iterator jit = mixjets.begin(); jit != mixjets.end(); jit++)
		{
		  //mixed lepton-jet system
		  LorentzVector sum = *lit + *jit;
		  double mlj=sum.M();
		  for(size_t icateg=0; icateg<categs.size(); icateg++)
		    {
		      TString ctf=categs[icateg];
		      controlHistos.fillHisto("swapmlj",ctf,mlj,1./ntrials,true);
		    }	  
		}
	    }

	  //random rotation
	  for(PhysicsObjectLeptonCollection::iterator lit = rotLeptons.begin(); lit != rotLeptons.end(); lit++)
	    {
	      for(PhysicsObjectJetCollection::iterator jit = jetColl.begin(); jit != jetColl.end(); jit++)
		{
		  //randomly rotated lepton-jet system
		  LorentzVector sum = *lit + *jit;
		  double mlj=sum.M();
		  for(size_t icateg=0; icateg<categs.size(); icateg++)
		    {
		      TString ctf=categs[icateg];
		      controlHistos.fillHisto("rotmlj",ctf,mlj,1./ntrials,true);
		    }
		}
	    }
	}
    }

  //
  // finalize model and estimate the misassignments
  //
  /*FIXME
    SelectionMonitor::StepMonitor_t &mons=controlHistos.getAllMonitors();
  for(SelectionMonitor::StepMonitor_t::iterator it=mons.begin(); it!=mons.end(); it++)
    {
      TString ctf = it->first;

      //the model (average individual estimates)
      TH1 *mljWrongModelH=controlHistos.getHisto("wrongmodelmlj",ctf);
      mljWrongModelH->Add( controlHistos.getHisto("swapmlj",ctf), 0.5);
      mljWrongModelH->Add( controlHistos.getHisto("rotmlj",ctf), 0.5);

      //inclusive spectrum
      TH1 *mljInclusiveH=controlHistos.getHisto("inclusivemlj", ctf);
      TH1 *mljCorrectH=controlHistos.getHisto("correctmlj", ctf);
      TH1 *mljWrongH=controlHistos.getHisto("wrongmlj", ctf);
      
      //estimate of alpha from the model
      double massCut = mcut;
      int nBins  = mljWrongModelH->GetXaxis()->GetNbins();
      double nPairsTotal(0), nPairsAboveCut(0),nPairsBelowCut(0), nPairsTotalModel(0), nPairsModelAboveCut(0), nPairsTotalSigOnly(0);
      double nPairsTotalErr(0), nPairsAboveCutErr(0), nPairsBelowCutErr(0), nPairsTotalModelErr(0), nPairsModelAboveCutErr(0), nPairsTotalSigOnlyErr(0);
      for(int ibin=1; ibin<=nBins; ibin++)
	{
	  double imass=mljWrongModelH->GetXaxis()->GetBinLowEdge(ibin);
	  double dM=mljWrongModelH->GetXaxis()->GetBinWidth(ibin);
	  if(imass<minMlj) 
	    {
	      nPairsBelowCut    +=  mljInclusiveH->GetBinContent(ibin)*dM;
	      nPairsBelowCutErr +=  pow(mljInclusiveH->GetBinError(ibin)*dM,2);
	      continue;
	    }
	      
	  nPairsTotal         += mljInclusiveH->GetBinContent(ibin)*dM;
	  nPairsTotalErr      += pow(mljInclusiveH->GetBinError(ibin)*dM,2);
	  nPairsTotalModel    += mljWrongModelH->GetBinContent(ibin)*dM;
	  nPairsTotalModelErr += pow(mljWrongModelH->GetBinError(ibin)*dM,2);
	  nPairsTotalSigOnly  += (mljCorrectH->GetBinContent(ibin)+mljWrongH->GetBinContent(ibin))*dM;
	  nPairsTotalSigOnlyErr  += pow(mljCorrectH->GetBinError(ibin)*dM,2)+pow(mljWrongH->GetBinError(ibin)*dM,2);
	  
	  if(imass<massCut) continue;
	  nPairsAboveCut         += mljInclusiveH->GetBinContent(ibin)*dM;
	  nPairsAboveCutErr      += pow(mljInclusiveH->GetBinError(ibin)*dM,2);
	  nPairsModelAboveCut    += mljWrongModelH->GetBinContent(ibin)*dM;
	  nPairsModelAboveCutErr += pow(mljWrongModelH->GetBinError(ibin)*dM,2);
	}
      nPairsTotalErr         = sqrt(nPairsTotalErr);
      nPairsBelowCutErr      = sqrt(nPairsBelowCutErr);
      nPairsAboveCutErr      = sqrt(nPairsAboveCutErr);
      nPairsTotalModelErr    = sqrt(nPairsTotalModelErr);
      nPairsModelAboveCutErr = sqrt(nPairsModelAboveCutErr);
      nPairsTotalSigOnlyErr  = sqrt(nPairsTotalSigOnlyErr);
            
      double a = nPairsTotal;
      double b = totalEventsUsed[ctf];
      double c = nPairsTotalModel;
      double d = nPairsAboveCut;
      double e = nPairsModelAboveCut;
	  
      //re-scale misassignment models
      kNorm[ctf] = c*d/(a*e);
      mljWrongModelH->Scale(kNorm[ctf]);
      controlHistos.getHisto("swapmlj",ctf)->Scale(kNorm[ctf]);
      controlHistos.getHisto("rotmlj",ctf)->Scale(kNorm[ctf]);

      //determine corect pairs fraction
      double nCorrectPairsEst = a-c*kNorm[ctf];
      double nCorrectPairsEstErr = sqrt( pow(1-c*c*d/(a*a*e),2)*a
					 + pow(-2*c*d/(a*e),2)*c
					 + pow(-c*c/(a*e),2)*d
					 + pow(c*c/(a*e*e),2)*e )/sqrt(2.);      
      
      alphaEst[ctf]            = nCorrectPairsEst/(2*b);
      alphaEstErr[ctf]         = nCorrectPairsEstErr/(2*b);
      //      fCorrectPairsEst[ctf]    = nCorrectPairsEst/nPairsTotal;//-bias[ctf]; //debug me
      // fCorrectPairsEstErr[ctf] = nCorrectPairsEstErr/nPairsTotal;
      fCorrectPairsEst[ctf]    = 1.-nPairsAboveCut/nPairsModelAboveCut-bias[ctf];
      fCorrectPairsEstErr[ctf] = sqrt(pow(nPairsAboveCut*nPairsModelAboveCutErr,2)+pow(nPairsModelAboveCut*nPairsBelowCutErr,2))/pow(nPairsModelAboveCut,2);
      double sfmc              = double(nCorrectAssignmentsFullSpectrum[ctf])/double(nCorrectAssignments[ctf]);
      
      //average models for posterity
      TString prefix( isData ? "data" : "avg" );
      controlHistos.getHisto( prefix+"inclusivemlj",ctf )->Add( mljInclusiveH );
      controlHistos.getHisto( prefix+"correctmlj",ctf )->Add( controlHistos.getHisto( "correctmlj", ctf) );
      controlHistos.getHisto( prefix+"wrongmlj",ctf )->Add( controlHistos.getHisto( "wrongmlj", ctf) );
      controlHistos.getHisto( prefix+"rotmlj",ctf )->Add( controlHistos.getHisto( "rotmlj", ctf) );
      controlHistos.getHisto( prefix+"swapmlj", ctf )->Add( controlHistos.getHisto( "swapmlj", ctf) );
      controlHistos.getHisto( prefix+"wrongmodelmlj", ctf )->Add( mljWrongModelH );

      //estimate monitoring distributions
      if(!isData)
	{
	  double fCorrectMCTruth=double(nCorrectAssignments[ctf])/double(nPairsTotal);
	  double fCorrectMCTruthSigOnly=double(nCorrectAssignments[ctf])/double(nPairsTotalSigOnly);
	  controlHistos.fillHisto("truefcorr", ctf, fCorrectMCTruth,1);
	  controlHistos.fillHisto("sigonlytruefcorr", ctf, fCorrectMCTruthSigOnly,1);
	  controlHistos.fillHisto("fcorrsf", ctf, sfmc,1);
	  controlHistos.fillHisto("fcorr",ctf,fCorrectPairsEst[ctf],1);
	  controlHistos.fillHisto("knorm",ctf,kNorm[ctf],1);
	  controlHistos.fillHisto("bias",ctf, fCorrectPairsEst[ctf]-fCorrectMCTruth,1);
	  controlHistos.fillHisto("pull",ctf, (fCorrectPairsEst[ctf]-fCorrectMCTruth)/fCorrectPairsEstErr[ctf],1);
	  controlHistos.fillHisto("staterr",ctf, fCorrectPairsEstErr[ctf]/fCorrectPairsEst[ctf],1);
	}
    }
  */
}


//
void MisassignmentMeasurement::fitCurrentModelsToData()
{
  
  TString cats[]={"ee","mumu","emu"};
  std::vector<Double_t> minll;
  for(size_t icat=0; icat<sizeof(cats)/sizeof(TString); icat++)
    {
      TString ctf=cats[icat];

      TH1F *dataH=(TH1F *)controlHistos.getHisto( "datainclusivemlj",ctf )->Clone("datatmp");
      dataH->Add( (TH1F *) controlHistos.getHisto("datawrongmodelmlj",ctf) , -1);

      TH1F *mcH=(TH1F *)controlHistos.getHisto( "avgcorrectmlj",ctf )->Clone("tmp");
      mcH->Scale( dataH->Integral(1,dataH->GetXaxis()->GetNbins()+1,"width")/mcH->Integral(1,dataH->GetXaxis()->GetNbins()+1,"width") ); 
      
      double nll(0);
      for(int ibin=1; ibin<=dataH->GetXaxis()->GetNbins(); ibin++)
	{
	  double x=dataH->GetBinContent(ibin);
	  double y=mcH->GetBinContent(ibin);
	  double width=dataH->GetBinWidth(ibin);
	  double prob=TMath::Poisson(x*width,y*width);
	  //cout << x << " " << y << " " << prob << endl;
	  if(prob>0) nll += -TMath::Log(prob);
	}
      minll.push_back(nll);
      
      delete dataH;
      delete mcH;
    }
  
  Double_t totalMinll(0);
  cout << "Likelihood results" << endl;
  for(size_t icat=0; icat<sizeof(cats)/sizeof(TString); icat++)
    {
      cout << cats[icat] << " " << minll[icat] << endl;
      totalMinll += minll[icat];
    }
  cout << "Inclusive: " << totalMinll << endl;

}
