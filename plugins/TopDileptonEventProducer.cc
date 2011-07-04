#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"
#include "DataFormats/PatCandidates/interface/EventHypothesisLooper.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "CMGTools/HtoZZ2l2nu/interface/TSelectionMonitor.h"

#include "LIP/Top/interface/GenTopEvent.h"

#include "TH1D.h"
#include "TString.h"
 

class TopDileptonEventProducer : public edm::EDProducer {
public:
  TopDileptonEventProducer(const edm::ParameterSet &iConfig) ;
  virtual void produce( edm::Event &iEvent, const edm::EventSetup &iSetup) ;
private:
  std::map<std::string, edm::ParameterSet > objConfig;
  TSelectionMonitor controlHistos_;
  gen::top::Event genEvent_;
};


using namespace std;


//
TopDileptonEventProducer::TopDileptonEventProducer(const edm::ParameterSet &iConfig)
  : controlHistos_("top") 
{
  produces<std::vector<pat::EventHypothesis> >("selectedEvent");
  produces<reco::VertexCollection>("selectedVertices");
  produces<std::vector<int> >("selectionInfo");
  std::string objs[]={"Generator", "Vertices", "Electrons", "Muons", "Dileptons", "Jets", "MET" };
  for(size_t iobj=0; iobj<sizeof(objs)/sizeof(string); iobj++)
    objConfig[ objs[iobj] ] = iConfig.getParameter<edm::ParameterSet>( objs[iobj] );

  //create the control histograms
  controlHistos_.addHistogram("ngenpu", "; Generated pileup; Events", 25, 0.,25.);
  controlHistos_.addHistogram("vertexmult", "; Reconstructed vertices; Events", 100, 0.,25.);
  controlHistos_.addHistogram("rho", "; #rho; Events", 100, 0.,10.);
  controlHistos_.addHistogram("ecaliso", ";Photon isolation; Events", 100, 0.,10.);
  controlHistos_.addHistogram("hcaliso", ";Neutral hadron isolation; Events", 100, 0.,10.);
  controlHistos_.addHistogram("caloiso", ";Charged+neutral hadron isolation; Events", 100, 0.,10.);
  controlHistos_.addHistogram("trackiso", ";Charged hadron isolation; Events", 100, 0.,10.);
  controlHistos_.addHistogram("reliso", ";Relative isolation; Events", 100, 0.,10.);
  TString subcat[]={"eb","ee"};
  for(size_t isubcat=0; isubcat<2; isubcat++)
    {
      controlHistos_.addHistogram(subcat[isubcat]+"eoverp",";E/p;Electrons",100,0,5);
      controlHistos_.addHistogram(subcat[isubcat]+"einoverp",";E_{seed}/p;Electrons",100,0,5);
      controlHistos_.addHistogram(subcat[isubcat]+"fbrem",";f_{bremm};Electrons",100,0,1.1);
      controlHistos_.addHistogram(subcat[isubcat]+"hovere",";H/E;Electrons",100,0,0.5);
      controlHistos_.addHistogram(subcat[isubcat]+"sigmaietaieta",";#sigma_{i-#eta,i-#eta};Electrons",100,0,0.04);
      controlHistos_.addHistogram(subcat[isubcat]+"dphiin",";#Delta#phi(SC,track_{in});Electrons",100,-0.1,0.1);
      controlHistos_.addHistogram(subcat[isubcat]+"detain",";#Delta#eta(SC,track_{in});Electrons",100,-0.1,0.1);
      controlHistos_.addHistogram(subcat[isubcat]+"misshits",";N_{miss hits};Electrons",5,0,5);
    }
  
  controlHistos_.addHistogram("jetmult", "; Jet multiplicity; Events", 100, 0.,10.);
  controlHistos_.addHistogram("jetchconst", "; Charged hadron multiplicity in jets; Jets", 25, 0.,25.);
  controlHistos_.addHistogram("jetneutconst", "; Neutral hadron multiplicity in jets; Jets", 25, 0.,25.);
  controlHistos_.addHistogram("jetpt", "; p_{T} [GeV/c]; Jet", 100, 0.,250.);
  controlHistos_.addHistogram("jetbeta", "; #beta(dilepton vertex); Jets", 100, 0.,1.);
 
  controlHistos_.initMonitorForStep("electron");
  controlHistos_.initMonitorForStep("muon");
  controlHistos_.initMonitorForStep("ee");
  controlHistos_.initMonitorForStep("emu");
  controlHistos_.initMonitorForStep("mumu");
}

//
void TopDileptonEventProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
  using namespace std;
  using namespace edm;
  using namespace pat::eventhypothesis;
  using reco::Candidate; 
  using reco::CandidatePtr;
  
  pat::EventHypothesis hyp;
  int selStep(0),selPath(0);
  
  //get the weight for the MC event
  float weight=1;
  if(!iEvent.isRealData())
    {
      edm::Handle<float> puWeightHandle;
      iEvent.getByLabel("puWeights","puWeight",puWeightHandle);
      if(puWeightHandle.isValid()) weight = *(puWeightHandle.product());

      int npuIT(0); 
      edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
      iEvent.getByType(puInfoH);
      if(puInfoH.isValid())
	{
	  for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++)
	    if(it->getBunchCrossing()==0) npuIT += it->getPU_NumInteractions();
	}
      controlHistos_.fillHisto("ngenpu","all",npuIT,weight);
    }
  
  //pre-select vertices
  Handle<reco::VertexCollection> hVtx;
  iEvent.getByLabel(objConfig["Vertices"].getParameter<edm::InputTag>("source"), hVtx);  
  std::vector<reco::VertexRef> selVertices = vertex::filter(hVtx,objConfig["Vertices"]);
  controlHistos_.fillHisto("vertexmult","all",selVertices.size(),weight);
  std::vector<reco::VertexRef> primaryVertexHyps;
  if(selVertices.size()>0) selStep=1;

  //beam spot
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel( objConfig["Vertices"].getParameter<edm::InputTag>("beamSpot"), beamSpot);

  //average energy density
  edm::Handle< double > rho;
  iEvent.getByLabel(objConfig["Jets"].getParameter<edm::InputTag>("rho"),rho);

  //select muons (id+very loose isolation)
  Handle<View<Candidate> > hMu; 
  iEvent.getByLabel(objConfig["Muons"].getParameter<edm::InputTag>("source"), hMu);
  CandidateWithVertexCollection selMuons = muon::filter(hMu, selVertices, *beamSpot, objConfig["Muons"]);

  //select electrons (id+conversion veto+very loose isolation)
  Handle<View<Candidate> > hEle; 
  iEvent.getByLabel(objConfig["Electrons"].getParameter<edm::InputTag>("source"), hEle);
  CandidateWithVertexCollection selElectrons = electron::filter(hEle, hMu, selVertices, *beamSpot, objConfig["Electrons"]);
  for(View<Candidate>::const_iterator eIt=hEle->begin(); eIt != hEle->end(); eIt++)
    {
      const pat::Electron &electron=dynamic_cast<const pat::Electron &>(*eIt);
      TString cat(electron.isEB() ? "eb":"ee");
      controlHistos_.fillHisto(cat+"eoverp","electron",electron.eSuperClusterOverP(),weight);
      controlHistos_.fillHisto(cat+"einoverp","electron", electron.eSeedClusterOverP(),weight);
      controlHistos_.fillHisto(cat+"fbrem","electron",electron.fbrem(),weight);
      controlHistos_.fillHisto(cat+"hovere","electron",electron.hadronicOverEm(),weight);
      controlHistos_.fillHisto(cat+"sigmaietaieta","electron",electron.sigmaIetaIeta(),weight); 
      controlHistos_.fillHisto(cat+"dphiin","electron",electron.deltaPhiSuperClusterTrackAtVtx(),weight);
      controlHistos_.fillHisto(cat+"detain","electron",electron.deltaEtaSuperClusterTrackAtVtx(),weight);
      if(electron.gsfTrack().isNull()) continue;
      controlHistos_.fillHisto("electron_"+cat+"misshits","electron",electron.gsfTrack()->trackerExpectedHitsInner().numberOfHits(),weight);
    }

  //build inclusive collection
  CandidateWithVertexCollection selLeptons = selMuons;
  selLeptons.insert(selLeptons.end(), selElectrons.begin(), selElectrons.end());
  if(selLeptons.size()>0) selStep=2;

  //if event is MC filter out the genparticle collection also
  int tteventcode=gen::top::Event::UNKNOWN;
  if(!iEvent.isRealData())
    {
      genEvent_.genLabel_=objConfig["Generator"].getParameter<edm::InputTag>("source");
      tteventcode = genEvent_.assignTTEvent(iEvent,iSetup);
      
      if(objConfig["Generator"].getParameter<bool>("filterSignal") 
	 && tteventcode!=gen::top::Event::EE 
	 && tteventcode != gen::top::Event::EMU 
	 && tteventcode!= gen::top::Event::MUMU) 
	selStep=0;

      //save the generator level event
      std::map<std::string, std::list<reco::CandidatePtr> > genParticles;
      genParticles["top"] = genEvent_.tops;
      genParticles["quarks"] = genEvent_.quarks;
      genParticles["leptons"] = genEvent_.leptons;
      genParticles["neutrinos"]  = genEvent_.neutrinos;
      for(std::map<std::string,std::list<reco::CandidatePtr> >::iterator it = genParticles.begin();
	  it != genParticles.end(); it++)
	{
       	  for(std::list<reco::CandidatePtr>::iterator itt = it->second.begin();
       	      itt != it->second.end();
       	      itt++)
	    hyp.add( *itt, it->first );
       	}
    }
  

  //build the dilepton (all tightly isolated leptons will be returned)
  if(selStep==2)
    {
      //control histos for leptons
      for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
	{
	  using namespace lepton;
	  int id = getLeptonId(selLeptons[ilep].first);
	  std::vector<double> isol=getLeptonIso(selLeptons[ilep].first,objConfig["Dileptons"].getParameter<double>("minPt"));
	  TString ptype(fabs(id)==ELECTRON ? "electron" : "muon");
	  controlHistos_.fillHisto("rho",ptype,*rho,weight);
	  controlHistos_.fillHisto("ecaliso",ptype,isol[ECAL_ISO],weight);
	  controlHistos_.fillHisto("hcaliso",ptype,isol[HCAL_ISO],weight);
	  controlHistos_.fillHisto("trackiso",ptype,isol[TRACKER_ISO],weight);
	  controlHistos_.fillHisto("caloiso",ptype,isol[ECAL_ISO]+isol[HCAL_ISO]+isol[TRACKER_ISO],weight);
	  controlHistos_.fillHisto("reliso",ptype,isol[REL_ISO],weight);
	}
      
      //search for dileptons
      CandidateWithVertexCollection isolLeptons;
      std::pair<CandidateWithVertex,CandidateWithVertex> dileptonWithVertex = dilepton::filter(selLeptons,
											       isolLeptons,
											       *rho,
											       objConfig["Dileptons"],
											       iSetup);
      selPath = dilepton::classify(dileptonWithVertex);
      if(selPath!= dilepton::UNKNOWN)
	{
	  std::vector<TString> dilCats;
	  dilCats.push_back("all"); 
	  if(selPath==dilepton::EE)   dilCats.push_back("ee"); 
	  if(selPath==dilepton::EMU)  dilCats.push_back("emu"); 
	  if(selPath==dilepton::MUMU) dilCats.push_back("mumu"); 

	  //add to the event
	  selStep=3;
	  hyp.add(dileptonWithVertex.first.first,"leg1");
	  primaryVertexHyps.push_back(dileptonWithVertex.first.second);
	  hyp.add(dileptonWithVertex.second.first,"leg2");
	  primaryVertexHyps.push_back(dileptonWithVertex.second.second);

	  //add the remaining isolated leptons now
	  for(CandidateWithVertexCollection::iterator lIt = isolLeptons.begin(); lIt != isolLeptons.end(); lIt++)
	    {
	      if(lIt->first.get()== dileptonWithVertex.first.first.get()) continue;
	      if(lIt->first.get()== dileptonWithVertex.second.first.get()) continue;

	      //check if lepton is associated to the same vertex as the dilepton
	      if(lIt->second.get()== dileptonWithVertex.first.second.get() 
		 || lIt->second.get() == dileptonWithVertex.second.second.get())
		hyp.add( lIt->first , fabs(lepton::getLeptonId(lIt->first))==lepton::MUON ? "muon" : "electron" );
	      else
		hyp.add( lIt->first , fabs(lepton::getLeptonId(lIt->first))==lepton::MUON ? "pumuon" : "puelectron" );
	      
	    }
	  
	  //add also the jets
	  Handle<View<Candidate> > hJet; 
	  iEvent.getByLabel(objConfig["Jets"].getParameter<edm::InputTag>("source"), hJet);
	  CandidateWithVertexCollection selJets = jet::filter(hJet, isolLeptons, selVertices, objConfig["Jets"]);
	  for(size_t icat=0; icat<dilCats.size(); icat++)
	    {
	      controlHistos_.fillHisto("jetmult",dilCats[icat],selJets.size(),weight);
	      for(CandidateWithVertexCollection::iterator jit = selJets.begin(); jit!=selJets.end(); jit++)
		{
		  const pat::Jet *j = dynamic_cast<const pat::Jet *>(jit->first.get());
		  controlHistos_.fillHisto("jetpt",dilCats[icat],j->pt(),weight);
		  controlHistos_.fillHisto("jetchconst",dilCats[icat],j->chargedHadronMultiplicity(),weight);
		  controlHistos_.fillHisto("jetneutconst",dilCats[icat],j->neutralHadronMultiplicity(),weight);
		  if(dileptonWithVertex.first.second.isNonnull())
		    controlHistos_.fillHisto("jetbeta",dilCats[icat],jet::fAssoc(j,dileptonWithVertex.first.second.get()),weight);
		}
	    }
	  CandidateWithVertexCollection assocJets, puJets;
	  jet::classifyJetsForDileptonEvent(selJets,dileptonWithVertex,assocJets,puJets,objConfig["Dileptons"].getParameter<double>("maxDz"));
	  for(CandidateWithVertexCollection::iterator jIt = assocJets.begin(); jIt != assocJets.end(); jIt++) hyp.add(jIt->first,"jet");
	  for(CandidateWithVertexCollection::iterator jIt = puJets.begin(); jIt != puJets.end(); jIt++) hyp.add(jIt->first,"pujet");
	  
	  //add the met
	  Handle<View<Candidate> > hMET; 
	  iEvent.getByLabel(objConfig["MET"].getParameter<edm::InputTag>("source"), hMET);
	  CandidatePtr met = hMET->ptrAt(0);
	  hyp.add(met, "met");
	}
    }
     
  // work done, save results
  auto_ptr<std::vector<pat::EventHypothesis> > hyps(new std::vector<pat::EventHypothesis>() );
  hyps->push_back(hyp);
  iEvent.put(hyps,"selectedEvent");

  auto_ptr<std::vector<int> > selectionInfo(new std::vector<int>());
  selectionInfo->push_back( selPath );
  selectionInfo->push_back( selStep );
  selectionInfo->push_back( tteventcode );
  iEvent.put(selectionInfo,"selectionInfo");

  auto_ptr<reco::VertexCollection> primVertices(new reco::VertexCollection() );
  for(std::vector<reco::VertexRef>::iterator vit=primaryVertexHyps.begin(); vit != primaryVertexHyps.end(); vit++)
    primVertices->push_back( *(vit->get()) );
  iEvent.put(primVertices,"selectedVertices");
}

//  
DEFINE_FWK_MODULE(TopDileptonEventProducer);
