#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "LIP/Top/interface/EventSummaryHandler.h"

#include <vector>
#include <map>
#include <string>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TString.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"
#include "DataFormats/PatCandidates/interface/EventHypothesisLooper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CMGTools/HtoZZ2l2nu/interface/TSelectionMonitor.h"
#include "LIP/Top/interface/GenTopEvent.h"

#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

using namespace std;
using namespace reco;
using namespace top;

class TopDileptonEventAnalyzer : public edm::EDAnalyzer 
{
public:
  explicit TopDileptonEventAnalyzer(const edm::ParameterSet& cfg);
  ~TopDileptonEventAnalyzer(){};
  virtual void analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup) ;
  void endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup);

private:

  void saveEvent(const edm::Event& event, int evCat, std::vector<reco::CandidatePtr> &leptons, std::vector<const pat::Jet *> &jets, const pat::MET *met, 
		 int nvertices, int npuIT, float rho, float weight, float normWeight);

  std::map<std::string, edm::ParameterSet> objConfig_;
  
  EventSummaryHandler summaryHandler_;
  TSelectionMonitor controlHistos_;
  gen::top::Event genEvent_;
};

using namespace std;

/// default constructor
TopDileptonEventAnalyzer::TopDileptonEventAnalyzer(const edm::ParameterSet& cfg)
  : controlHistos_("top")
{
  try{

    //summary tree 
    edm::Service<TFileService> fs;
    summaryHandler_.initTree( fs->make<TTree>("data","Event Summary") );

    //configure the selection
    std::string objs[]={"Generator", "Trigger", "Vertices", "Electrons", "LooseElectrons", "Muons", "LooseMuons", "Dileptons", "Jets", "MET" };
    for(size_t iobj=0; iobj<sizeof(objs)/sizeof(string); iobj++)
      objConfig_[ objs[iobj] ] = cfg.getParameter<edm::ParameterSet>( objs[iobj] );

    //monitoring histograms
    TString selSteps[]={"Reco","2 leptons","M_{ll}","#geq 2 jets","MET>30,0","OS","=0 b-tags","=1 b-tags", "#geq 2 b-tags"};
    const size_t nselsteps=sizeof(selSteps)/sizeof(TString);
    controlHistos_.addHistogram("cutflow", ";Step; Events",nselsteps,0,nselsteps);
    for(int ibin=1; ibin<=controlHistos_.getHisto("cutflow","all")->GetXaxis()->GetNbins(); ibin++)
      controlHistos_.getHisto("cutflow","all")->GetXaxis()->SetBinLabel(ibin,selSteps[ibin-1]);
    
    //vertex control
    controlHistos_.addHistogram("ngoodvertex", ";Vertices; Events", 25, 0.,25.); 
    controlHistos_.addHistogram("vertex_sumpt", ";#Sigma_{tracks} p_{T} [GeV/c]; Events", 100, 0.,300.);
    controlHistos_.addHistogram("othervertex_sumpt", ";#Sigma_{tracks} p_{T} [GeV/c]; Events", 100, 0.,300.);
    controlHistos_.addHistogram("vertex_ndof", ";NDOF; Events", 100, 0.,100.);
    controlHistos_.addHistogram("othervertex_ndof", ";NDOF; Events", 100, 0.,100.);

    //lepton control 
    controlHistos_.addHistogram("egammaiso", "; Isolation; Events", 100, 0.,0.5);
    controlHistos_.addHistogram("echhadroniso", "; Isolation; Events", 100, 0.,0.5);
    controlHistos_.addHistogram("eneuhadroniso", "; Isolation; Events", 100, 0.,0.5);
    controlHistos_.addHistogram("ereliso", "; Isolation; Events", 100, 0.,1.);
    controlHistos_.addHistogram("mugammaiso", "; Isolation; Events", 100, 0.,0.5);
    controlHistos_.addHistogram("muchhadroniso", "; Isolation; Events", 100, 0.,0.5);
    controlHistos_.addHistogram("muneuhadroniso", "; Isolation; Events", 100, 0.,0.5);
    controlHistos_.addHistogram("mureliso", "; Isolation; Events", 100, 0.,1.);
    
    //dilepton control
    controlHistos_.addHistogram("dilepton_mass", ";Invariant Mass(l,l') [GeV/c^{2}]; Events", 100, 0.,300.);
    controlHistos_.addHistogram("dilepton_sumpt", ";#Sigma |#vec{p}_{T}| [GeV/c]; Events", 100, 0.,300.);
    controlHistos_.addHistogram("dilepton_pt", ";|#Sigma #vec{p}_{T}| [GeV/c]; Events", 100, 0.,300.);
	
    //jets
    controlHistos_.addHistogram("jetchhadenfrac",";f_{charged hadrons}; Jets",50,0,1);
    controlHistos_.addHistogram("jetneuthadenfrac",";f_{neutral hadrons}; Jets",50,0,1);
    controlHistos_.addHistogram("jetchemenfrac",";f_{charged electromagnetic}; Jets",50,0,1);
    controlHistos_.addHistogram("jetneutemenfrac",";f_{neutral electromagnetic}; Jets",50,0,1);
    controlHistos_.addHistogram("jetphoenfrac",";f_{photons}; Jets",50,0,1);
    controlHistos_.addHistogram("jetmuenfrac",";f_{muons}; Jets",50,0,1);
    controlHistos_.addHistogram("jetpt",";p_{T} [GeV/c]; Jets",100,0,200);
    controlHistos_.addHistogram("jeteta",";#eta; Jets",100,-2.5,2.5);
    controlHistos_.addHistogram("njets",";Jet multiplicity; Events",4,0,4);
    controlHistos_.addHistogram("njetsfinal",";Jet multiplicity; Events",4,0,4);
    controlHistos_.addHistogram("btags",";b tags (TCHE); Jets",100,-1,50);
    controlHistos_.addHistogram("bmult",";b tag multiplicity (TCHEL); Events",4,0,4);
    controlHistos_.addHistogram("bmultfinal",";b tag multiplicity (TCHEL); Events",4,0,4);
    for(int ibin=1; ibin<=controlHistos_.getHisto("njets","all")->GetXaxis()->GetNbins(); ibin++)
      {
	TString ilabel(""); ilabel+=(ibin-1);
	if(ibin==controlHistos_.getHisto("njets","all")->GetXaxis()->GetNbins()) ilabel="#geq"+ilabel;
	controlHistos_.getHisto("njets","all")->GetXaxis()->SetBinLabel(ibin,ilabel);
	controlHistos_.getHisto("njetsfinal","all")->GetXaxis()->SetBinLabel(ibin,ilabel);
	controlHistos_.getHisto("bmult","all")->GetXaxis()->SetBinLabel(ibin,ilabel);
	controlHistos_.getHisto("bmultfinal","all")->GetXaxis()->SetBinLabel(ibin,ilabel);
      }

    //MET
    controlHistos_.addHistogram("met", ";#slash{E}_{T} [GeV]; Events", 30,  0.,300.);
    controlHistos_.addHistogram("neutralhadetfrac", ";f_{neutral had} E_{T} [GeV]; Events", 100,  0.,1.);
    controlHistos_.addHistogram("neutralemetfrac", ";f_{neutral em} E_{T} [GeV]; Events", 100,  0.,1.);
    controlHistos_.addHistogram("chargedememetfrac", ";f_{charged em} E_{T} [GeV]; Events", 100,  0.,1.);
    controlHistos_.addHistogram("chargedhadetfrac", ";f_{charged had} E_{T} [GeV]; Events", 100,  0.,1.);
    controlHistos_.addHistogram("muonetfrac", ";f_{muons} E_{T} [GeV]; Events", 100,  0.,1.);

    //replicate histograms for each channel
    controlHistos_.initMonitorForStep("ee");
    controlHistos_.initMonitorForStep("emu");
    controlHistos_.initMonitorForStep("mumu");

  }catch(std::exception &e){
  }  
}

/// everything that needs to be done during the event loop
void TopDileptonEventAnalyzer::analyze(const edm::Event& event,const edm::EventSetup &iSetup)
{

  std::vector<std::string> selStreams;
  selStreams.push_back("all");
  selStreams.push_back("mumu");
  selStreams.push_back("emu");
  selStreams.push_back("ee");
  
  try{

    //
    //MC TRUTH, WEIGHT, ETC.
    //
    float weight(1),normWeight(1);
    int npuOOT(0),npuIT(0);
    int gentteventcode=gen::top::Event::UNKNOWN;
    float dyMass(0.);
    summaryHandler_.evSummary_.nmcparticles=0;
    if(!event.isRealData())
      {
	edm::Handle<float> puWeightHandle;
	event.getByLabel("puWeights","puWeight",puWeightHandle);
	if(puWeightHandle.isValid()) weight = *(puWeightHandle.product());

	edm::Handle<float> normPuWeightHandle;
	event.getByLabel("puWeights","normPuReweight", normPuWeightHandle );
	if(normPuWeightHandle.isValid()) normWeight = *(normPuWeightHandle.product());

	edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
	event.getByType(puInfoH);
	for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++)
	  {
	    if(it->getBunchCrossing()==0) npuIT += it->getPU_NumInteractions();
	    else npuOOT += it->getPU_NumInteractions();
	  }
	
	
	genEvent_.genLabel_ = objConfig_["Generator"].getParameter<edm::InputTag>("source");

	//filter on DY mass if required
	dyMass=genEvent_.getDYMass(event,iSetup);
	if(objConfig_["Generator"].getParameter<bool>("filterDYmassWindow") && dyMass>50) return;

	gentteventcode = genEvent_.assignTTEvent(event,iSetup);
	summaryHandler_.evSummary_.isSignal = ( gentteventcode == gen::top::Event::EE ||
						gentteventcode == gen::top::Event::EMU ||
						gentteventcode == gen::top::Event::MUMU ); 
	
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
	      {
		int ipart=summaryHandler_.evSummary_.nmcparticles;
		summaryHandler_.evSummary_.mcpx[ipart]= itt->get()->px();
		summaryHandler_.evSummary_.mcpy[ipart]=itt->get()->py();
		summaryHandler_.evSummary_.mcpz[ipart]=itt->get()->pz();
		summaryHandler_.evSummary_.mcen[ipart]=itt->get()->energy();
		summaryHandler_.evSummary_.mcid[ipart]=itt->get()->pdgId();
		summaryHandler_.evSummary_.nmcparticles++;
	    }
	  }
      }

    //
    //TRIGGER
    //
    summaryHandler_.evSummary_.hasTrigger=true;
    if(!event.isRealData())
      {
	edm::Handle<edm::TriggerResults> allTriggerBits_;
	event.getByLabel( objConfig_["Trigger"].getParameter<edm::InputTag>("source"), allTriggerBits_);
	const edm::TriggerNames &triggerNames = event.triggerNames( *allTriggerBits_);
	std::vector<std::string> triggerPaths = objConfig_["Trigger"].getParameter<std::vector<std::string> >("triggerPaths"); 
	summaryHandler_.evSummary_.hasTrigger=false;
	for (size_t itrig = 0; itrig != allTriggerBits_->size(); ++itrig)
	  {
	    std::string trigName = triggerNames.triggerName(itrig);
	    if( !allTriggerBits_->wasrun(itrig) ) continue;
	    if( allTriggerBits_->error(itrig) ) continue;
	    if( !allTriggerBits_->accept(itrig) ) continue;
	    if( find(triggerPaths.begin(), triggerPaths.end(), trigName) == triggerPaths.end() ) continue;
	    summaryHandler_.evSummary_.hasTrigger=true;
	    break;
	  }
      }

    //
    // VERTEX SELECTION
    //
    edm::Handle<reco::VertexCollection> hVtx;
    event.getByLabel(objConfig_["Vertices"].getParameter<edm::InputTag>("source"), hVtx);
    edm::Handle<reco::BeamSpot> beamSpot;
    event.getByLabel( objConfig_["Vertices"].getParameter<edm::InputTag>("beamSpot"), beamSpot);
    std::vector<reco::VertexRef> selVertices = getGoodVertices(hVtx,objConfig_["Vertices"]);
    edm::Handle< double > rhoH;
    event.getByLabel(objConfig_["Jets"].getParameter<edm::InputTag>("rho"),rhoH);
    float rho =*rhoH;
    for(size_t is=0; is<selStreams.size(); is++)
      {
	TString ctf=selStreams[is];
	controlHistos_.fillHisto("ngoodvertex", ctf,selVertices.size());
	for(size_t ivtx=0; ivtx<selVertices.size(); ivtx++)
	  {
	    TString vtype=(ivtx ==0 ? "" : "other");
	    controlHistos_.fillHisto(vtype+"vertex_ndof",ctf,selVertices[ivtx]->ndof() , weight);
	    try{
	      controlHistos_.fillHisto(vtype+"vertex_sumpt",ctf,getVertexMomentumFlux(selVertices[ivtx].get()) , weight);
	    }catch(std::exception &e){
	      //tracks might not have been saved
	    }
	  }
      }
    if(selVertices.size()==0) return;


    //
    // LEPTON SELECTION
    //
    edm::Handle<edm::View<reco::Candidate> > hMu;
    event.getByLabel(objConfig_["Muons"].getParameter<edm::InputTag>("source"), hMu);
    edm::Handle<edm::View<reco::Candidate> > hEle;
    event.getByLabel(objConfig_["Electrons"].getParameter<edm::InputTag>("source"), hEle);

    std::vector<CandidatePtr> selLooseMuons     = getGoodMuons(hMu, *beamSpot, rho, objConfig_["LooseMuons"]);
    std::vector<CandidatePtr> selLooseElectrons = getGoodElectrons(hEle, hMu, *beamSpot, rho, objConfig_["LooseElectrons"]);
    std::vector<CandidatePtr> selLooseLeptons   = selLooseMuons;
    selLooseLeptons.insert(selLooseLeptons.end(), selLooseElectrons.begin(), selLooseElectrons.end());

    std::vector<CandidatePtr> selMuons     = getGoodMuons(hMu, *beamSpot, rho, objConfig_["Muons"]);
    std::vector<CandidatePtr> selElectrons = getGoodElectrons(hEle, hMu, *beamSpot, rho, objConfig_["Electrons"]);
    std::vector<CandidatePtr> selLeptons   = selMuons;
    selLeptons.insert(selLeptons.end(), selElectrons.begin(), selElectrons.end());
    
    for(std::vector<CandidatePtr>::iterator lepIt = selLooseLeptons.begin(); lepIt != selLooseLeptons.end(); lepIt++)
      {
	reco::CandidatePtr &lep=*lepIt;
	float pt=lep->pt();
	if(pt<20) continue;
	std::vector<double> isol = getLeptonIso(lep,20.);
	TString lepType( fabs(getLeptonId(lep))==MUON ? "mu":"e" );
	for(size_t is=0; is<selStreams.size(); is++)
	  {
	    TString ctf=selStreams[is];
	    controlHistos_.fillHisto(lepType+"gammaiso", ctf,isol[G_ISO]/pt , weight);
	    controlHistos_.fillHisto(lepType+"chhadroniso",ctf,isol[C_ISO]/pt , weight);
	    controlHistos_.fillHisto(lepType+"neuhadroniso",ctf,isol[N_ISO]/pt , weight);
	    controlHistos_.fillHisto(lepType+"reliso",ctf,isol[PFREL_ISO] , weight);
	  }
      }    
    if(selLeptons.size()<2) return;

    selStreams.clear();
    selStreams.push_back("all");
    if(selElectrons.size()>1) selStreams.push_back("ee");
    if(selMuons.size()>1)      selStreams.push_back("mumu");
    if(selElectrons.size()>0 && selMuons.size()>0) selStreams.push_back("emu");
    for(size_t is=0; is<selStreams.size(); is++)
      {
	TString ctf=selStreams[is];
	controlHistos_.fillHisto("cutflow",ctf,1.,weight);
      }

    //
    // DILEPTON
    //
    std::vector<CandidatePtr> dilepton = getDileptonCandidate(selLeptons, objConfig_["Dileptons"], iSetup);
    int selPath =  getDileptonId(dilepton);
    if(selPath==UNKNOWN) return;
    reco::CandidatePtr lepton1 = dilepton[0];
    reco::CandidatePtr lepton2 = dilepton[1];
    LorentzVector lepton1P = lepton1->p4();
    LorentzVector lepton2P = lepton2->p4();
    LorentzVector dileptonP=lepton1P+lepton2P;
    float dilsumpt=lepton1P.pt()+lepton2P.pt();
    float dilpt=dileptonP.pt();
    float dilmass=dileptonP.mass();
    bool isZCand( ( (selPath==EE || selPath==MUMU) && fabs(dilmass-91)<15) );
    bool isOS(lepton1->charge()*lepton2->charge()<0);
    double minDileptonMass = objConfig_["Dileptons"].getParameter<double>("minDileptonMass");
    double maxDileptonMass = objConfig_["Dileptons"].getParameter<double>("maxDileptonMass");
    if(dilmass<minDileptonMass || dilmass>maxDileptonMass) return;

    selStreams.clear();
    selStreams.push_back("all");
    if(selPath==EE)   selStreams.push_back("ee");
    if(selPath==EMU)  selStreams.push_back("emu");
    if(selPath==MUMU) selStreams.push_back("mumu");
    for(size_t is=0; is<selStreams.size(); is++)
      {
	TString ctf=selStreams[is];
	controlHistos_.fillHisto("dilepton_sumpt",ctf,dilsumpt,weight);
	controlHistos_.fillHisto("dilepton_pt",ctf,dilpt,weight);
	controlHistos_.fillHisto("dilepton_mass",ctf,dilmass,weight);

	if(!isZCand) controlHistos_.fillHisto("cutflow",ctf,2.,weight);
      }

    //leptons to save in the tree (dilepton+all selected loose leptons) 
    std::vector<reco::CandidatePtr> leptons;
    leptons.push_back(lepton1);
    leptons.push_back(lepton2);
    for(vector<CandidatePtr>::iterator lit = selLooseLeptons.begin(); lit != selLooseLeptons.end(); lit++)
      {
	if(deltaR((*lit)->eta(),(*lit)->phi(),lepton1->eta(),lepton1->phi())<0.1 
	   || deltaR((*lit)->eta(),(*lit)->phi(),lepton2->eta(),lepton2->phi())<0.1 ) continue;
	leptons.push_back(*lit);
      }
    
    //
    // JETS
    //
    //add also the jets                                                                                                                                                             
    edm::Handle<edm::View<reco::Candidate> > hJet;
    event.getByLabel(objConfig_["Jets"].getParameter<edm::InputTag>("source"), hJet);
    std::vector<CandidatePtr> jets = getGoodJets(hJet, selLeptons, objConfig_["Jets"]);
    int njets(0), nbjets(0);
    std::vector<const pat::Jet *> selJets;
    for(std::vector<CandidatePtr>::iterator jit = jets.begin(); jit!=jets.end(); jit++)
      {
	const pat::Jet *j=dynamic_cast<const pat::Jet*>(jit->get());
	selJets.push_back(j);
	
	if(j->pt()>30 && fabs(j->eta())<2.5) njets++;
	float btag=j->bDiscriminator("trackCountingHighEffBJetTags");
	if(btag>1.74) nbjets+=1; //loose point
	
	//monitor jets
	if(!isZCand)
	  {       
	    for(size_t is=0; is<selStreams.size(); is++)
	      {
		TString ctf=selStreams[is];
		{
		  controlHistos_.fillHisto("jetpt",ctf,j->pt(),weight);
		  controlHistos_.fillHisto("jeteta",ctf,j->eta(),weight);
		  controlHistos_.fillHisto("btags",ctf,btag,weight);
		  
		  controlHistos_.fillHisto("jetchhadenfrac",ctf, j->chargedHadronEnergyFraction(), weight );
		  controlHistos_.fillHisto("jetneuthadenfrac",ctf, j->neutralHadronEnergyFraction(), weight );
		  controlHistos_.fillHisto("jetchemenfrac",ctf, j->chargedEmEnergyFraction(),weight );
		  controlHistos_.fillHisto("jetneutemenfrac",ctf, j->neutralEmEnergyFraction(),weight );
		  controlHistos_.fillHisto("jetphoenfrac",ctf, j->photonEnergyFraction(),weight );
		  controlHistos_.fillHisto("jetmuenfrac",ctf, j->muonEnergyFraction(),weight );
		}
	      }
	  }
      }

    //monitor jet multiplicity
    if(!isZCand)
      {
	for(size_t is=0; is<selStreams.size(); is++)
	  {
	    TString ctf=selStreams[is];
	    controlHistos_.fillHisto("njets",ctf,njets,weight);
	  }
      }
    bool passJets(njets>=2);
    int btagbin(nbjets);
    if(btagbin>2) btagbin=2;
    
    //monitor b-tag multiplicity for events with >= 2 jets
    if(!isZCand && passJets)
      {
	for(size_t is=0; is<selStreams.size(); is++)
	  {
	    TString ctf=selStreams[is];
	    controlHistos_.fillHisto("cutflow",ctf,3.,weight);
	    controlHistos_.fillHisto("bmult",ctf,nbjets,weight);
	  }
      }

    //
    // MET
    //
    edm::Handle<edm::View<pat::MET> > hMET;
    event.getByLabel(objConfig_["MET"].getParameter<edm::InputTag>("source"), hMET);
    const pat::MET &pfmet = (*hMET)[0];
    LorentzVector met(pfmet.px(),pfmet.py(),0,pfmet.pt());
    float neutralhadetfrac=pfmet.NeutralHadEtFraction();
    float neutralemetfrac=pfmet.NeutralEMFraction();
    float chargedemetfrac=pfmet.ChargedEMEtFraction();
    float chargedhadetfrac=pfmet.ChargedHadEtFraction();
    float muonetfrac=pfmet.MuonEtFraction();
    bool passMET( selPath==EMU || met.pt()>30);  
    if(!isZCand)
      {
	for(size_t is=0; is<selStreams.size(); is++)
	  {
	    TString ctf=selStreams[is];
	    if(passJets)
	      {
		controlHistos_.fillHisto("met",ctf,met.pt(),weight);
		controlHistos_.fillHisto("neutralhadetfrac",ctf,neutralhadetfrac,weight);
		controlHistos_.fillHisto("neutralemetfrac",ctf,neutralemetfrac,weight);
		controlHistos_.fillHisto("chargedemetfrac",ctf,chargedemetfrac,weight);
		controlHistos_.fillHisto("chargedhadetfrac",ctf,chargedhadetfrac,weight);
		controlHistos_.fillHisto("muonetfrac",ctf,muonetfrac,weight);
	      }
	    
	    if(!passMET) continue;
	    if(passJets) controlHistos_.fillHisto("cutflow",ctf,4,weight);
	    if(!isOS) continue;
	
	    if(passJets) controlHistos_.fillHisto("cutflow",ctf,5,weight);
	    controlHistos_.fillHisto("njetsfinal",ctf,njets,weight);
	    if(passJets) 
	      {
		controlHistos_.fillHisto("bmultfinal",ctf,nbjets,weight);
		controlHistos_.fillHisto("cutflow",ctf,6.+btagbin,weight);
	      }
	  }
      }
    
    //
    // ALL DONE: SAVE TO NTUPLE
    //    
    saveEvent(event,selPath,leptons,selJets,&pfmet,selVertices.size(),npuIT,rho,weight,normWeight);
    
  } catch(std::exception &e) {
    std::cout << "[TopDileptonEventAnalyzer][analyze] failed with " << e.what() << std::endl;
  }
  
}

//
void TopDileptonEventAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup)
{
  TString streams[]={"all","ee","mumu","emu"};
  edm::Handle<edm::MergeableCounter> ctrHandle;
  iLumi.getByLabel("startCounter", ctrHandle);
  for(size_t istream=0; istream<sizeof(streams)/sizeof(TString); istream++)
    controlHistos_.fillHisto("cutflow",streams[istream],0.,ctrHandle->value);
}

//
void TopDileptonEventAnalyzer::saveEvent(const edm::Event& event, int evCat, std::vector<reco::CandidatePtr> &leptons, std::vector<const pat::Jet *> &jets, const pat::MET *met, 
				     int nvertices, int npuIT, float rho, float weight, float normWeight)
{
  //save event header
  summaryHandler_.evSummary_.run=event.id().run();
  summaryHandler_.evSummary_.lumi=event.luminosityBlock();
  summaryHandler_.evSummary_.event=event.id().event();
  summaryHandler_.evSummary_.cat=evCat; 
  summaryHandler_.evSummary_.weight=weight; 
  summaryHandler_.evSummary_.normWeight=normWeight; 
  summaryHandler_.evSummary_.nvtx=nvertices;
  summaryHandler_.evSummary_.ngenpu=npuIT;
  summaryHandler_.evSummary_.rho=rho;
  summaryHandler_.evSummary_.nparticles=leptons.size()+jets.size()+1;

  //save the leptons
  for(size_t ilepton=0; ilepton<leptons.size(); ilepton++)
    {
      summaryHandler_.evSummary_.px[ilepton]=leptons[ilepton].get()->px();
      summaryHandler_.evSummary_.py[ilepton]=leptons[ilepton].get()->py();
      summaryHandler_.evSummary_.pz[ilepton]=leptons[ilepton].get()->pz();
      summaryHandler_.evSummary_.en[ilepton]=leptons[ilepton].get()->energy();

      std::vector<double> liso=getLeptonIso(leptons[ilepton],0.);
      summaryHandler_.evSummary_.info1[ilepton]=getLeptonPtError(leptons[ilepton]);
      summaryHandler_.evSummary_.info2[ilepton]=liso[C_ISO];
      summaryHandler_.evSummary_.info3[ilepton]=liso[G_ISO];
      summaryHandler_.evSummary_.info4[ilepton]=liso[N_ISO];

      int id = getLeptonId(leptons[ilepton]);
      const reco::GenParticle *gen=getLeptonGenMatch(leptons[ilepton]);
      int genid( gen ? gen->pdgId() : -9999);  
      summaryHandler_.evSummary_.id[ilepton] = id;
      summaryHandler_.evSummary_.genid[ilepton] = genid;
      summaryHandler_.evSummary_.genflav[ilepton] = genid;
      
      if(fabs(id)==11)
	{
	  const pat::Electron *ele = dynamic_cast<const pat::Electron *>( leptons[ilepton].get() );
	  summaryHandler_.evSummary_.info5[ilepton]= ( (int(ele->electronID("eidVBTF70")) & 0x1) )
            | ( (int(ele->electronID("eidVBTF80")) & 0x1) << 1)
            | ( (int(ele->electronID("eidVBTF85")) & 0x1) << 2)	    
            | ( (int(ele->electronID("eidVBTF90")) & 0x1) << 3)
            | ( (int(ele->electronID("eidVBTF95")) & 0x1) << 4);
	}
      else
	{
	  const pat::Muon *mu = dynamic_cast<const pat::Muon *>( leptons[ilepton].get() );
	  summaryHandler_.evSummary_.info5[ilepton]=  ( (int(mu->muonID("GlobalMuonPromptTight")) & 0x1) )
            | ( (int(mu->muonID("TMLastStationLoose")) & 0x1) << 1)
            | ( (int(mu->muonID("TMLastStationTight")) & 0x1) << 2)
            | ( (int(mu->muonID("TMLastStationAngTight")) & 0x1) << 3);
	}
    }

  //save the jets
  for(size_t ijet=0; ijet<jets.size(); ijet++)
    {
      int pidx = leptons.size()+ijet;
      summaryHandler_.evSummary_.px[pidx]=jets[ijet]->px();
      summaryHandler_.evSummary_.py[pidx]=jets[ijet]->py();
      summaryHandler_.evSummary_.pz[pidx]=jets[ijet]->pz();
      summaryHandler_.evSummary_.en[pidx]=jets[ijet]->energy();
      summaryHandler_.evSummary_.id[pidx] = 1;
      const reco::Candidate *genParton = jets[ijet]->genParton();
      summaryHandler_.evSummary_.genid[pidx] = genParton ? genParton->pdgId() : -9999;
      summaryHandler_.evSummary_.genflav[pidx] = jets[ijet]->partonFlavour();
      summaryHandler_.evSummary_.info1[pidx]=jets[ijet]->bDiscriminator("trackCountingHighEffBJetTags");
      summaryHandler_.evSummary_.info2[pidx]=jets[ijet]->bDiscriminator("trackCountingHighPurBJetTags");
      summaryHandler_.evSummary_.info3[pidx]=jets[ijet]->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");	  
      summaryHandler_.evSummary_.info4[pidx]=jets[ijet]->bDiscriminator("jetBProbabilityBJetTags");
      summaryHandler_.evSummary_.info5[pidx]=jets[ijet]->bDiscriminator("jetProbabilityBJetTags");
      summaryHandler_.evSummary_.info6[pidx]=jets[ijet]->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
      summaryHandler_.evSummary_.info7[pidx]=jets[ijet]->bDiscriminator("combinedSecondaryVertexBJetTags");
      summaryHandler_.evSummary_.info8[pidx]=jets[ijet]->neutralHadronEnergyFraction();
      summaryHandler_.evSummary_.info9[pidx]=jets[ijet]->chargedEmEnergyFraction()+jets[ijet]->neutralEmEnergyFraction();
    }
  
  //save met
  int pidx=leptons.size()+jets.size();
  summaryHandler_.evSummary_.px[pidx]=met->px();
  summaryHandler_.evSummary_.py[pidx]=met->py();
  summaryHandler_.evSummary_.pz[pidx]=met->pz();
  summaryHandler_.evSummary_.en[pidx]=met->energy();
  summaryHandler_.evSummary_.id[pidx] = 0;
  summaryHandler_.evSummary_.genid[pidx] = -9999;
  summaryHandler_.evSummary_.genflav[pidx] = -9999;
  
  //no further measurements for now
  summaryHandler_.evSummary_.nmeasurements=0;

  //all done
  summaryHandler_.fillTree();
}



DEFINE_FWK_MODULE(TopDileptonEventAnalyzer);


