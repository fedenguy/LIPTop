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

#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"

#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class DileptonEventCleaner : public edm::EDAnalyzer 
{
public:
  explicit DileptonEventCleaner(const edm::ParameterSet& cfg);
  ~DileptonEventCleaner(){};
  virtual void analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup) ;
  void endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup);

private:

  void saveEvent(const edm::Event& event, int evCat, std::vector<reco::CandidatePtr> &leptons, std::vector<const pat::Jet *> &jets, const pat::MET *met, 
		 int nvertices, int npuIT, float rho, float weight);

  std::map<std::string, edm::ParameterSet> objConfig_;
  
  EventSummaryHandler summaryHandler_;
  SelectionMonitor controlHistos_;

};

using namespace std;

/// default constructor
DileptonEventCleaner::DileptonEventCleaner(const edm::ParameterSet& cfg)
{
  try{

    edm::Service<TFileService> fs;

    summaryHandler_.initTree( fs->make<TTree>("data","Event Summary") );

    objConfig_["Vertices"] = cfg.getParameter<edm::ParameterSet>("Vertices");
    objConfig_["Jets"] = cfg.getParameter<edm::ParameterSet>("Jets");

    TFileDirectory baseDir=fs->mkdir(cfg.getParameter<std::string>("dtag"));    
    TString streams[]={"ee","mumu","emu"};
    TString selSteps[]={"Reco","2 leptons","Z-veto","#geq 2 jets","MET>40,0","=0 b-tags","=1 b-tags", "#geq 2 b-tags"};
    const size_t nselsteps=sizeof(selSteps)/sizeof(TString);
    controlHistos_.addHistogram( baseDir.make<TH1F>("cutflow", ";Step; Events",nselsteps,0,nselsteps),false ); 
    for(size_t istream=0; istream<sizeof(streams)/sizeof(TString); istream++)
      {
	TString cat=streams[istream];
	TFileDirectory newDir=baseDir.mkdir(cat.Data());

	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_cutflow", ";Step; Events",nselsteps,0,nselsteps),false ); 
	for(int ibin=1; ibin<=controlHistos_.getHisto(cat+"_cutflow","all")->GetXaxis()->GetNbins(); ibin++)
	  {
	    controlHistos_.getHisto(cat+"_cutflow","all")->GetXaxis()->SetBinLabel(ibin,selSteps[ibin-1]);
	    controlHistos_.getHisto("cutflow","all")->GetXaxis()->SetBinLabel(ibin,selSteps[ibin-1]);
	  }
      
	//dilepton control
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_dilepton_mass", ";Invariant Mass(l,l') [GeV/c^{2}]; Events", 100, 0.,300.),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_dilepton_mass", ";Invariant Mass(l,l') [GeV/c^{2}]; Events", 100, 0.,300.),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_dilepton_sumpt", ";#Sigma |#vec{p}_{T}| [GeV/c]; Events", 100, 0.,300.),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_dilepton_pt", ";|#Sigma #vec{p}_{T}| [GeV/c]; Events", 100, 0.,300.),false ); 

	//vertex control
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_vertex_sumpt", ";#Sigma_{tracks} p_{T} [GeV/c]; Events", 100, 0.,300.),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_othervertex_sumpt", ";#Sigma_{tracks} p_{T} [GeV/c]; Events", 100, 0.,300.),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_vertex_pt", ";|#Sigma_{tracks} #vec{p}_{T}| [GeV/c]; Events", 100, 0.,300.),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_othervertex_pt", ";|#Sigma_{tracks} #vec{p}_{T}| [GeV/c]; Events", 100, 0.,300.),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_ngoodvertex", ";Vertices; Events", 25, 0.,25.) ); 
      
	//jets
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_jetpt",";p_{T} [GeV/c]; Jets",100,0,200),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_pujetpt",";p_{T} [GeV/c]; Jets",100,0,200),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_jeteta",";#eta; Jets",100,-2.5,2.5),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_jetfassoc",";f_{assoc}; Jets",100,0,1),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_pujetfassoc",";f_{assoc}; Jets",100,0,1),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_njets",";Jet multiplicity; Events",4,0,4),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_npujets",";Jet multiplicity; Events",4,0,4),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_bmult",";b tag multiplicity (TCHEL); Events",4,0,4),false ); 
	for(int ibin=1; ibin<=controlHistos_.getHisto(cat+"_njets","all")->GetXaxis()->GetNbins(); ibin++)
	  {
	    TString ilabel(""); ilabel+=(ibin-1);
	    if(ibin==controlHistos_.getHisto(cat+"_npujets","all")->GetXaxis()->GetNbins()) ilabel="#geq"+ilabel;
	    controlHistos_.getHisto(cat+"_npujets","all")->GetXaxis()->SetBinLabel(ibin,ilabel);
	    controlHistos_.getHisto(cat+"_bmult","all")->GetXaxis()->SetBinLabel(ibin,ilabel);
	  }
	
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_btags",";b tags (TCHE); Jets",100,-1,50),false ); 
	controlHistos_.addHistogram( newDir.make<TH1F>(cat+"_met", ";#slash{E}_{T} [GeV]; Events", 30,  0.,300.),false ); 

	//meant inclusive validation of JET/MET
	if(istream==0)
	  {
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetpt",";p_{T} [GeV/c]; Jets",100,0,200),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jeteta",";#eta; Jets",100,-2.5,2.5),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetchhadenfrac",";f_{charged hadrons}; Jets",50,0,1),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetneuthadenfrac",";f_{neutral hadrons}; Jets",50,0,1),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetchemenfrac",";f_{charged electromagnetic}; Jets",50,0,1),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetneutemenfrac",";f_{neutral electromagnetic}; Jets",50,0,1),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetphoenfrac",";f_{photons}; Jets",50,0,1),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetmuenfrac",";f_{muons}; Jets",50,0,1),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetchhaden",";Charged hadron components in jets [GeV/c]; Jets",50,0,200),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetneuthaden",";Neutral hadrons component in jets [GeV/c]; Jets",50,0,200),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetchemen",";Charged electromagnetic component in jets [GeV/c]; Jets",50,0,200),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetneutemen",";Neutral electromagnetic component in jets [GeV/c]; Jets",50,0,200),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetphoen",";Photon component in jets [GeV/c]; Jets",50,0,200),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("jetmuen",";Muon component in jets [GeV/c]; Jets",50,0,200),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("met", ";#slash{E}_{T} [GeV]; Events", 30,  0.,300.),false ); 
	    //controlHistos_.addHistogram( baseDir.make<TH1F>("sumet", ";#sum E_{T} [GeV]; Events", 100,  0.,1000.),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("neutralhadetfrac", ";f_{neutral had} E_{T} [GeV]; Events", 100,  0.,1.),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("neutralemetfrac", ";f_{neutral em} E_{T} [GeV]; Events", 100,  0.,1.),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("chargedememetfrac", ";f_{charged em} E_{T} [GeV]; Events", 100,  0.,1.),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("chargedhadetfrac", ";f_{charged had} E_{T} [GeV]; Events", 100,  0.,1.),false ); 
	    controlHistos_.addHistogram( baseDir.make<TH1F>("muonetfrac", ";f_{muons} E_{T} [GeV]; Events", 100,  0.,1.),false ); 
	  }
      }
  }catch(std::exception &e){
  }  
}

/// everything that needs to be done during the event loop
void DileptonEventCleaner::analyze(const edm::Event& event,const edm::EventSetup &iSetup)
{
  
  try{

    //get the weight for the event
    float weight=1;
    if(!event.isRealData())
      {
	edm::Handle<float> puWeightHandle;
	event.getByLabel("puWeights","puWeight",puWeightHandle);
	if(puWeightHandle.isValid()) weight = *(puWeightHandle.product());
      }

    //get objects for this event
    edm::Handle<std::vector<pat::EventHypothesis> > evHandle;
    event.getByLabel(edm::InputTag("cleanEvent:selectedEvent"),evHandle);
    if(!evHandle.isValid()) 
      {
	cout << "No valid event hypothesis" << endl;
	return;
      }
    const pat::EventHypothesis &evhyp = (*(evHandle.product()))[0];
    
    edm::Handle<std::vector<reco::Vertex> > vertexHandle;
    event.getByLabel(edm::InputTag("cleanEvent:selectedVertices"),vertexHandle);
    if(!vertexHandle.isValid()) 
      {
	cout << "No vertex selected" << endl;
	return;
      }

    //get other vertices
    std::vector<reco::VertexRef> selVertices;
    try{
      edm::Handle<reco::VertexCollection> hVtx;
      event.getByLabel(objConfig_["Vertices"].getParameter<edm::InputTag>("source"), hVtx);
      selVertices = vertex::filter(hVtx,objConfig_["Vertices"]);
    }catch(std::exception &e){
    }
        
    edm::Handle< std::vector<int> > selInfo;
    event.getByLabel(edm::InputTag("cleanEvent:selectionInfo"),selInfo);
    if(!selInfo.isValid()) 
      {
	cout << "No valid selection info" << endl;
	return;
      }

    edm::Handle< double > rhoH;
    event.getByLabel(objConfig_["Jets"].getParameter<edm::InputTag>("rho"),rhoH);
    float rho =*rhoH;
   
    int selPath = (*(selInfo.product()))[0];
    int selStep = (*(selInfo.product()))[1];
    int tteventcode = (*(selInfo.product()))[2]; 

    //require that a dilepton has benn selected
    if(selPath==0 or selStep<3) return;
    std::string istream="mumu";
    if(selPath==2) istream="ee";
    if(selPath==3) istream="emu";
    controlHistos_.fillHisto(istream+"_cutflow","all",1,weight);
    controlHistos_.fillHisto("cutflow","all",1,weight);
    
    //vertex quantities
    controlHistos_.fillHisto(istream+"_ngoodvertex","all",selVertices.size(),weight);
    const reco::Vertex &primVertex = (*(vertexHandle.product()))[0];
    //MC truth on pileup (if available)
    int npuOOT(0),npuIT(0);
    if(!event.isRealData())
      {
	edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
	event.getByType(puInfoH);
	if(puInfoH.isValid())
	  {
	    for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++)
	      {
		if(it->getBunchCrossing()==0) npuIT += it->getPU_NumInteractions();
		else npuOOT += it->getPU_NumInteractions();
	      }
	  }
      }

    //basic dilepton kinematics
    reco::CandidatePtr lepton1 = evhyp["leg1"];
    LorentzVector lepton1P = lepton1->p4();
    reco::CandidatePtr lepton2 = evhyp["leg2"];
    LorentzVector lepton2P = lepton2->p4();
    LorentzVector dileptonP=lepton1P+lepton2P;
    controlHistos_.fillHisto(istream+"_dilepton_sumpt","all",lepton1P.pt()+lepton2P.pt(),weight);
    controlHistos_.fillHisto(istream+"_dilepton_pt","all",dileptonP.pt(),weight);
    controlHistos_.fillHisto(istream+"_dilepton_mass","all",dileptonP.mass(),weight);

    //Z+quarkonia veto
    bool isZCand(false);
    if(dileptonP.mass()<20) return;
    if( (istream=="ee" || istream=="mumu") && fabs(dileptonP.mass()-91)<15) isZCand=true;
    if(!isZCand)
      {
	controlHistos_.fillHisto(istream+"_cutflow","all",2,weight);
	controlHistos_.fillHisto("cutflow","all",2,weight);
      }
    
    //count the jets in the event
    std::vector<reco::CandidatePtr> seljets= evhyp.all("jet");
    int njets(0), nbjets(0);
    std::vector<const pat::Jet *> selJets,assocJets;
    for (pat::eventhypothesis::Looper<pat::Jet> jet = evhyp.loopAs<pat::Jet>("jet"); jet; ++jet) 
      {
	
	assocJets.push_back( jet.get() );
	if(jet->pt()<30 || fabs(jet->eta())>2.5) continue;
	selJets.push_back( jet.get() );
	
	float btag=jet->bDiscriminator("trackCountingHighEffBJetTags");
	if(btag>1.74) nbjets+=1; //loose point

	if(!isZCand)
	  {
	    controlHistos_.fillHisto(istream+"_jetpt","all",jet->pt(),weight);
	    controlHistos_.fillHisto(istream+"_jeteta","all",jet->eta(),weight);
	    controlHistos_.fillHisto(istream+"_jetfassoc","all", jet::fAssoc( jet.get(), &primVertex), weight );
	    controlHistos_.fillHisto(istream+"_btags","all",btag,weight);

	    //simulation validation
	    controlHistos_.fillHisto("jetpt","all",jet->pt(),weight);
	    controlHistos_.fillHisto("jeteta","all",jet->eta(),weight);
	    controlHistos_.fillHisto("jetchhadenfrac","all", jet->chargedHadronEnergyFraction(), weight );
	    controlHistos_.fillHisto("jetneuthadenfrac","all", jet->neutralHadronEnergyFraction(), weight );
	    controlHistos_.fillHisto("jetchemenfrac","all", jet->chargedEmEnergyFraction(),weight );
	    controlHistos_.fillHisto("jetneutemenfrac","all", jet->neutralEmEnergyFraction(),weight );
	    controlHistos_.fillHisto("jetphoenfrac","all", jet->photonEnergyFraction(),weight );
	    controlHistos_.fillHisto("jetmuenfrac","all", jet->muonEnergyFraction(),weight );
	    controlHistos_.fillHisto("jetchhaden","all", jet->chargedHadronEnergy(), weight );
	    controlHistos_.fillHisto("jetneuthaden","all", jet->neutralHadronEnergy(), weight );
	    controlHistos_.fillHisto("jetchemen","all", jet->chargedEmEnergy(),weight );
	    controlHistos_.fillHisto("jetneutemen","all", jet->neutralEmEnergy(),weight );
	    controlHistos_.fillHisto("jetphoen","all", jet->photonEnergy(),weight );
	    controlHistos_.fillHisto("jetmuen","all", jet->muonEnergy(),weight );
	  }
      }
    njets=selJets.size();
    if(!isZCand)
      {
	controlHistos_.fillHisto(istream+"_njets","all",njets,weight);
	controlHistos_.fillHisto(istream+"_bmult","all",nbjets,weight);
      }

    //count the pu jets
    std::vector<reco::CandidatePtr> pujets= evhyp.all("pujet");
    int npujets(0);
    std::vector<const pat::Jet *> puJets;
    for (pat::eventhypothesis::Looper<pat::Jet> jet = evhyp.loopAs<pat::Jet>("pujet"); jet; ++jet) 
      {
	npujets++;
	puJets.push_back(jet.get());
	if(!isZCand)
	  {
	    controlHistos_.fillHisto(istream+"_pujetpt","all",jet->pt(),weight);
	    controlHistos_.fillHisto(istream+"_pujetfassoc","all",jet::fAssoc(jet.get(),&primVertex),weight);
	  }
      }
    if(!isZCand) controlHistos_.fillHisto(istream+"_npujets","all",npujets,weight);


    //require two jets
    if(njets<2) return;
    if(!isZCand)
      {
	controlHistos_.fillHisto(istream+"_cutflow","all",3,weight);
	controlHistos_.fillHisto("cutflow","all",3,weight);
      }

    //base met kinematics
    const pat::MET *evmet = dynamic_cast<const pat::MET *>(evhyp["met"].get());
    LorentzVector met(evmet->px(),evmet->py(),0,evmet->pt());
    //float sumet = evmet->sumEt();
    float neutralhadetfrac=evmet->NeutralHadEtFraction();
    float neutralemetfrac=evmet->NeutralEMFraction();
    float chargedemetfrac=evmet->ChargedEMEtFraction();
    float chargedhadetfrac=evmet->ChargedHadEtFraction();
    float muonetfrac=evmet->MuonEtFraction();

    //control  
    if(!isZCand)
      {
	controlHistos_.fillHisto("met","all",met.pt(),weight);
	controlHistos_.fillHisto(istream+"_met","all",met.pt(),weight);
	//controlHistos_.fillHisto("sumet","all",sumet,weight);
	controlHistos_.fillHisto("neutralhadetfrac","all",neutralhadetfrac,weight);
	controlHistos_.fillHisto("neutralemetfrac","all",neutralemetfrac,weight);
	controlHistos_.fillHisto("chargedemetfrac","all",chargedemetfrac,weight);
	controlHistos_.fillHisto("chargedhadetfrac","all",chargedhadetfrac,weight);
	controlHistos_.fillHisto("muonetfrac","all",muonetfrac,weight);
      }

    //require met for same flavor channels
    //if(metP.pt()<20) return;
    if( (istream=="ee" || istream=="mumu") && met.pt()<40) return;
    if(!isZCand)
      {
	controlHistos_.fillHisto(istream+"_cutflow","all",4,weight);
	controlHistos_.fillHisto("cutflow","all",4,weight);
      }


    //b-tagged sample
    int btagbin= nbjets;
    if(btagbin>2) btagbin=2;
    if(!isZCand)
      {
	controlHistos_.fillHisto(istream+"_cutflow","all",5+btagbin,weight);
	controlHistos_.fillHisto("cutflow","all",5+btagbin,weight);
      }

    //MC truth on the event
    if(!event.isRealData())
      {
	std::string mcparticles[]={"top","quarks","leptons","neutrinos"};
        summaryHandler_.evSummary_.nmcparticles=0;
	for(size_t imcpart=0; imcpart<sizeof(mcparticles)/sizeof(std::string); imcpart++)
	  {
	    for (pat::eventhypothesis::Looper<reco::GenParticle> genpart = evhyp.loopAs<reco::GenParticle>(mcparticles[imcpart]); genpart; ++genpart)
	      {
		summaryHandler_.evSummary_.mcpx[summaryHandler_.evSummary_.nmcparticles] = genpart->px();
		summaryHandler_.evSummary_.mcpy[summaryHandler_.evSummary_.nmcparticles] = genpart->py();
		summaryHandler_.evSummary_.mcpz[summaryHandler_.evSummary_.nmcparticles] = genpart->pz();
		summaryHandler_.evSummary_.mcen[summaryHandler_.evSummary_.nmcparticles]=genpart->energy();
		summaryHandler_.evSummary_.mcid[summaryHandler_.evSummary_.nmcparticles]=genpart->pdgId();
		summaryHandler_.evSummary_.nmcparticles++;
	      }
	  }
      }

    std::vector<reco::CandidatePtr> leptons;
    leptons.push_back(lepton1);
    leptons.push_back(lepton2);
    saveEvent(event,selPath,leptons,selJets,evmet,selVertices.size(),npuIT,rho,weight);
    
  }catch(std::exception &e){
    std::cout << "[DileptonEventCleaner][analyze] failed with " << e.what() << std::endl;
  }

}

//
void DileptonEventCleaner::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup)
{
  cout << "[DileptonEventCleaner][endLuminosityBlock]" << endl;
  TString streams[]={"ee","mumu","emu"};
  edm::Handle<edm::MergeableCounter> ctrHandle;
  iLumi.getByLabel("startCounter", ctrHandle);
  for(size_t istream=0; istream<sizeof(streams)/sizeof(TString); istream++)
    {
      controlHistos_.fillHisto(streams[istream]+"_cutflow","all",0.,ctrHandle->value);
      if(istream==0) controlHistos_.fillHisto("cutflow","all",0.,ctrHandle->value);
    }
}

//
void DileptonEventCleaner::saveEvent(const edm::Event& event, int evCat, std::vector<reco::CandidatePtr> &leptons, std::vector<const pat::Jet *> &jets, const pat::MET *met, 
				     int nvertices, int npuIT, float rho, float weight)
{
  //save event header
  summaryHandler_.evSummary_.run=event.id().run();
  summaryHandler_.evSummary_.lumi=event.luminosityBlock();
  summaryHandler_.evSummary_.event=event.id().event();
  summaryHandler_.evSummary_.cat=evCat; 
  summaryHandler_.evSummary_.weight=weight; 
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

      int id(11);
      const reco::GenParticle *gen=0;
      const pat::Electron *ele = dynamic_cast<const pat::Electron *>(leptons[ilepton].get());
      if(ele)
	{
	  id *= ele->charge();
	  gen=ele->genLepton();
	}
      else
	{
	  const pat::Muon *mu = dynamic_cast<const pat::Muon *>(leptons[ilepton].get());
	  if(mu)
	    {
	      id=13;
	      id *= mu->charge();
	      gen=mu->genLepton();	 
	    }
	}
      int genid( gen? gen->pdgId() : -9999);  

      summaryHandler_.evSummary_.id[ilepton] = id;
      summaryHandler_.evSummary_.genid[ilepton] = genid;
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
      summaryHandler_.evSummary_.info1[pidx]=jets[ijet]->bDiscriminator("trackCountingHighEffBJetTags");
      summaryHandler_.evSummary_.info2[pidx]=jets[ijet]->bDiscriminator("trackCountingHighPurBJetTags");
      summaryHandler_.evSummary_.info3[pidx]=jets[ijet]->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");	  
      summaryHandler_.evSummary_.info4[pidx]=jets[ijet]->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
      summaryHandler_.evSummary_.info5[pidx]=jets[ijet]->bDiscriminator("combinedSecondaryVertexBJetTags");
    }

  //save met
  int pidx=leptons.size()+jets.size();
  summaryHandler_.evSummary_.px[pidx]=met->px();
  summaryHandler_.evSummary_.py[pidx]=met->py();
  summaryHandler_.evSummary_.pz[pidx]=met->pz();
  summaryHandler_.evSummary_.en[pidx]=met->energy();
  summaryHandler_.evSummary_.id[pidx] = 0;
  summaryHandler_.evSummary_.genid[pidx] = -9999;

  //all done
  summaryHandler_.fillTree();
}



DEFINE_FWK_MODULE(DileptonEventCleaner);


