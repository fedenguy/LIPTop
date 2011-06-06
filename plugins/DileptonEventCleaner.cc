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

  void saveEvent(const edm::Event& event, int evCat, std::vector<reco::CandidatePtr> &leptons, std::vector<const pat::Jet *> &jets, const reco::PFMET *met, 
		 int nvertices, int npuIT, float rho, float weight);

  inline TH1 *getHist(TString key)
  {
    if(results_.find(key)==results_.end()) return 0;
    return (TH1 *)results_[key];
  }
  std::map<TString, TObject *>  results_;
  std::map<std::string, edm::ParameterSet> objConfig_;
  
  EventSummaryHandler summaryHandler_;
};

using namespace std;

/// default constructor
DileptonEventCleaner::DileptonEventCleaner(const edm::ParameterSet& cfg)
{
  try{

    edm::Service<TFileService> fs;

    summaryHandler_.initTree( fs->make<TTree>("data","Event Summary") );

    objConfig_["Vertices"] = cfg.getParameter<edm::ParameterSet>("Vertices");
       
    TFileDirectory baseDir=fs->mkdir(cfg.getParameter<std::string>("dtag"));    
    TString streams[]={"ee","mumu","emu"};
    TString selSteps[]={"Reco","2 leptons","Z-veto","#geq 2 jets","MET>40,0","=0 b-tags","=1 b-tags", "#geq 2 b-tags"};
    const size_t nselsteps=sizeof(selSteps)/sizeof(TString);
    results_["cutflow"]=formatPlot( baseDir.make<TH1F>("cutflow", ";Step; Events",nselsteps,0,nselsteps), 1,1,1,20,0,false,true,1,1,1);
    for(size_t istream=0; istream<sizeof(streams)/sizeof(TString); istream++)
      {
	TString cat=streams[istream];
	TFileDirectory newDir=baseDir.mkdir(cat.Data());

	results_[cat+"_cutflow"]=formatPlot( newDir.make<TH1F>(cat+"_cutflow", ";Step; Events",nselsteps,0,nselsteps), 1,1,1,20,0,false,true,1,1,1);
	for(int ibin=1; ibin<=((TH1F *)results_[cat+"_cutflow"])->GetXaxis()->GetNbins(); ibin++)
	  {
	    ((TH1F *)results_[cat+"_cutflow"])->GetXaxis()->SetBinLabel(ibin,selSteps[ibin-1]);
	    if(istream==0) ((TH1F *)results_["cutflow"])->GetXaxis()->SetBinLabel(ibin,selSteps[ibin-1]);
	  }
      
	//dilepton control
	results_[cat+"_dilepton_mass"]=formatPlot( newDir.make<TH1F>(cat+"_dilepton_mass", ";Invariant Mass(l,l') [GeV/c^{2}]; Events", 100, 0.,300.), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_dilepton_mass"]=formatPlot( newDir.make<TH1F>(cat+"_dilepton_mass", ";Invariant Mass(l,l') [GeV/c^{2}]; Events", 100, 0.,300.), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_dilepton_sumpt"]= formatPlot( newDir.make<TH1F>(cat+"_dilepton_sumpt", ";#Sigma |#vec{p}_{T}| [GeV/c]; Events", 100, 0.,300.), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_dilepton_pt"] = formatPlot( newDir.make<TH1F>(cat+"_dilepton_pt", ";|#Sigma #vec{p}_{T}| [GeV/c]; Events", 100, 0.,300.), 1,1,1,20,0,false,true,1,1,1);

	//vertex control
	results_[cat+"_vertex_sumpt"] = formatPlot( newDir.make<TH1F>(cat+"_vertex_sumpt", ";#Sigma_{tracks} p_{T} [GeV/c]; Events", 100, 0.,300.), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_othervertex_sumpt"] = formatPlot( newDir.make<TH1F>(cat+"_othervertex_sumpt", ";#Sigma_{tracks} p_{T} [GeV/c]; Events", 100, 0.,300.), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_vertex_pt"] = formatPlot( newDir.make<TH1F>(cat+"_vertex_pt", ";|#Sigma_{tracks} #vec{p}_{T}| [GeV/c]; Events", 100, 0.,300.), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_othervertex_pt"] = formatPlot( newDir.make<TH1F>(cat+"_othervertex_pt", ";|#Sigma_{tracks} #vec{p}_{T}| [GeV/c]; Events", 100, 0.,300.), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_ngoodvertex"] = formatPlot( newDir.make<TH1F>(cat+"_ngoodvertex", ";Vertices; Events", 25, 0.,25.), 1,1,1,20,0,false,true,1,1,1);
      
	//jets
	results_[cat+"_jetpt"]    = formatPlot( newDir.make<TH1F>(cat+"_jetpt",";p_{T} [GeV/c]; Jets",100,0,200), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_pujetpt"]    = formatPlot( newDir.make<TH1F>(cat+"_pujetpt",";p_{T} [GeV/c]; Jets",100,0,200), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_jeteta"]    = formatPlot( newDir.make<TH1F>(cat+"_jeteta",";#eta; Jets",100,-2.5,2.5), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_jetfassoc"]    = formatPlot( newDir.make<TH1F>(cat+"_jetfassoc",";f_{assoc}; Jets",100,0,1), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_pujetfassoc"]    = formatPlot( newDir.make<TH1F>(cat+"_pujetfassoc",";f_{assoc}; Jets",100,0,1), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_njets"]    = formatPlot( newDir.make<TH1F>(cat+"_njets",";Jet multiplicity; Events",4,0,4), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_npujets"]    = formatPlot( newDir.make<TH1F>(cat+"_npujets",";Jet multiplicity; Events",4,0,4), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_bmult"]    = formatPlot( newDir.make<TH1F>(cat+"_bmult",";b tag multiplicity (TCHEL); Events",4,0,4), 1,1,1,20,0,false,true,1,1,1);
	for(int ibin=1; ibin<=((TH1F *)results_[cat+"_njets"])->GetXaxis()->GetNbins(); ibin++)
	  {
	    TString ilabel(""); ilabel+=(ibin-1);
	    if(ibin==((TH1F *)results_[cat+"_njets"])->GetXaxis()->GetNbins()) ilabel="#geq"+ilabel;
	    ((TH1F *)results_[cat+"_njets"])->GetXaxis()->SetBinLabel(ibin,ilabel);
	    ((TH1F *)results_[cat+"_bmult"])->GetXaxis()->SetBinLabel(ibin,ilabel);
	  }

	results_[cat+"_btags"]  = formatPlot( newDir.make<TH1F>(cat+"_btags",";b tags (TCHE); Jets",100,-1,50), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_met"]               = formatPlot( newDir.make<TH1F>(cat+"_met", ";#slash{E}_{T} [GeV]; Events", 30,  0.,300.), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_jesmet"]            = formatPlot( newDir.make<TH1F>(cat+"_jesmet", ";#slash{E}_{T} (JES corrected) [GeV/c]; Events", 100,  0.,500.), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_jesmetnopu"]        = formatPlot( newDir.make<TH1F>(cat+"_jesmetnopu", ";#slash{E}_{T} (JES corrected, no PU) [GeV/c]; Events", 100,  0.,500.), 1,1,1,20,0,false,true,1,1,1);
	results_[cat+"_metsig"]            = formatPlot( newDir.make<TH1F>(cat+"_metsig", ";#slash{E}_{T} significance; Events", 100,  0.,100.), 1,1,1,20,0,false,true,1,1,1);

	//meant inclusive validation of JET/MET
	if(istream==0)
	  {
	    results_["jetpt"]    = formatPlot( baseDir.make<TH1F>("jetpt",";p_{T} [GeV/c]; Jets",100,0,200), 1,1,1,20,0,false,true,1,1,1);
	    results_["jeteta"]    = formatPlot( baseDir.make<TH1F>("jeteta",";#eta; Jets",100,-2.5,2.5), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetchhadenfrac"]    = formatPlot( baseDir.make<TH1F>("jetchhadenfrac",";f_{charged hadrons}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetneuthadenfrac"]    = formatPlot( baseDir.make<TH1F>("jetneuthadenfrac",";f_{neutral hadrons}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetchemenfrac"]    = formatPlot( baseDir.make<TH1F>("jetchemenfrac",";f_{charged electromagnetic}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetneutemenfrac"]    = formatPlot( baseDir.make<TH1F>("jetneutemenfrac",";f_{neutral electromagnetic}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetphoenfrac"]    = formatPlot( baseDir.make<TH1F>("jetphoenfrac",";f_{photons}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetmuenfrac"]    = formatPlot( baseDir.make<TH1F>("jetmuenfrac",";f_{muons}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetchhaden"]    = formatPlot( baseDir.make<TH1F>("jetchhaden",";Charged hadron components in jets [GeV/c]; Jets",50,0,200), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetneuthaden"]    = formatPlot( baseDir.make<TH1F>("jetneuthaden",";Neutral hadrons component in jets [GeV/c]; Jets",50,0,200), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetchemen"]    = formatPlot( baseDir.make<TH1F>("jetchemen",";Charged electromagnetic component in jets [GeV/c]; Jets",50,0,200), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetneutemen"]    = formatPlot( baseDir.make<TH1F>("jetneutemen",";Neutral electromagnetic component in jets [GeV/c]; Jets",50,0,200), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetphoen"]    = formatPlot( baseDir.make<TH1F>("jetphoen",";Photon component in jets [GeV/c]; Jets",50,0,200), 1,1,1,20,0,false,true,1,1,1);
	    results_["jetmuen"]    = formatPlot( baseDir.make<TH1F>("jetmuen",";Muon component in jets [GeV/c]; Jets",50,0,200), 1,1,1,20,0,false,true,1,1,1);
	    results_["met"]               = formatPlot( baseDir.make<TH1F>("met", ";#slash{E}_{T} [GeV]; Events", 30,  0.,300.), 1,1,1,20,0,false,true,1,1,1);
	    results_["sumet"]             = formatPlot( baseDir.make<TH1F>("sumet", ";#sum E_{T} [GeV]; Events", 100,  0.,1000.), 1,1,1,20,0,false,true,1,1,1);
	    results_["photonet"]          = formatPlot( baseDir.make<TH1F>("photonet", ";#sum_{photons} E_{T} [GeV]; Events", 100,  0.,500.), 1,1,1,20,0,false,true,1,1,1);
	    results_["neutralet"]         = formatPlot( baseDir.make<TH1F>("neutralet", ";#sum_{neutral} E_{T} [GeV]; Events", 100,  0.,500.), 1,1,1,20,0,false,true,1,1,1);
	    results_["electronet"]        = formatPlot( baseDir.make<TH1F>("electronet", ";#sum_{electron} E_{T} [GeV]; Events", 100,  0.,500.), 1,1,1,20,0,false,true,1,1,1);
	    results_["chhadet"]        = formatPlot( baseDir.make<TH1F>("chhadet", ";#sum_{charged hadron} E_{T} [GeV]; Events", 100,  0.,500.), 1,1,1,20,0,false,true,1,1,1);
	    results_["muonet"]        = formatPlot( baseDir.make<TH1F>("muonet", ";#sum_{muons} E_{T} [GeV]; Events", 100,  0.,500.), 1,1,1,20,0,false,true,1,1,1);
	    results_["photonetfrac"]          = formatPlot( baseDir.make<TH1F>("photonetfrac", ";f_{photons} E_{T} [GeV]; Events", 100,  0.,1.), 1,1,1,20,0,false,true,1,1,1);
	    results_["neutraletfrac"]         = formatPlot( baseDir.make<TH1F>("neutraletfrac", ";f_{neutral} E_{T} [GeV]; Events", 100,  0.,1.), 1,1,1,20,0,false,true,1,1,1);
	    results_["electronetfrac"]        = formatPlot( baseDir.make<TH1F>("electronetfrac", ";f_{electron} E_{T} [GeV]; Events", 100,  0.,1.), 1,1,1,20,0,false,true,1,1,1);
	    results_["chhadetfrac"]        = formatPlot( baseDir.make<TH1F>("chhadetfrac", ";f_{charged hadron} E_{T} [GeV]; Events", 100,  0.,1.), 1,1,1,20,0,false,true,1,1,1);
	    results_["muonetfrac"]        = formatPlot( baseDir.make<TH1F>("muonetfrac", ";f_{muons} E_{T} [GeV]; Events", 100,  0.,1.), 1,1,1,20,0,false,true,1,1,1);
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
    event.getByLabel(edm::InputTag("kt6PFJets:rho"),rhoH);
    float rho =*rhoH;

    edm::Handle<std::vector<reco::PFMET> > hPfMET;
    event.getByLabel(edm::InputTag("pfMet"), hPfMET);
    
    int selPath = (*(selInfo.product()))[0];
    int selStep = (*(selInfo.product()))[1];

    //require that a dilepton has benn selected
    if(selPath==0 or selStep<3) return;
    std::string istream="mumu";
    if(selPath==2) istream="ee";
    if(selPath==3) istream="emu";
    getHist(istream+"_cutflow")->Fill(1,weight);
    getHist("cutflow")->Fill(1,weight);
    
    //vertex quantities
    getHist(istream+"_ngoodvertex")->Fill(selVertices.size(),weight);
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
    getHist(istream+"_dilepton_sumpt")->Fill(lepton1P.pt()+lepton2P.pt(),weight);
    getHist(istream+"_dilepton_pt")->Fill(dileptonP.pt(),weight);
    getHist(istream+"_dilepton_mass")->Fill(dileptonP.mass(),weight);

    //Z+quarkonia veto
    bool isZCand(false);
    if(dileptonP.mass()<20) return;
    if( (istream=="ee" || istream=="mumu") && fabs(dileptonP.mass()-91)<15) isZCand=true;
    if(!isZCand)
      {
	getHist(istream+"_cutflow")->Fill(2,weight);
	getHist("cutflow")->Fill(2,weight);
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
	    getHist(istream+"_jetpt")->Fill(jet->pt(),weight);
	    getHist(istream+"_jeteta")->Fill(jet->eta(),weight);
	    getHist(istream+"_jetfassoc")->Fill( jet::fAssoc( jet.get(), &primVertex), weight );
	    getHist(istream+"_btags")->Fill(btag,weight);

	    //simulation validation
	    getHist("jetpt")->Fill(jet->pt(),weight);
	    getHist("jeteta")->Fill(jet->eta(),weight);
	    getHist("jetchhadenfrac")->Fill( jet->chargedHadronEnergyFraction(), weight );
	    getHist("jetneuthadenfrac")->Fill( jet->neutralHadronEnergyFraction(), weight );
	    getHist("jetchemenfrac")->Fill( jet->chargedEmEnergyFraction(),weight );
	    getHist("jetneutemenfrac")->Fill( jet->neutralEmEnergyFraction(),weight );
	    getHist("jetphoenfrac")->Fill( jet->photonEnergyFraction(),weight );
	    getHist("jetmuenfrac")->Fill( jet->muonEnergyFraction(),weight );
	    getHist("jetchhaden")->Fill( jet->chargedHadronEnergy(), weight );
	    getHist("jetneuthaden")->Fill( jet->neutralHadronEnergy(), weight );
	    getHist("jetchemen")->Fill( jet->chargedEmEnergy(),weight );
	    getHist("jetneutemen")->Fill( jet->neutralEmEnergy(),weight );
	    getHist("jetphoen")->Fill( jet->photonEnergy(),weight );
	    getHist("jetmuen")->Fill( jet->muonEnergy(),weight );
	  }
      }
    njets=selJets.size();
    if(!isZCand)
      {
	getHist(istream+"_njets")->Fill(njets,weight);
	getHist(istream+"_bmult")->Fill(nbjets,weight);
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
	    getHist(istream+"_pujetpt")->Fill(jet->pt(),weight);
	    getHist(istream+"_pujetfassoc")->Fill(jet::fAssoc(jet.get(),&primVertex),weight);
	  }
      }
    if(!isZCand) getHist(istream+"_npujets")->Fill(npujets,weight);


    //require two jets
    if(njets<2) return;
    if(!isZCand)
      {
	getHist(istream+"_cutflow")->Fill(3,weight);
	getHist("cutflow")->Fill(3,weight);
      }

    //base met kinematics
    const reco::PFMET &pfmet = (*hPfMET)[0];
    LorentzVectorCollection corMets = met::filter(pfmet,assocJets,puJets);

    LorentzVector metP = corMets[0];
    float metsig(-1),sumet(-1);
    float photonet(-1),neutralet(-1),electronet(-1),chhadet(-1),muonet(-1),
      photonetfrac(-1),neutraletfrac(-1),electronetfrac(-1),chhadetfrac(-1),muonetfrac(-1);
    metsig = pfmet.significance();
    sumet = pfmet.sumEt();
    photonet =pfmet.photonEt();
    neutralet=pfmet.neutralHadronEt();
    electronet=pfmet.electronEt();
    chhadet=pfmet.chargedHadronEt();
    muonet=pfmet.muonEt();
    photonetfrac =pfmet.photonEtFraction();
    neutraletfrac=pfmet.neutralHadronEtFraction();
    electronetfrac=pfmet.electronEtFraction();
    chhadetfrac=pfmet.chargedHadronEtFraction();
    muonetfrac=pfmet.muonEtFraction();

    //control  
    if(!isZCand)
      {
	getHist(istream+"_met")->Fill(metP.pt(),weight);
	getHist(istream+"_jesmet")->Fill(corMets[met::TYPEIMET].pt(),weight);
	getHist(istream+"_jesmetnopu")->Fill(corMets[met::CORRECTED_TYPEIMET].pt(),weight);
	getHist(istream+"_metsig")->Fill(metsig,weight);

	getHist("met")->Fill(metP.pt(),weight);
	getHist("sumet")->Fill(sumet,weight);
	getHist("photonet")->Fill(photonet,weight);
	getHist("neutralet")->Fill(neutralet,weight);
	getHist("electronet")->Fill(electronet,weight);
	getHist("chhadet")->Fill(chhadet,weight);
	getHist("muonet")->Fill(muonet,weight);
	getHist("photonetfrac")->Fill(photonetfrac,weight);
	getHist("neutraletfrac")->Fill(neutraletfrac,weight);
	getHist("electronetfrac")->Fill(electronetfrac,weight);
	getHist("chhadetfrac")->Fill(chhadetfrac,weight);
	getHist("muonetfrac")->Fill(muonetfrac,weight);
      }

    //require met for same flavor channels
    //if(metP.pt()<20) return;
    if( (istream=="ee" || istream=="mumu") && metP.pt()<40) return;
    if(!isZCand)
      {
	getHist(istream+"_cutflow")->Fill(4,weight);
	getHist("cutflow")->Fill(4,weight);
      }


    //b-tagged sample
    int btagbin= nbjets;
    if(btagbin>2) btagbin=2;
    if(!isZCand)
      {
	getHist(istream+"_cutflow")->Fill(5+btagbin,weight);
	getHist("cutflow")->Fill(5+btagbin,weight);
      }

    std::vector<reco::CandidatePtr> leptons;
    leptons.push_back(lepton1);
    leptons.push_back(lepton2);
    saveEvent(event,selPath,leptons,selJets,&pfmet,selVertices.size(),npuIT,rho,weight);
    
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
      ((TH1F *)getHist(streams[istream]+"_cutflow"))->Fill(0.,ctrHandle->value);
      if(istream==0) ((TH1F *)getHist("cutflow"))->Fill(0.,ctrHandle->value);
    }
}

//
void DileptonEventCleaner::saveEvent(const edm::Event& event, int evCat, std::vector<reco::CandidatePtr> &leptons, std::vector<const pat::Jet *> &jets, const reco::PFMET *met, 
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


