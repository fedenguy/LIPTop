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
#include "TLorentzVector.h"

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
#include "CMGTools/HtoZZ2l2nu/interface/ReducedMETComputer.h"

class DileptonEventCleaner : public edm::EDAnalyzer 
{
public:
  explicit DileptonEventCleaner(const edm::ParameterSet& cfg);
  ~DileptonEventCleaner(){};
  virtual void analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup) ;
  void endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup);

private:

  void saveEvent(const edm::Event& event, int evCat, std::vector<reco::CandidatePtr> &leptons, std::vector<const pat::Jet *> &jets, const pat::MET *met, float weight);

  inline TH1 *getHist(TString key)
  {
    if(results_.find(key)==results_.end()) return 0;
    return (TH1 *)results_[key];
  }
  std::map<TString, TObject *>  results_;
  std::map<std::string, edm::ParameterSet> objConfig_;
  
  EventSummaryHandler summaryHandler_;
  ReducedMETComputer rmet_;
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
    TString selSteps[]={"Reco","2 leptons","Z-veto","=0 jets","=1 jet","#geq 2 jets","MET>30,0","=0 b-tags","=1 b-tags", "#geq 2 b-tags"};
    for(size_t istream=0; istream<sizeof(streams)/sizeof(TString); istream++)
    {
      TString cat=streams[istream];
      TFileDirectory newDir=baseDir.mkdir(cat.Data());

      const size_t nselsteps=sizeof(selSteps)/sizeof(TString);
      results_[cat+"_cutflow"]=formatPlot( newDir.make<TH1F>(cat+"_cutflow", ";Step; Events",nselsteps,0,nselsteps), 1,1,1,20,0,false,true,1,1,1);
      for(int ibin=1; ibin<=((TH1F *)results_[cat+"_cutflow"])->GetXaxis()->GetNbins(); ibin++)
	{
	  ((TH1F *)results_[cat+"_cutflow"])->GetXaxis()->SetBinLabel(ibin,selSteps[ibin-1]);
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
      results_[cat+"_jeteta"]    = formatPlot( newDir.make<TH1F>(cat+"_jeteta",";#eta; Jets",100,-2.5,2.5), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_jetfassoc"]    = formatPlot( newDir.make<TH1F>(cat+"_jetfassoc",";f_{assoc}; Jets",100,0,1), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_njets"]    = formatPlot( newDir.make<TH1F>(cat+"_njets",";Jet multiplicity; Events",4,0,4), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_chhadenfrac"]    = formatPlot( newDir.make<TH1F>(cat+"_chhadenfrac",";f_{charged hadrons}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_neuthadenfrac"]    = formatPlot( newDir.make<TH1F>(cat+"_neuthadenfrac",";f_{neutral hadrons}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_chemenfrac"]    = formatPlot( newDir.make<TH1F>(cat+"_chemenfrac",";f_{charged electromagnetic}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_neutemenfrac"]    = formatPlot( newDir.make<TH1F>(cat+"_neutemenfrac",";f_{neutral electromagnetic}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_phoenfrac"]    = formatPlot( newDir.make<TH1F>(cat+"_phoenfrac",";f_{_photons}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_muenfrac"]    = formatPlot( newDir.make<TH1F>(cat+"_muenfrac",";f_{muons}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_hfhadenfrac"]    = formatPlot( newDir.make<TH1F>(cat+"_hfhadenfrac",";f_{HF hadrons}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_hfemenfacr"]    = formatPlot( newDir.make<TH1F>(cat+"",";f_{HF electromagnetic}; Jets",50,0,1), 1,1,1,20,0,false,true,1,1,1);

      results_[cat+"_bmult"]    = formatPlot( newDir.make<TH1F>(cat+"_bmult",";b tag multiplicity (TCHEL); Events",4,0,4), 1,1,1,20,0,false,true,1,1,1);
      for(int ibin=1; ibin<=((TH1F *)results_[cat+"_njets"])->GetXaxis()->GetNbins(); ibin++)
	{
	  TString ilabel(""); ilabel+=(ibin-1);
	  if(ibin==((TH1F *)results_[cat+"_njets"])->GetXaxis()->GetNbins()) ilabel="#geq"+ilabel;
	  ((TH1F *)results_[cat+"_njets"])->GetXaxis()->SetBinLabel(ibin,ilabel);
	  ((TH1F *)results_[cat+"_bmult"])->GetXaxis()->SetBinLabel(ibin,ilabel);
	}
      results_[cat+"_btags"]             = formatPlot( newDir.make<TH1F>(cat+"_btags",";b tags (TCHE); Jets",100,-1,50), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_met"]               = formatPlot( newDir.make<TH1F>(cat+"_met", ";#slash{E}_{T} [GeV/c]; Events", 30,  0.,300.), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_metsig"]            = formatPlot( newDir.make<TH1F>(cat+"_metsig", ";#slash{E}_{T} significance; Events", 100,  0.,100.), 1,1,1,20,0,false,true,1,1,1);
    
      results_[cat+"_met"]               = formatPlot( newDir.make<TH1F>(cat+"_met", ";#slash{E}_{T} [GeV/c]; Events", 100,  0.,300.), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_rmet"]               = formatPlot( newDir.make<TH1F>(cat+"_rmet", ";red-#slash{E}_{T} [GeV/c]; Events", 100,  0.,300.), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_metsig"]            = formatPlot( newDir.make<TH1F>(cat+"_metsig", ";#slash{E}_{T} significance; Events", 100,  0.,100.), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_mT_individual"]     = formatPlot( newDir.make<TH1F>(cat+"_mT_individual",";Transverse mass(lepton,MET) [GeV/c^{2}]; Events",100,0,500), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_mT_corr"]           = formatPlot( newDir.make<TH2F>(cat+"_mT_corr",";Transverse mass(leading lepton,MET) [GeV/c^{2}];Transverse mass(trailer lepton,MET) [GeV/c^{2}]; Events",50,0,500,50,0,500), 1,1,1,20,0,false,true,1,1,1);
      results_[cat+"_mT_individualsum"]  = formatPlot( newDir.make<TH1F>(cat+"_mT_individualsum",";#Sigma Transverse mass(lepton,MET) [GeV/c^{2}]; Events",100,0,500), 1,1,1,20,0,false,true,1,1,1);
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
    
    int selPath = (*(selInfo.product()))[0];
    int selStep = (*(selInfo.product()))[1];

      //require that a dilepton has benn selected
    if(selPath==0 or selStep<3) return;
    std::string istream="mumu";
    if(selPath==2) istream="ee";
    if(selPath==3) istream="emu";
    getHist(istream+"_cutflow")->Fill(1,weight);
    
    //vertex quantities
    getHist(istream+"_ngoodvertex")->Fill(selVertices.size(),weight);
    const reco::Vertex &primVertex = (*(vertexHandle.product()))[0];

    //basic dilepton kinematics
    reco::CandidatePtr lepton1 = evhyp["leg1"];
    TLorentzVector lepton1P(lepton1->px(),lepton1->py(),lepton1->pz(), lepton1->energy());
    double lepton1pterr=dilepton::getPtErrorFor(lepton1);
    reco::CandidatePtr lepton2 = evhyp["leg2"];
    TLorentzVector lepton2P(lepton2->px(),lepton2->py(),lepton2->pz(), lepton2->energy());
    double lepton2pterr=dilepton::getPtErrorFor(lepton2);
    TLorentzVector dileptonP=lepton1P+lepton2P;
    getHist(istream+"_dilepton_sumpt")->Fill(lepton1P.Pt()+lepton2P.Pt(),weight);
    getHist(istream+"_dilepton_pt")->Fill(dileptonP.Pt(),weight);
    getHist(istream+"_dilepton_mass")->Fill(dileptonP.M(),weight);

    //Z+quarkonia veto
    if(dileptonP.M()<20) return;
    if( (istream=="ee" || istream=="mumu") && fabs(dileptonP.M()-91)<15) return;
    getHist(istream+"_cutflow")->Fill(2,weight);
    
    //count the jets in the event
    std::vector<reco::CandidatePtr> seljets= evhyp.all("jet");
    int njets( seljets.size() ), nbjets(0);
    std::vector<const pat::Jet *> selJets;
    std::vector<LorentzVector> jetmomenta;
    for (pat::eventhypothesis::Looper<pat::Jet> jet = evhyp.loopAs<pat::Jet>("jet"); jet; ++jet) {

      getHist(istream+"_jetpt")->Fill(jet->pt(),weight);
      getHist(istream+"_jeteta")->Fill(jet->eta(),weight);
      if(jet->pt()<30 || fabs(jet->eta())>2.5) continue;

      float fassoc=jet::fAssoc( jet.get(), &primVertex);
      getHist(istream+"_jetfassoc")->Fill( jet::fAssoc( jet.get(), &primVertex), weight );
      if(fassoc<0.1) continue;

      jetmomenta.push_back(jet->p4());

      float btag=jet->bDiscriminator("trackCountingHighEffBJetTags");
      if(btag>1.74) nbjets+=1; //loose point
      getHist(istream+"_btags")->Fill(btag,weight);

      getHist(istream+"_chhadenfrac")->Fill( jet->chargedHadronEnergyFraction(), weight );
      getHist(istream+"_neuthadenfrac")->Fill( jet->neutralHadronEnergyFraction(), weight );
      getHist(istream+"_chemenfrac")->Fill( jet->chargedEmEnergyFraction(),weight );
      getHist(istream+"_neutemenfrac")->Fill( jet->neutralEmEnergyFraction(),weight );
      getHist(istream+"_phoenfrac")->Fill( jet->photonEnergyFraction(),weight );
      getHist(istream+"_muenfrac")->Fill( jet->muonEnergyFraction(),weight );
      getHist(istream+"_hfhadenfrac")->Fill( jet->HFHadronEnergyFraction() ,weight);
      getHist(istream+"_hfemenfacr")->Fill( jet->HFEMEnergyFraction(),weight );
      
      selJets.push_back( jet.get() );
    }
    getHist(istream+"_njets")->Fill(njets,weight);
    getHist(istream+"_bmult")->Fill(nbjets,weight);

    //require two jets
    if(njets==0) getHist(istream+"_cutflow")->Fill(3,weight);
    if(njets==1) getHist(istream+"_cutflow")->Fill(4,weight);
    if(njets<2) return;
    getHist(istream+"_cutflow")->Fill(5,weight);

    //base met kinematics
    const pat::MET *themet=evhyp.getAs<pat::MET>("met");
    TLorentzVector metP(themet->px(),themet->py(),0,themet->pt());
    float metsig(-1);
    try{
      const reco::MET *origMet = dynamic_cast<const reco::MET *>(themet->originalObject());
      metsig=origMet->significance();
    }catch(std::exception &e){
      metsig=themet->significance();
    }
    float dphil2met[]={ fabs(metP.DeltaPhi(lepton1P)), fabs(metP.DeltaPhi(lepton2P)) };
    float mTlmet[]={ TMath::Sqrt(2*metP.Pt()*lepton1P.Pt()*(1-TMath::Cos(dphil2met[0]))) ,   TMath::Sqrt(2*metP.Pt()*lepton2P.Pt()*(1-TMath::Cos(dphil2met[1]))) };

    //reduced met                                                                                                                                                                                             
    rmet_.compute(lepton1->p4(),lepton1pterr,
                  lepton2->p4(),lepton2pterr,
                  jetmomenta,
                  themet->p4()
                  );
    float reducedMET=rmet_.reducedMET();
    
    getHist(istream+"_met")->Fill(metP.Pt(),weight);
    getHist(istream+"_rmet")->Fill(reducedMET,weight);
    getHist(istream+"_metsig")->Fill(metsig,weight);
    getHist(istream+"_mT_individual")->Fill(mTlmet[0],weight);
    getHist(istream+"_mT_individual")->Fill(mTlmet[1],weight);
    if(lepton1P.Pt()>lepton2P.Pt()) ((TH2 *)getHist(istream+"_mT_corr"))->Fill(mTlmet[0],mTlmet[1],weight);
    else ((TH2 *)getHist(istream+"_mT_corr"))->Fill(mTlmet[1],mTlmet[0],weight);
    getHist(istream+"_mT_individualsum")->Fill(mTlmet[0]+mTlmet[1],weight);

    //require met for same flavor channels
    //if(metP.Pt()<20) return;
    if( (istream=="ee" || istream=="mumu") && metP.Pt()<30) return;
    getHist(istream+"_cutflow")->Fill(6,weight);

    //b-tagged sample
    int btagbin= nbjets;
    if(btagbin>2) btagbin=2;
    getHist(istream+"_cutflow")->Fill(7+btagbin,weight);

    std::vector<reco::CandidatePtr> leptons;
    leptons.push_back(lepton1);
    leptons.push_back(lepton2);
    saveEvent(event,selPath,leptons,selJets,themet,weight);
    
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
    ((TH1F *)getHist(streams[istream]+"_cutflow"))->Fill(0.,ctrHandle->value);
}

//
void DileptonEventCleaner::saveEvent(const edm::Event& event, int evCat, std::vector<reco::CandidatePtr> &leptons, std::vector<const pat::Jet *> &jets, const pat::MET *met, float weight)
{
  //save event header
  summaryHandler_.evSummary_.run=event.id().run();
  summaryHandler_.evSummary_.lumi=event.luminosityBlock();
  summaryHandler_.evSummary_.event=event.id().event();
  summaryHandler_.evSummary_.cat=evCat; 
  summaryHandler_.evSummary_.weight=weight; 
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


