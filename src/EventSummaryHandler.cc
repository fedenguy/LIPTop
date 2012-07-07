#include "LIP/Top/interface/EventSummaryHandler.h"

using namespace std;

namespace top
{

  //
  EventSummaryHandler::EventSummaryHandler()
  {

  }
 
  //
  bool EventSummaryHandler::initTree(TTree *t, bool needsToRecreate)
  {
    if(t==0) return false;
    t_ = t;
    if(!needsToRecreate) return true;

    t_->Branch("run",        &evSummary_.run,    "run/I");
    t_->Branch("lumi",       &evSummary_.lumi,   "lumi/I");
    t_->Branch("event",      &evSummary_.event,  "event/I");
    t_->Branch("cat",        &evSummary_.cat,   "cat/I");
    t_->Branch("hasTrigger",   &evSummary_.hasTrigger,   "hasTrigger/O");
    t_->Branch("isSignal",   &evSummary_.isSignal,   "isSignal/O");
    t_->Branch("nvtx",       &evSummary_.nvtx,   "nvtx/I");
    t_->Branch("ngenpu",     &evSummary_.ngenpu,   "ngenpu/I");
    t_->Branch("ngenootpum1",     &evSummary_.ngenootpum1,   "ngenootpum1/I");
    t_->Branch("ngenootpup1",     &evSummary_.ngenootpup1,   "ngenootpup1/I");
    t_->Branch("rho",        &evSummary_.rho,    "rho/F");
    t_->Branch("weight",     &evSummary_.weight,   "weight/F");
    t_->Branch("normWeight",     &evSummary_.normWeight,   "normWeight/F");
    t_->Branch("xsecWeight",     &evSummary_.xsecWeight,   "xsecWeight/F");
    t_->Branch("nparticles", &evSummary_.nparticles, "nparticles/I");
    t_->Branch("px",         evSummary_.px,          "px[nparticles]/F");
    t_->Branch("py",         evSummary_.py,          "py[nparticles]/F");
    t_->Branch("pz",         evSummary_.pz,          "pz[nparticles]/F");
    t_->Branch("en",         evSummary_.en,          "en[nparticles]/F");
    t_->Branch("id",         evSummary_.id,          "id[nparticles]/I");
    t_->Branch("genid",      evSummary_.genid,       "genid[nparticles]/I"); 
    t_->Branch("genflav",    evSummary_.genflav,     "genflav[nparticles]/I"); 
    t_->Branch("info1",      evSummary_.info1,       "info1[nparticles]/F");
    t_->Branch("info2",      evSummary_.info2,       "info2[nparticles]/F");
    t_->Branch("info3",      evSummary_.info3,       "info3[nparticles]/F");
    t_->Branch("info4",      evSummary_.info4,       "info4[nparticles]/F");
    t_->Branch("info5",      evSummary_.info5,       "info5[nparticles]/F");
    t_->Branch("info6",      evSummary_.info6,       "info6[nparticles]/F");
    t_->Branch("info7",      evSummary_.info7,       "info7[nparticles]/F");
    t_->Branch("info8",      evSummary_.info8,       "info8[nparticles]/F");
    t_->Branch("info9",      evSummary_.info9,       "info9[nparticles]/F");
    
    t_->Branch("pthat",      &evSummary_.pthat,      "pthat/F");
    t_->Branch("genWeight",  &evSummary_.genWeight,  "genWeight/F");
    t_->Branch("qscale",     &evSummary_.qscale,     "qscale/F");
    t_->Branch("x1",         &evSummary_.x1,         "x1/F");
    t_->Branch("x2",         &evSummary_.x2,         "x2/F");
    t_->Branch("id1",        &evSummary_.id1,        "id1/I");
    t_->Branch("id2",        &evSummary_.id2,        "id2/I");
    t_->Branch("pdfWgts",    evSummary_.pdfWgts,     "pdfWgts[44]/F");
    
    t_->Branch("nmcparticles", &evSummary_.nmcparticles, "nmcparticles/I");
    t_->Branch("mcpx",         evSummary_.mcpx,          "mcpx[nmcparticles]/F");
    t_->Branch("mcpy",         evSummary_.mcpy,          "mcpy[nmcparticles]/F");
    t_->Branch("mcpz",         evSummary_.mcpz,          "mcpz[nmcparticles]/F");
    t_->Branch("mcen",         evSummary_.mcen,          "mcen[nmcparticles]/F");
    t_->Branch("mcid",         evSummary_.mcid,          "mcid[nmcparticles]/I");

    t_->Branch("nmeasurements",  &evSummary_.nmeasurements, "nmeasurements/I");
    t_->Branch("evmeasurements", evSummary_.evmeasurements, "evmeasurements[nmeasurements]/F");

    return true;
  }

  //
  bool EventSummaryHandler::attachToTree(TTree *t)
  {
    if(t==0) return false;
    t_ = t;

    t_->GetBranch("run")->SetAddress(&evSummary_.run);
    t_->GetBranch("lumi")->SetAddress(&evSummary_.lumi);
    t_->GetBranch("event")->SetAddress(&evSummary_.event);
    t_->GetBranch("cat")->SetAddress(&evSummary_.cat);
    if(t_->GetBranch("hasTrigger")) t_->GetBranch("hasTrigger")->SetAddress(&evSummary_.hasTrigger);
    t_->GetBranch("isSignal")->SetAddress(&evSummary_.isSignal);
    if(t_->GetBranch("nvtx"))   t_->GetBranch("nvtx")->SetAddress(&evSummary_.nvtx);
    t_->GetBranch("ngenpu")->SetAddress(&evSummary_.ngenpu);
    if(t_->GetBranch("ngenootpup1")) t_->GetBranch("ngenootpup1")->SetAddress(&evSummary_.ngenootpup1);
    if(t_->GetBranch("ngenootpum1")) t_->GetBranch("ngenootpum1")->SetAddress(&evSummary_.ngenootpum1);
    if(t_->GetBranch("rho"))    t_->GetBranch("rho")->SetAddress(&evSummary_.rho);
    t_->GetBranch("weight")->SetAddress(&evSummary_.weight);
    if(t_->GetBranch("normWeight") ) t_->GetBranch("normWeight")->SetAddress(&evSummary_.normWeight);
    if(t_->GetBranch("xsecWeight") ) t_->GetBranch("xsecWeight")->SetAddress(&evSummary_.xsecWeight);
    t_->GetBranch("nparticles")->SetAddress(&evSummary_.nparticles);
    t_->GetBranch("px")->SetAddress(evSummary_.px);
    t_->GetBranch("py")->SetAddress(evSummary_.py);
    t_->GetBranch("pz")->SetAddress(evSummary_.pz);
    t_->GetBranch("en")->SetAddress(evSummary_.en);
    t_->GetBranch("id")->SetAddress(evSummary_.id);
    t_->GetBranch("genid")->SetAddress(evSummary_.genid);
    if(t_->GetBranch("genflav")) t_->GetBranch("genflav")->SetAddress(evSummary_.genflav);
    t_->GetBranch("info1")->SetAddress(evSummary_.info1);
    t_->GetBranch("info2")->SetAddress(evSummary_.info2);
    t_->GetBranch("info3")->SetAddress(evSummary_.info3);
    t_->GetBranch("info4")->SetAddress(evSummary_.info4);
    t_->GetBranch("info5")->SetAddress(evSummary_.info5);
    if(t_->GetBranch("info6")) t_->GetBranch("info6")->SetAddress(evSummary_.info6);
    if(t_->GetBranch("info7")) t_->GetBranch("info7")->SetAddress(evSummary_.info7);
    if(t_->GetBranch("info8")) t_->GetBranch("info8")->SetAddress(evSummary_.info8);
    if(t_->GetBranch("info9")) t_->GetBranch("info9")->SetAddress(evSummary_.info9);

    if(t_->GetBranch("pthat"))     t_->GetBranch("pthat")->SetAddress(&evSummary_.pthat);
    if(t_->GetBranch("genWeight")) t_->GetBranch("genWeight")->SetAddress(&evSummary_.genWeight);
    if(t_->GetBranch("qscale"))    t_->GetBranch("qscale")->SetAddress(&evSummary_.qscale);
    if(t_->GetBranch("x1"))        t_->GetBranch("x1")->SetAddress(&evSummary_.x1);
    if(t_->GetBranch("x2"))        t_->GetBranch("x2")->SetAddress(&evSummary_.x2);
    if(t_->GetBranch("id1"))       t_->GetBranch("id1")->SetAddress(&evSummary_.id1);
    if(t_->GetBranch("id2"))       t_->GetBranch("id2")->SetAddress(&evSummary_.id2);
    if(t_->GetBranch("pdfWgts"))   t_->GetBranch("pdfWgts")->SetAddress(evSummary_.pdfWgts);

    if(t_->GetBranch("nmcparticles"))
      {
	t_->GetBranch("nmcparticles")->SetAddress(&evSummary_.nmcparticles);
	t_->GetBranch("mcpx")->SetAddress(evSummary_.mcpx);
	t_->GetBranch("mcpy")->SetAddress(evSummary_.mcpy);
	t_->GetBranch("mcpz")->SetAddress(evSummary_.mcpz);
	t_->GetBranch("mcen")->SetAddress(evSummary_.mcen);
	t_->GetBranch("mcid")->SetAddress(evSummary_.mcid);
      }

    if( t_->GetBranch("nmeasurements") )
      {
	t->GetBranch("nmeasurements")->SetAddress(&evSummary_.nmeasurements);
	t_->GetBranch("evmeasurements")->SetAddress(evSummary_.evmeasurements);
      }

    return true;
  }


  //
  void EventSummaryHandler::fillTree()
  {
    if(t_) t_->Fill();
  }

  //
  void EventSummaryHandler::fillTreeWithEvent(const EventSummary_t &ev)
  {
    std::vector<float> measurements;
    for(int i=0; i< ev.nmeasurements; i++) measurements.push_back(ev.evmeasurements[i]);
    fillTreeWithEvent(ev, measurements);
  }


  //
  void EventSummaryHandler::fillTreeWithEvent(const EventSummary_t &ev, std::vector<float> &addMeasurements)
  {
    evSummary_.run=ev.run; evSummary_.lumi=ev.lumi;     evSummary_.event=ev.event;
    evSummary_.cat=ev.cat; evSummary_.nvtx=ev.nvtx;     evSummary_.ngenpu=ev.ngenpu;
    evSummary_.rho=ev.rho; evSummary_.nparticles=ev.nparticles;
    evSummary_.isSignal = ev.isSignal;
    evSummary_.weight = ev.weight;
    evSummary_.normWeight = ev.normWeight;
    evSummary_.xsecWeight = ev.xsecWeight;

    evSummary_.nparticles=0;
    for(Int_t ipart=0; ipart<ev.nparticles; ipart++)
      {
	evSummary_.px[ipart]=ev.px[ipart];    
	evSummary_.py[ipart]=ev.py[ipart];
	evSummary_.pz[ipart]=ev.pz[ipart];     
	evSummary_.en[ipart]=ev.en[ipart];
	evSummary_.id[ipart]=ev.id[ipart];      
	evSummary_.genid[ipart]=ev.genid[ipart];
	evSummary_.genflav[ipart]=ev.genflav[ipart];
	evSummary_.info1[ipart]=ev.info1[ipart]; 
	evSummary_.info2[ipart]=ev.info2[ipart];
	evSummary_.info3[ipart]=ev.info3[ipart]; 
	evSummary_.info4[ipart]=ev.info4[ipart];
	evSummary_.info5[ipart]=ev.info5[ipart];
	evSummary_.info6[ipart]=ev.info6[ipart];
	evSummary_.info7[ipart]=ev.info7[ipart];
	evSummary_.info8[ipart]=ev.info8[ipart];     
	evSummary_.info9[ipart]=ev.info9[ipart];
	evSummary_.nparticles++;
      }

    evSummary_.nmcparticles=0;
    for(Int_t ipart=0; ipart<ev.nmcparticles; ipart++)
      {
	evSummary_.mcpx[ipart]=ev.mcpx[ipart];       evSummary_.mcpy[ipart]=ev.mcpy[ipart];
	evSummary_.mcpz[ipart]=ev.mcpz[ipart];       evSummary_.mcen[ipart]=ev.mcen[ipart];
	evSummary_.mcid[ipart]=ev.mcid[ipart];
	evSummary_.nmcparticles++;
      }

    evSummary_.nmeasurements=0;
    for(size_t i=0; i< addMeasurements.size(); i++) 
      {
	evSummary_.evmeasurements[i]=addMeasurements[i];
	evSummary_.nmeasurements++;
      }

    evSummary_.x1=ev.x1;    evSummary_.id1=ev.id1;
    evSummary_.x2=ev.x2;    evSummary_.id2=ev.id2;
    evSummary_.qscale=ev.qscale;
    for(size_t i=0; i<44; i++) evSummary_.pdfWgts[i-1]=ev.pdfWgts[i-1];
    
    fillTree();
  }


  //
  EventSummaryHandler::~EventSummaryHandler()
  {
  }


  //                                                                                                                                                                                                                                                                                                                        
  PhysicsEvent_t getPhysicsEventFrom(EventSummary_t &ev)
  {
    PhysicsEvent_t newev;
    newev.jets.clear(); newev.leptons.clear(); newev.topdecay.clear(); newev.antitopdecay.clear();
    newev.nbjets=0; newev.nljets=0;

    //get particles from the event
    for(Int_t ipart=0; ipart<ev.nparticles; ipart++)
      {
	LorentzVector p4(ev.px[ipart],ev.py[ipart],ev.pz[ipart],ev.en[ipart]);
	float pt(p4.pt());
	if(std::isnan(pt) || std::isinf(pt))continue;
	switch( ev.id[ipart] )
	  {
	  case 0:
	    newev.met = p4;
	    break;
	  case 1:
	    {
	      newev.jets.push_back( PhysicsObject_Jet(p4,ev.genid[ipart],ev.genflav[ipart],ev.info1[ipart],ev.info2[ipart],ev.info3[ipart],ev.info4[ipart], ev.info5[ipart], ev.info6[ipart], ev.info7[ipart], ev.info8[ipart],ev.info9[ipart]) );
	      newev.nbjets +=(fabs(ev.genflav[ipart])==5);
	      newev.nljets +=(fabs(ev.genflav[ipart])!=5);
	    }
	    break;
	  default:
	    newev.leptons.push_back( PhysicsObject_Lepton(p4,ev.id[ipart],ev.genid[ipart],ev.info1[ipart],ev.info2[ipart],ev.info3[ipart],ev.info4[ipart],ev.info5[ipart]) );
	    break;
	  }
      }
    if(newev.leptons.size()>1)  newev.dil = newev.leptons[0]+newev.leptons[1]; 


    //get gen level info
    for(Int_t ipart=0; ipart<ev.nmcparticles; ipart++)
      {
	LorentzVector p4(ev.mcpx[ipart],ev.mcpy[ipart],ev.mcpz[ipart],ev.mcen[ipart]);
	int id=ev.mcid[ipart];
	if(id==6)  newev.top=p4;
	if(fabs(id)<6 && id<0)  newev.topdecay.push_back( PhysicsObject(p4,id) );
	if(id<0 && (fabs(id)==11 || fabs(id)==13 || fabs(id)==15))  newev.topdecay.push_back( PhysicsObject(p4,id) );
      
	if(id==-6) newev.antitop=p4;
	if(fabs(id)<6 && id>0)  newev.antitopdecay.push_back( PhysicsObject(p4,id) );
	if(id>0 && (fabs(id)==11 || fabs(id)==13 || fabs(id)==15))  newev.antitopdecay.push_back( PhysicsObject(p4,id) );
      }

    return newev;
  }
  
}

//
top::EventSummaryHandler *eventHandlerFactory()
{
  return new top::EventSummaryHandler;
}
