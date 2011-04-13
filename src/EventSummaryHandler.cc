#include "LIP/Top/interface/EventSummaryHandler.h"

//
EventSummaryHandler::EventSummaryHandler()
{

}
 
//
void EventSummaryHandler::initTree(TTree *t)
{
  if(t==0) return;
  t_ = t;

  t_->Branch("run",        &evSummary_.run,    "run/I");
  t_->Branch("lumi",       &evSummary_.lumi,   "lumi/I");
  t_->Branch("event",      &evSummary_.event,  "event/I");
  t_->Branch("cat",        &evSummary_.cat,   "cat/I");
  t_->Branch("weight",     &evSummary_.weight,   "weight/F");
  t_->Branch("nparticles", &evSummary_.nparticles, "nparticles/I");
  t_->Branch("px",         evSummary_.px,          "px[nparticles]/F");
  t_->Branch("py",         evSummary_.py,          "py[nparticles]/F");
  t_->Branch("pz",         evSummary_.pz,          "pz[nparticles]/F");
  t_->Branch("en",         evSummary_.en,          "en[nparticles]/F");
  t_->Branch("id",         evSummary_.id,          "id[nparticles]/I");
  t_->Branch("genid",      evSummary_.genid,       "genid[nparticles]/I"); 
  t_->Branch("info1",      evSummary_.info1,       "info1[nparticles]/F");
  t_->Branch("info2",      evSummary_.info2,       "info2[nparticles]/F");
  t_->Branch("info3",      evSummary_.info3,       "info3[nparticles]/F");
  t_->Branch("info4",      evSummary_.info4,       "info4[nparticles]/F");
  t_->Branch("info5",      evSummary_.info5,       "info5[nparticles]/F");
}

//
void EventSummaryHandler::fillTree()
{
  if(t_) t_->Fill();
}

//
EventSummaryHandler::~EventSummaryHandler()
{
}
