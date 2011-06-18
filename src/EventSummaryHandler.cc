#include "LIP/Top/interface/EventSummaryHandler.h"

//
EventSummaryHandler::EventSummaryHandler()
{

}
 
//
bool EventSummaryHandler::initTree(TTree *t)
{
  if(t==0) return false;
  t_ = t;

  t_->Branch("run",        &evSummary_.run,    "run/I");
  t_->Branch("lumi",       &evSummary_.lumi,   "lumi/I");
  t_->Branch("event",      &evSummary_.event,  "event/I");
  t_->Branch("cat",        &evSummary_.cat,   "cat/I");
  t_->Branch("nvtx",       &evSummary_.nvtx,   "nvtx/I");
  t_->Branch("ngenpu",     &evSummary_.ngenpu,   "ngenpu/I");
  t_->Branch("rho",        &evSummary_.rho,    "rho/F");
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

  t_->Branch("nmcparticles", &evSummary_.nmcparticles, "nmcparticles/I");
  t_->Branch("mcpx",         evSummary_.mcpx,          "mcpx[nmcparticles]/F");
  t_->Branch("mcpy",         evSummary_.mcpy,          "mcpy[nmcparticles]/F");
  t_->Branch("mcpz",         evSummary_.mcpz,          "mcpz[nmcparticles]/F");
  t_->Branch("mcen",         evSummary_.mcen,          "mcen[nmcparticles]/F");
  t_->Branch("mcid",         evSummary_.mcid,          "mcid[nmcparticles]/I");

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
  if(t_->GetBranch("nvtx"))   t_->GetBranch("nvtx")->SetAddress(&evSummary_.nvtx);
  if(t_->GetBranch("ngenpu")) t_->GetBranch("ngenpu")->SetAddress(&evSummary_.ngenpu);
  if(t_->GetBranch("rho"))    t_->GetBranch("rho")->SetAddress(&evSummary_.rho);
  t_->GetBranch("weight")->SetAddress(&evSummary_.weight);
  t_->GetBranch("nparticles")->SetAddress(&evSummary_.nparticles);
  t_->GetBranch("px")->SetAddress(evSummary_.px);
  t_->GetBranch("py")->SetAddress(evSummary_.py);
  t_->GetBranch("pz")->SetAddress(evSummary_.pz);
  t_->GetBranch("en")->SetAddress(evSummary_.en);
  t_->GetBranch("id")->SetAddress(evSummary_.id);
  t_->GetBranch("genid")->SetAddress(evSummary_.genid);
  t_->GetBranch("info1")->SetAddress(evSummary_.info1);
  t_->GetBranch("info2")->SetAddress(evSummary_.info2);
  t_->GetBranch("info3")->SetAddress(evSummary_.info3);
  t_->GetBranch("info4")->SetAddress(evSummary_.info4);
  t_->GetBranch("info5")->SetAddress(evSummary_.info5);

  if(t_->GetBranch("nmcparticles"))
    {
      t_->GetBranch("nmcparticles")->SetAddress(&evSummary_.nmcparticles);
      t_->GetBranch("mcpx")->SetAddress(evSummary_.mcpx);
      t_->GetBranch("mcpy")->SetAddress(evSummary_.mcpy);
      t_->GetBranch("mcpz")->SetAddress(evSummary_.mcpz);
      t_->GetBranch("mcen")->SetAddress(evSummary_.mcen);
      t_->GetBranch("mcid")->SetAddress(evSummary_.mcid);
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
  evSummary_.run=ev.run; evSummary_.lumi=ev.lumi; evSummary_.event=ev.event;
  evSummary_.cat=ev.cat; evSummary_.nvtx=ev.nvtx; evSummary_.ngenpu=ev.ngenpu;
  evSummary_.rho=ev.rho; evSummary_.weight=ev.weight; evSummary_.nparticles=ev.nparticles;
  for(Int_t ipart=0; ipart<evSummary_.nparticles; ipart++)
    {
      evSummary_.px[ipart]=ev.px[ipart];       evSummary_.py[ipart]=ev.py[ipart];
      evSummary_.pz[ipart]=ev.pz[ipart];       evSummary_.en[ipart]=ev.en[ipart];
      evSummary_.id[ipart]=ev.id[ipart];       evSummary_.genid[ipart]=ev.genid[ipart];
      evSummary_.info1[ipart]=ev.info1[ipart]; evSummary_.info2[ipart]=ev.info2[ipart];
      evSummary_.info3[ipart]=ev.info3[ipart]; evSummary_.info4[ipart]=ev.info4[ipart];
      evSummary_.info5[ipart]=ev.info5[ipart];
    }

  for(Int_t ipart=0; ipart<evSummary_.nmcparticles; ipart++)
    {
      evSummary_.mcpx[ipart]=ev.mcpx[ipart];       evSummary_.mcpy[ipart]=ev.mcpy[ipart];
      evSummary_.mcpz[ipart]=ev.mcpz[ipart];       evSummary_.mcen[ipart]=ev.mcen[ipart];
      evSummary_.mcid[ipart]=ev.mcid[ipart];
    }

  fillTree();
}


//
EventSummaryHandler::~EventSummaryHandler()
{
}
