#include <iostream>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/KinAnalysis.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

//
int main(int argc, char* argv[])
{

  TString url="/afs/cern.ch/user/p/psilva/scratch0/CMSSW_4_1_3_patch2/src/LIP/Top/data/TTJets_madgraph_Spring11.root";
  TString dirname="evAnalyzer/data";
  int evStart=0;
  int evEnd=-1;
  TString scheme="std";

  TString jetDataUrl("${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data");
  gSystem->ExpandPathName(jetDataUrl);
  TString etaFileName = jetDataUrl + "/Spring10_EtaResolution_AK5PF.txt";
  TString phiFileName = jetDataUrl + "/Spring10_PhiResolution_AK5PF.txt";
  TString ptFileName  = jetDataUrl + "/Spring10_PtResolution_AK5PF.txt";
  TString uncFile = jetDataUrl+ "/Spring10_Uncertainty_AK5PF.txt";
  JetResolution stdEtaResol(etaFileName.Data(),false);
  JetResolution stdPhiResol(phiFileName.Data(),false);
  JetResolution stdPtResol(ptFileName.Data(),true); 
  JetCorrectionUncertainty jecUnc(uncFile.Data());

  //open the file and get directory
  TFile *file = TFile::Open(url);
  if(file==0) return -1;
  if(file->IsZombie()) return -1;
  EventSummaryHandler evSummaryHandler;
  if( !evSummaryHandler.attachToTree( (TTree *)file->Get(dirname) ) ) 
    {
      file->Close();
      return -1;
    }

  //check run range
  if(evEnd<0 || evEnd>evSummaryHandler.getEntries() ) 
    evEnd=evSummaryHandler.getEntries();
  if(evStart > evEnd ) 
    {
      file->Close();
      return -1;
    }
  
  //run the kin analysis
  //kin = KinAnalysis(scheme)
  for( int iev=evStart; iev<evEnd; iev++)
    {
      evSummaryHandler.getEntry(iev);
      cout << evSummaryHandler.evSummary_.run << " " <<  evSummaryHandler.evSummary_.lumi << " " << evSummaryHandler.evSummary_.event << endl;
      //kin.runOn(evSummaryHandler.getEvent(),ptresol,etaresol,phiresol,jecunc);
    }
  
  //return results;
}  
