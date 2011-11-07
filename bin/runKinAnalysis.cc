#include <iostream>
#include <boost/shared_ptr.hpp>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/KinResultsHandler.h"
#include "LIP/Top/interface/KinAnalysis.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;
using namespace top;

//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  //check arguments
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }
  
  const edm::ParameterSet &kinProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("kinProcess");
  TString url=kinProcess.getParameter<std::string>("input");
  TString output=kinProcess.getParameter<std::string>("output");
  int evStart=kinProcess.getParameter<int>("evStart");
  int evEnd=kinProcess.getParameter<int>("evEnd");
  TString dirname = kinProcess.getParameter<std::string>("dirName");
  TString scheme = kinProcess.getParameter<std::string>("kinScheme");
  int maxTries = kinProcess.getParameter<int>("maxTries");
  int maxJetMult = kinProcess.getParameter<int>("maxJetMult");
  
  float mw=kinProcess.getParameter<double>("mw");
  float mb=kinProcess.getParameter<double>("mb");

  TString etaFileName = kinProcess.getParameter<std::string>("etaResolFileName"); gSystem->ExpandPathName(etaFileName);
  JetResolution stdEtaResol(etaFileName.Data(),false);

  TString phiFileName = kinProcess.getParameter<std::string>("phiResolFileName"); gSystem->ExpandPathName(phiFileName);
  JetResolution stdPhiResol(phiFileName.Data(),false);

  TString ptFileName  = kinProcess.getParameter<std::string>("ptResolFileName");  gSystem->ExpandPathName(ptFileName);
  JetResolution stdPtResol(ptFileName.Data(),true); 
  
  TString uncFile =  kinProcess.getParameter<std::string>("jesUncFileName"); gSystem->ExpandPathName(uncFile);
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
  //TString localFile("/tmp/KinAnalysis.root");  
  TString localFile(argv[1]);
  localFile.ReplaceAll("_cfg.py",".root"); 
  KinAnalysis kin(scheme,maxTries,maxJetMult,mw,mb,localFile,true);

  for( int iev=evStart; iev<evEnd; iev++)
    {
      evSummaryHandler.getEntry(iev);
      EventSummary_t &ev = evSummaryHandler.getEvent();
      kin.runOn(ev, &stdPtResol,&stdEtaResol,&stdPtResol,&jecUnc);
    }
  kin.endAnalysis();

  //move result to the output
  TString inputDir = gSystem->DirName( url );
  url.ReplaceAll(inputDir,"");
  url.ReplaceAll(".root","");
  url.ReplaceAll("/","");
  TString mkcmd( output.Contains("/castor") ? "rfmkdir -p" : "mkdir -p");
  gSystem->Exec(mkcmd + " " + output);
  output += "/" +  url + "/" + scheme;
  gSystem->Exec(mkcmd + " " + output);
  TString cpcmd( output.Contains("/castor") ? "rfcp" : "cp -v");
  TString outUrl(output+"/KinAnalysis_");  
  outUrl += evStart;
  outUrl += ".root";
  cpcmd += " " + localFile + " " + outUrl;
  gSystem->Exec(cpcmd);
  gSystem->Exec("rm " + localFile);
  std::cout  << "Result can be found @ :" << outUrl << std::endl;
}  
