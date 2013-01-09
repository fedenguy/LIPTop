#include <iostream>
#include <vector>

#include "TObjArray.h"
#include "TList.h"
#include "TH1D.h"
#include "TFile.h"
#include "TSystem.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

using namespace std;

//
void printHelp();

//                                                                                                                                                                                          
void printHelp()
{
  printf("--help    --> print this\n");
  printf("--in      --> plotter files to be used (csv list)\n");
}


//
int main(int argc, char* argv[])
{
  // load framework libraries                                                                                          
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  //configure
  std::vector<TString> url;
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos)              { printHelp();    return 0; }
      if(arg.find("--in")!=string::npos && i+1<argc)    { 
	TString val=argv[i+1];   
	TObjArray *tkns=val.Tokenize(",");
	for(int i=0; i<tkns->GetEntriesFast(); i++) url.push_back( gSystem->ExpandPathName( tkns->At(i)->GetName()) );
      }
    }
  if(url.size()==0) { printHelp(); return 0;}

  TString ch[]={"ee","emu","mumu"};
  TString dists[]={"jetlxy","jetmass"};
  TString flavs[]={"","udsg","c","b"};
  TString instrSysts[]={"jesup","jesdown","jerup","jerdown","puup","pudown"};
  TFile *outF=TFile::Open("lxydil.root","RECREATE");
  for(size_t ich=0; ich<sizeof(ch)/sizeof(TString); ich++)
    {
      outF->cd();
      TDirectory *outDir=outF->mkdir(ch[ich]);
      
      for(size_t i=0; i<url.size(); i++)
	{
	  TFile *inF=TFile::Open(url[i]);
	  if(inF==0) continue;
	  if(inF->IsZombie()) continue;
	  TIter next( inF->GetListOfKeys() );
	  TObject *obj=0;
	  while(true)
	    {
	      obj=next();
	      if(obj==0) break;

	      TString proc(obj->GetName());
	      proc.ReplaceAll("{","");
	      proc.ReplaceAll("}","");
	      proc.ReplaceAll("rightarrow","to");
	      proc.ReplaceAll("#","");
	      proc.ReplaceAll(" ","");
	      proc.ReplaceAll("tbart","ttbar");
	      proc.ToLower();
	      if(proc.Contains("ztotautausystdata")) continue;
	      TString mtop("0");
	      if(proc=="ttbar172.5") { cout << "Skipping " << proc << endl; continue; }
	      if(proc.Contains("ttbar")) 
		{
		  mtop="1725";
		  if(proc.Contains("163.5")) { mtop="1635"; proc.ReplaceAll("163.5",""); }
		  if(proc.Contains("166.5")) { mtop="1665"; proc.ReplaceAll("166.5",""); }
		  if(proc.Contains("169.5")) { mtop="1695"; proc.ReplaceAll("169.5",""); }
		  if(proc.Contains("172.5")) { mtop="1725"; proc.ReplaceAll("172.5",""); }
		  if(proc.Contains("175.5")) { mtop="1755"; proc.ReplaceAll("175.5",""); }
		  if(proc.Contains("178.5")) { mtop="1785"; proc.ReplaceAll("178.5",""); }
		  if(proc.Contains("181.5")) { mtop="1815"; proc.ReplaceAll("181.5",""); }
		}
	      	            
	      for(size_t idist=0; idist<sizeof(dists)/sizeof(TString); idist++)
		{
		  for(size_t iflav=0; iflav<sizeof(flavs)/sizeof(TString); iflav++)
		    {
		      TH1 *templ=(TH1 *)inF->Get( TString(obj->GetName()) + "/" + ch[ich]+"_"+flavs[iflav]+dists[idist] );
		      if(templ==0) continue;
		      TString outName( ch[ich]+"_"+flavs[iflav]+dists[idist]+"_"+proc+"_"+mtop );
		      templ->SetName(outName);
		      templ->SetDirectory(outDir);
		      outDir->cd();
		      templ->Write();

		      //get variations
		      for(size_t iSyst=0; iSyst<sizeof(instrSysts)/sizeof(TString); iSyst++)
			{
			  templ=(TH1 *)inF->Get( TString(obj->GetName()) + "/" + ch[ich]+"_"+flavs[iflav]+dists[idist]+instrSysts[iSyst] );
			  if(templ==0) continue;
			  TString outName( ch[ich]+"_"+flavs[iflav]+dists[idist]+"_"+proc+"syst"+instrSysts[iSyst]+"_"+mtop );
			  templ->SetName(outName);
			  templ->SetDirectory(outDir);
			  outDir->cd();
			  templ->Write();
			}

		    }
		}
	    }
	  
	  inF->Close();
	}
    }
  
  outF->Close();
}
