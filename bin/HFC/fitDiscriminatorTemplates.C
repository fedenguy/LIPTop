#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "RooNLLVar.h"
#include "RooSimultaneous.h"
#include "RooMinuit.h"

#include "TList.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TSystem.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TFile.h"
#include "TObjArray.h"

#include<sstream>
#include<vector>
#include<map>

using namespace RooFit;
using namespace std;

struct Shape_t
{
  TH2F *nominalH;
  std::map<TString,TH2F *> varH; 
  std::map<TString, std::vector<TH1F *> > qcdTemplatesH;
};

void fitData(TH1 *preTagData, TH1 *postTagData, TH1 *vetoData,
	     std::vector<TH1 *> preTagTemplates, std::vector<TH1 *> postTagTemplates, std::vector<TH1 *> vetoTemplates,
	     TString tag, TString tagLabel);

struct FitResults_t
{
  Double_t pretag, pretagerr,tag,tagerr;
  Double_t mcpretag, mcpretagerr,mctag,mctagerr;
  Double_t sfb,sfberr;
};

std::map<TString, FitResults_t> results;

void showShape(Shape_t &data,TString tag, TString tagLabel, int firstBin,int lastBin);
TGraphAsymmErrors *getUncertaintyBandFrom(std::vector<TH1 *> &varColl);


//
void fitDiscriminatorTemplates(TString url="plotter.root",TString btagger="jp",TString postTagWP="csvL",TString qcdUrl="QCDTemplates_PtHatBins.root")
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TString systVars[]={"","jerdown","jerup","pudown","puup","jesdown","jesdown"}; //,"umetdown","umetup"};
  size_t nSystVars=sizeof(systVars)/sizeof(TString); 
  TString flavors[]={"","b","c","udsg"};
  const size_t nFlavors=sizeof(flavors)/sizeof(TString);
  TString procs[]={"t#bar{t}","data"};
  //  {"Di-boson","other t#bar{t}","t#bar{t}","Single top","QCD","W#rightarrow l#nu","Z#rightarrow ll","data"};
  const size_t nProcs=sizeof(procs)/sizeof(TString);
  TString postTag[]={"",postTagWP,postTagWP+"v"};
  const size_t npostTags=sizeof(postTag)/sizeof(TString);

  //get the histos from the file
  TFile *fIn    = TFile::Open(url);
  TFile *fQCDIn = TFile::Open(qcdUrl);

  Shape_t              preTagData,      postTagData,       vetoData;
  std::vector<Shape_t> preTagTemplates, postTagTemplates,  vetoTemplates;
  for(size_t iflav=0; iflav<nFlavors; iflav++)
    {
      for(size_t ipt=0; ipt<npostTags; ipt++)
	{
	  Shape_t templateShape,dataShape;
	  
	  //get qcd templates for this flavor (temporary will be fixed)
	  /*
	  if(iflav>0)
	    {
	      TString flavName("B");
	      if(flavors[iflav]=="c") flavName="C";
	      if(flavors[iflav]=="udsg") flavName="L";
	      TString hname(postTag[ipt]);
	      if(ipt>0) hname += "_";  
	      hname+=flavName; 
	      hname.ToUpper();
	      hname = "Discr_"+hname+"_Pt";
	      TString qcdcats[]={"3050","5080","80120","120"};
	      for(size_t iqcdcat=0; iqcdcat<sizeof(qcdcats)/sizeof(TString); iqcdcat++)
		{
		  TH1F *h= (TH1F *)fQCDIn->Get(hname+qcdcats[iqcdcat]+"_Eta24");
		  h->Sumw2();
		  h->SetDirectory(0);
		  h->Scale(1./h->Integral());
		  h->SetLineWidth(2);
		  h->SetLineColor(8);
		  templateShape.qcdTemplatesH[flavors[iflav]].push_back(h);
		}
	    }
	  */

	  //get templates
	  for(size_t isyst=0; isyst<nSystVars; isyst++)
	    {
	      //integrate over all the processes
	      TH2F *mcFlavH=0,*dataFlavH=0;
	      for(size_t iproc=0; iproc<nProcs; iproc++)
		{
		  TH2F *ih=(TH2F *) fIn->Get(procs[iproc]+"/"+btagger+flavors[iflav]+systVars[isyst]+postTag[ipt]);
		  if(ih==0) continue;
		  if(procs[iproc]!="data")
		    {
		      //if non-existing clone it, otherwise add it
		      if(mcFlavH==0) { 
			mcFlavH=(TH2F *)ih->Clone(btagger+flavors[iflav]+systVars[isyst]+postTag[ipt]);
			mcFlavH->SetDirectory(0); 
		      }
		      else  { 
			mcFlavH->Add(ih);
		      }
		    }
		  else if(isyst==0) {
		    dataFlavH = (TH2F *)ih->Clone("data"+btagger+flavors[iflav]+systVars[isyst]+postTag[ipt]); 
		    dataFlavH->SetDirectory(0); 
		    dataShape.nominalH=dataFlavH;
		  }
		}

	      //save the template
	      if(isyst==0) templateShape.nominalH=mcFlavH; 
	      else         templateShape.varH[systVars[isyst]]=mcFlavH;                            
	    }
	  
	  //add the MC statistics also as systematic for completion
	  TH2F *mcStatUpH  =(TH2F *)dataShape.nominalH->Clone(btagger+flavors[iflav]+"mcstatup"+postTag[ipt]);   
	  mcStatUpH->SetDirectory(0);   mcStatUpH->Reset("ICE");
	  TH2F *mcStatDownH=(TH2F *)dataShape.nominalH->Clone(btagger+flavors[iflav]+"mcstatdown"+postTag[ipt]); 
	  mcStatDownH->SetDirectory(0); mcStatDownH->Reset("ICE");
	  for(int xbin=1; xbin<= templateShape.nominalH->GetNbinsX(); xbin++)
	    {
	      for(int ybin=1; ybin<= templateShape.nominalH->GetNbinsY(); ybin++)
		{
		  Double_t val=templateShape.nominalH->GetBinContent(xbin,ybin);
		  Double_t err=templateShape.nominalH->GetBinError(xbin,ybin);
		  mcStatUpH->SetBinContent(xbin,ybin,val+err);
		  mcStatDownH->SetBinContent(xbin,ybin,val-err);
		}
	    }
	  templateShape.varH["mcstatup"]   = mcStatUpH;
	  templateShape.varH["mcstatdown"] = mcStatDownH;
	  
	  //save the shape
	  if(ipt==0)       preTagTemplates.push_back(templateShape);
	  else if(ipt==1)  postTagTemplates.push_back(templateShape);
	  else             vetoTemplates.push_back(templateShape);
	  if(iflav>0) continue;
	  if(ipt==0)      preTagData=dataShape;
	  else if(ipt==1) postTagData=dataShape;
	  else            vetoData=dataShape;
	}
    }
  fIn->Close();
  fQCDIn->Close();
  
  //now do the fitting
  //inclusive
  TH1 *preTagDataH  = preTagData.nominalH->ProjectionX("pretagincdata",1,1);
  TH1 *postTagDataH = postTagData.nominalH->ProjectionX("posttagincdata",1,1);
  TH1 *vetoDataH    = vetoData.nominalH->ProjectionX("vetoincdata",1,1);
  std::vector<TH1 *> preTagTemplatesH, postTagTemplatesH, vetoTemplatesH;
  for(size_t iflav=1; iflav<preTagTemplates.size(); iflav++)
    {
      TH1 *mcH=preTagTemplates[iflav].nominalH->ProjectionX(flavors[iflav]+"pretag",1,1);
      preTagTemplatesH.push_back(mcH);
      mcH=postTagTemplates[iflav].nominalH->ProjectionX(flavors[iflav]+"posttag",1,1);
      postTagTemplatesH.push_back(mcH);
      mcH=vetoTemplates[iflav].nominalH->ProjectionX(flavors[iflav]+"veto",1,1);
      vetoTemplatesH.push_back(mcH);

      //display the shape 
      //showShape(allMCShapes[iflav],flavors[iflav]+"inclusive",flavors[iflav]+ " inclusive",1,1);
    }
  
  fitData(preTagDataH,postTagDataH, vetoDataH, preTagTemplatesH,postTagTemplatesH, vetoTemplatesH,"inclusive","inclusive");
  //fitData(preTagDataH,postTagDataH,0, preTagTemplatesH,postTagTemplatesH, vetoTemplatesH,"inclusive","inclusive");

  /*
  //per categories
  for(int ybin=1; ybin<=theDataShape.nominalH->GetYaxis()->GetNbins(); ybin++)
  {
  dataH=theDataShape.nominalH->ProjectionX("data",ybin,ybin);
  TString yrang=theDataShape.nominalH->GetYaxis()->GetBinLabel(ybin);
  TString yLabel=yrang; yLabel.ReplaceAll("to"," - "); yLabel.ReplaceAll("Inf","#infty"); 
  templH.clear();
  for(size_t iflav=1; iflav<allMCShapes.size(); iflav++)
  {
  TH1 *mcH=allMCShapes[iflav].nominalH->ProjectionX(flavors[iflav],ybin,ybin);
  templH.push_back(mcH);
  
  //display the shape 
  TString yFlavLabel=flavors[iflav]+ " ("+yLabel+")";
  showShape(allMCShapes[iflav],flavors[iflav]+postTag[ipt]+yrang,yFlavLabel,ybin,ybin);
  }
  //fitData(dataH,templH,yrang,yLabel);
  }
  }
  
  cout << "--------------------------------------" << endl
       << " Fit results for " << postTagWP << endl;
  TGraphErrors *ptgr=new TGraphErrors;
  TGraphErrors *etagr=new TGraphErrors;
  TGraphErrors *vtxgr=new TGraphErrors;
  Double_t incsfb,incsfberr;
  for(std::map<TString,FitResults_t>::iterator it= results.begin(); it != results.end(); it++)
    {
      float mceffb=it->second.mctag/it->second.mcpretag;
      float effb=it->second.tag/it->second.pretag;
      float effberr=sqrt(pow(it->second.tag*it->second.pretagerr,2)+pow(it->second.pretag*it->second.tagerr,2))/pow(it->second.pretag,2);
      float sfb=effb/mceffb;       it->second.sfb=sfb;
      float sfberr=effberr/mceffb; it->second.sfberr=sfberr;

      if(sfberr>sfb)
	{
	  cout << "[Warning] " << it->first << " did not converge properly" << endl;
	  continue;
	}

      printf("%s & \t %3.3lf & \t & %3.3lf $\\pm$ %3.3lf & \t %3.3lf $\\pm$ %3.3lf\n",
	     it->first.Data(), mceffb,effb,effberr,sfb,sfberr);
      
      TGraphErrors *gr=0;
      Double_t x,ex;
      if     (it->first=="30to50")   { gr=ptgr;  x=0.5*(30+50);   ex=50-x; }
      else if(it->first=="50to80")   { gr=ptgr;  x=0.5*(50+80);   ex=80-x; }
      else if(it->first=="80to120")  { gr=ptgr;  x=0.5*(80+120);  ex=120-x; }
      else if(it->first=="120toInf") { gr=ptgr;  x=0.5*(120+200); ex=200-x; }
      else if(it->first=="0to0.5")   { gr=etagr; x=0.5*(0+0.5);   ex=0.5-x; }
      else if(it->first=="0.5to1.0") { gr=etagr; x=0.5*(1.0+0.5); ex=1.0-x; }
      else if(it->first=="1.0to1.5") { gr=etagr; x=0.5*(1.0+1.5); ex=1.5-x; }
      else if(it->first=="1.5to2.5") { gr=etagr; x=0.5*(1.5+2.5); ex=2.5-x; }
      else if(it->first=="1.5to2.5") { gr=etagr; x=0.5*(1.5+2.5); ex=2.5-x; }
      else if(it->first=="0to10")    { gr=vtxgr; x=0.5*(10+0);    ex=10-x; }
      else if(it->first=="11to14")   { gr=vtxgr; x=0.5*(15+10);   ex=15-x; }
      else if(it->first=="15to18")   { gr=vtxgr; x=0.5*(19+15);   ex=19-x; }
      else if(it->first=="18toInf")  { gr=vtxgr; x=0.5*(25+19);   ex=25-x; }    
      else                           { incsfb=sfb; incsfberr=sfberr; }
      if(gr==0) continue;
      Int_t np=gr->GetN();
      gr->SetPoint(np,x,sfb); gr->SetPointError(np,ex,sfberr);
    }
  cout << "--------------------------------------" << endl;


  TCanvas *fc= new TCanvas("fc","fc",600,600);
  fc->Divide(1,3);
  
  TPad *p=(TPad *)fc->cd(1); 
  p->SetTopMargin(0.05); p->SetBottomMargin(0.2);
  ptgr->Draw("ap");  ptgr->SetMarkerStyle(20);
  ptgr->GetXaxis()->SetTitle("p_{T} [GeV/c]"); ptgr->GetXaxis()->SetTitleSize(0.08);  ptgr->GetXaxis()->SetLabelSize(0.07); 
  ptgr->GetYaxis()->SetTitle("SF_{b}");        ptgr->GetYaxis()->SetTitleSize(0.08);  ptgr->GetYaxis()->SetLabelSize(0.07);  ptgr->GetYaxis()->SetTitleOffset(0.5);
  ptgr->GetYaxis()->SetRangeUser(0.7,1.3);
  TLine *l= new TLine(30,incsfb+incsfberr,200,incsfb+incsfberr); l->SetLineColor(kRed); l->SetLineStyle(7); l->Draw();
  l= new TLine(30,incsfb-incsfberr,200,incsfb-incsfberr); l->SetLineColor(kRed); l->SetLineStyle(7); l->Draw();
  
  TPaveText *txt= new TPaveText(0.8,0.8,0.9,0.9,"brNDC");
  txt->AddText(postTagWP);
  txt->SetTextFont(42); txt->SetBorderSize(0); txt->SetFillStyle(0);
  txt->Draw();


  p=(TPad *)fc->cd(2);
  p->SetTopMargin(0); p->SetBottomMargin(0.2);
  etagr->Draw("ap"); etagr->SetMarkerStyle(20);
  etagr->GetXaxis()->SetTitle("|#eta|");       etagr->GetXaxis()->SetTitleSize(0.08);  etagr->GetXaxis()->SetLabelSize(0.07);  
  etagr->GetYaxis()->SetTitle("SF_{b}");       etagr->GetYaxis()->SetTitleSize(0.08);  etagr->GetYaxis()->SetLabelSize(0.07);  etagr->GetYaxis()->SetTitleOffset(0.5);
  etagr->GetYaxis()->SetRangeUser(0.7,1.3);
  l= new TLine(0,incsfb+incsfberr,2.5,incsfb+incsfberr); l->SetLineColor(kRed); l->SetLineStyle(7); l->Draw();
  l= new TLine(0,incsfb-incsfberr,2.5,incsfb-incsfberr); l->SetLineColor(kRed); l->SetLineStyle(7); l->Draw();

  p=(TPad *)fc->cd(3); 
  p->SetTopMargin(0); p->SetBottomMargin(0.2);
  vtxgr->Draw("ap"); vtxgr->SetMarkerStyle(20);
  vtxgr->GetXaxis()->SetTitle("Vertices");     vtxgr->GetXaxis()->SetTitleSize(0.08);  vtxgr->GetXaxis()->SetLabelSize(0.07); 
  vtxgr->GetYaxis()->SetTitle("SF_{b}");       vtxgr->GetYaxis()->SetTitleSize(0.08);  vtxgr->GetYaxis()->SetLabelSize(0.07); vtxgr->GetYaxis()->SetTitleOffset(0.5);
  vtxgr->GetYaxis()->SetRangeUser(0.7,1.3);
  l= new TLine(0,incsfb+incsfberr,25,incsfb+incsfberr); l->SetLineColor(kRed); l->SetLineStyle(7); l->Draw();
  l= new TLine(0,incsfb-incsfberr,25,incsfb-incsfberr); l->SetLineColor(kRed); l->SetLineStyle(7); l->Draw();

  */
  

}



//
void showShape(Shape_t &mc,TString tag, TString tagLabel, int firstBin,int lastBin)
{

  TCanvas *c = new TCanvas(tag+"c",tag+"c",600,600); c->cd();

  //project the inclusive shape
  c->cd();
  TPad* t1 = new TPad("t1","t1", 0.0, 0.20, 1.0, 1.0);  t1->Draw();  t1->cd();
  
  //project the shape
  TH1 *mcIncH=mc.nominalH->ProjectionX(tag+"incmc",firstBin,lastBin);       
  mcIncH->SetDirectory(0);  

  //build the total uncertainty band
  std::vector<TH1 *> allVarsH;
  for(std::map<TString,TH2F *>::iterator hIt = mc.varH.begin(); hIt != mc.varH.end(); hIt++)
    {
      TH1 *mcVarH = hIt->second->ProjectionX(TString(hIt->second->GetName())+"incmc",firstBin,lastBin); 
      mcVarH->Add(mcIncH,-1);
      mcVarH->SetDirectory(0);
      allVarsH.push_back(mcVarH);
    }
  TGraphAsymmErrors *mcErrGr=getUncertaintyBandFrom(allVarsH);  mcErrGr->SetMarkerStyle(1);    mcErrGr->SetFillStyle(3002);
  TGraphAsymmErrors *absMcErrGr = new TGraphAsymmErrors;        absMcErrGr->SetMarkerStyle(1); absMcErrGr->SetFillStyle(3002);
  TGraphAsymmErrors *relMcErrGr = new TGraphAsymmErrors;        relMcErrGr->SetMarkerStyle(1); relMcErrGr->SetFillStyle(3003);
  for(int ip=0; ip<mcErrGr->GetN(); ip++)
    {
      Double_t ix,ixerr,iy,iyerrhigh,iyerrlow;
      mcErrGr->GetPoint(ip,ix,iy);
      iy=mcIncH->GetBinContent(ip+1);
      ixerr=mcErrGr->GetErrorX(ip);
      iyerrhigh=mcErrGr->GetErrorYhigh(ip);
      iyerrlow=mcErrGr->GetErrorYlow(ip);

      mcErrGr->SetPoint(ip,ix,1);
      absMcErrGr->SetPoint(ip,ix,mcIncH->GetBinContent(ip+1));
      absMcErrGr->SetPointError(ip,ixerr,ixerr,iyerrlow,iyerrhigh);
     
      relMcErrGr->SetPoint(ip,ix,1.0);
      if(iy!=0)      relMcErrGr->SetPointError(ip,ixerr,ixerr,iyerrlow/iy,iyerrhigh/iy);
      else           relMcErrGr->SetPointError(ip,ixerr,ixerr,0,0);
    }

  //now draw
  mcIncH->Draw("hist");    mcIncH->GetYaxis()->SetTitle("Jets"); mcIncH->GetXaxis()->SetTitle("Discriminator");
  absMcErrGr->Draw("e2");
  
  TPaveText* T = new TPaveText(0.1,0.995,0.9,0.96, "NDC");
  T->SetFillColor(0);
  T->SetFillStyle(0);  T->SetLineColor(0);
  T->SetTextAlign(22);
  char Buffer[1024]; sprintf(Buffer, "CMS simulation, %s", tagLabel.Data());
  T->AddText(Buffer);
  T->Draw("same");

  //project the variations for the mc
  c->cd();
  TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, 0.2);     t2->Draw(); t2->cd();
  float yscale = (1.0-0.2)/(0.18-0);

  //draw the total uncertainty band
  relMcErrGr->Draw("ae2");
  relMcErrGr->SetMinimum(0.5);
  relMcErrGr->SetMaximum(1.5);
  relMcErrGr->GetXaxis()->SetTitleOffset(1.3);
  relMcErrGr->GetXaxis()->SetLabelSize(0.033*yscale);
  relMcErrGr->GetXaxis()->SetTitleSize(0.036*yscale);
  relMcErrGr->GetXaxis()->SetTickLength(0.03*yscale);
  relMcErrGr->GetYaxis()->SetTitleOffset(0.3);
  relMcErrGr->GetYaxis()->SetNdivisions(5);
  relMcErrGr->GetYaxis()->SetLabelSize(0.033*yscale);
  relMcErrGr->GetYaxis()->SetTitleSize(0.036*yscale);
  relMcErrGr->GetXaxis()->SetRangeUser(mcIncH->GetXaxis()->GetXmin(),mcIncH->GetXaxis()->GetXmax());

  //draw the breakup of the uncertainties
  int varColors[] = {1,kRed,kRed,kBlue,kBlue,kGreen+3,kGreen+3,kGreen,kGreen, kBlue+2, kBlue+2};
  int varStyles[] = {1,1,   7,   1,    7,    1,       7,       1,     7     , 1,       7};
  for(size_t ivar=0; ivar<allVarsH.size(); ivar++)
    {
      allVarsH[ivar]->Add(mcIncH);
      allVarsH[ivar]->Divide(mcIncH);
      TGraph *mcVarGr=new TGraph(allVarsH[ivar]);
      mcVarGr->SetLineColor(varColors[ivar+1]);
      mcVarGr->SetLineStyle(varStyles[ivar+1]);
      mcVarGr->Draw("l");
      delete allVarsH[ivar];
    }


  c->Modified();
  c->Update();
  c->SaveAs(tag+".png");


  //now compare with QCD templates
  if(firstBin!=1) return;
  TString flav("");
  if(tag.BeginsWith("b"))    flav="b";
  if(tag.BeginsWith("udsg")) flav="udsg";
  if(tag.BeginsWith("c"))    flav="c";
  if(mc.qcdTemplatesH.find(flav)==mc.qcdTemplatesH.end()) return;
  std::vector<TH1F *> &qcd=mc.qcdTemplatesH.find(flav)->second;
  TCanvas *cqcd=new TCanvas(tag+"qcdc",tag+"qcdc",800,800);
  cqcd->Divide(2,2);
  for(size_t i=0; i<qcd.size(); i++)
    {
      cqcd->cd(i+1);
      TString nn(tag+"incmc_"); nn+=i;

      TH1 *mcH=mc.nominalH->ProjectionX(nn,i+1,i+1); mcH->SetDirectory(0);
      TString kinLabel=mc.nominalH->GetYaxis()->GetBinLabel(i+1);
      mcH->Scale(1./mcH->Integral());
      qcd[i]->Draw("hist");
      qcd[i]->GetXaxis()->SetTitle("Discriminator");
      qcd[i]->GetYaxis()->SetTitle("Jets (a.u.)");
      TGraphErrors *qcdGr=new TGraphErrors(qcd[i]);
      qcdGr->SetFillStyle(3002); qcdGr->SetFillColor(qcd[i]->GetLineColor());
      qcdGr->Draw("e2");

      mcH->Draw("hist same");
      TGraphErrors *mcGr=new TGraphErrors(mcH);
      mcGr->SetFillStyle(3002); 
      mcGr->SetFillColor(mcH->GetLineColor());
      mcGr->Draw("e2");

      char buf[100];
      sprintf(buf,"#chi^{2}/ndof : %3.2f", qcd[i]->Chi2Test(mcH,"WWCHI2/NDF") );
      TPaveText *pave = new TPaveText(0.5,0.75,0.93,0.85,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      pave->SetTextAlign(12);
      pave->SetTextFont(42);
      pave->AddText(kinLabel);
      //pave->AddText(buf);
      pave->Draw();

      if(i>0) continue;
      pave = new TPaveText(0.2,0.93,0.8,0.99,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      //      pave->SetTextFont(42);
      pave->AddText("CMS simulation, " + flav);
      pave->Draw();
      
      TLegend *leg=new TLegend(0.5,0.50,0.93,0.75);
      //leg->SetHeader(kinLabel);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->SetTextFont(42);  
      leg->AddEntry(qcd[i],"QCD","L");
      leg->AddEntry(mcH,"t#bar{t} sample","L"); 
      leg->Draw();
    }
  cqcd->Modified();
  cqcd->Update();
  cqcd->SaveAs(tag+"_qcd.png");
}


//
TGraphAsymmErrors *getUncertaintyBandFrom(std::vector<TH1 *> &relVarH)
{

  if(relVarH.size()==0) return 0;
  TH1 *varUpH= (TH1 *)relVarH[0]->Clone("varup");     varUpH->Reset("ICE");
  TH1 *varDownH= (TH1 *)relVarH[0]->Clone("vardown"); varDownH->Reset("ICE");
  for(std::vector<TH1 *>::iterator hit=relVarH.begin(); hit!=relVarH.end(); hit++)
    {
      for(int xbin=1; xbin<=(*hit)->GetXaxis()->GetNbins(); xbin++)
	{
	  Double_t diff=(*hit)->GetBinContent(xbin);
	  TH1 *hToIncrement=(diff<0 ? varUpH : varDownH);
	  Double_t curVal=hToIncrement->GetBinContent(xbin);
	  curVal += pow(diff,2);
	  hToIncrement->SetBinContent(xbin,curVal);
	}
    }
  
  //finalize the graph
  TGraphAsymmErrors *gr=new TGraphAsymmErrors;
  for(int xbin=1; xbin<=varUpH->GetXaxis()->GetNbins(); xbin++)
    {
      Double_t errUp=sqrt(varUpH->GetBinContent(xbin));
      Double_t errDown=sqrt(varDownH->GetBinContent(xbin));
      int np=gr->GetN();
      gr->SetPoint(np,varUpH->GetBinCenter(xbin),0);
      gr->SetPointError(np,0.5*varUpH->GetBinWidth(xbin),0.5*varUpH->GetBinWidth(xbin),errDown,errUp);
    }

  delete varUpH;
  delete varDownH;
  return gr;
}

//
void fitData(TH1 *preTagData, TH1 *postTagData, TH1 *vetoData,
	     std::vector<TH1 *> preTagTemplates, std::vector<TH1 *> postTagTemplates, std::vector<TH1 *> vetoTemplates,
	     TString tag, TString tagLabel)
{
  RooRealVar x("x","Discriminator",preTagData->GetXaxis()->GetXmin(),preTagData->GetXaxis()->GetXmax());

  //create the data set
  RooCategory dataCategories("categ","categ") ;
  dataCategories.defineType("pretag");
  dataCategories.defineType("posttag");
  if(vetoData) dataCategories.defineType("veto");
  std::map<std::string, TH1 *> dataH;
  dataH["pretag"]=preTagData;
  dataH["posttag"]=postTagData;
  if(vetoData) dataH["veto"]=vetoData;
  RooDataHist combData("combData","combData",RooArgSet(x),dataCategories,dataH);

  //create the summed pdf
  RooArgSet preTagFlavorPdfs,   postTagFlavorPdfs,   vetoFlavorPdfs;
  RooArgSet preTagFlavorYields, postTagFlavorYields, vetoFlavorYields;
  RooArgSet tagEffs;
  for(size_t i=0; i<preTagTemplates.size(); i++)
    {
      TH1 *h=preTagTemplates[i];
      RooDataHist *itempl=new RooDataHist(h->GetName(),h->GetName(), RooArgList(x), Import(*h));
      RooHistPdf *ipdf=new RooHistPdf(itempl->GetName()+TString("pdf"),itempl->GetName()+TString("pdf"), RooArgSet(x), *itempl);
      preTagFlavorPdfs.add(*ipdf);
      double totalExpErr;
      double totalExp=h->IntegralAndError(1,h->GetXaxis()->GetNbins(),totalExpErr);
      RooRealVar *iSF = new RooRealVar(h->GetName()+TString("sf"),h->GetName()+TString("sf"),1,0,5);
      RooFormulaVar *ipreN = new RooFormulaVar(h->GetName()+TString("yields"),h->GetName()+TString("yields"),"@0*@1",RooArgList(*iSF,RooConst(totalExp)));
      preTagFlavorYields.add(*ipreN);
      
      h=postTagTemplates[i];
      itempl=new RooDataHist(h->GetName(),h->GetName(), RooArgList(x), Import(*h));
      ipdf=new RooHistPdf(itempl->GetName()+TString("pdf"),itempl->GetName()+TString("pdf"), RooArgSet(x), *itempl);
      postTagFlavorPdfs.add(*ipdf);
      double totalPostExpErr;
      double totalPostExp=h->IntegralAndError(1,h->GetXaxis()->GetNbins(),totalPostExpErr);
      RooRealVar *iEff = new RooRealVar(h->GetName()+TString("sf"),h->GetName()+TString("sf"),1,0,2);
      RooFormulaVar *ipostN = new RooFormulaVar(h->GetName()+TString("yields"),h->GetName()+TString("yields"),"@0*@1*@2",RooArgList(*iEff,RooConst(totalPostExp/totalExp),*ipreN));
      postTagFlavorYields.add(*ipostN);

      if(vetoData==0) continue;
      h=vetoTemplates[i];
      itempl=new RooDataHist(h->GetName(),h->GetName(), RooArgList(x), Import(*h));
      ipdf=new RooHistPdf(itempl->GetName()+TString("pdf"),itempl->GetName()+TString("pdf"), RooArgSet(x), *itempl);
      vetoFlavorPdfs.add(*ipdf);
      RooFormulaVar *ivetoN = new RooFormulaVar(h->GetName()+TString("yields"),h->GetName()+TString("yields"),"@0-@1",RooArgList(*ipreN,*ipostN));
      vetoFlavorYields.add(*ivetoN);
    }

  //create a simultaneous PDF for each category
  RooSimultaneous simPdf("simPdf","simultaneous pdf",dataCategories) ;
  RooAddPdf preTagShapeModel("pretagshapemodel",   "pretag flavor sum",  preTagFlavorPdfs,  preTagFlavorYields);
  simPdf.addPdf(preTagShapeModel,  "pretag");
  RooAddPdf postTagShapeModel("posttagshapemodel", "posttag flavor sum", postTagFlavorPdfs, postTagFlavorYields);
  simPdf.addPdf(postTagShapeModel, "posttag");
  if(vetoData) 
    {
      RooAddPdf *vetoShapeModel=new RooAddPdf("vetoshapemodel", "veto flavor sum", vetoFlavorPdfs, vetoFlavorYields);
      simPdf.addPdf(*vetoShapeModel, "veto");
    }

  //add all to a workspace
  RooWorkspace* w = new RooWorkspace("w") ;
  w->import(RooArgSet(x,dataCategories));
  w->import(combData);
  w->import(simPdf);
  w->Print("v");

  //fit it profiling the b-tagging efficiency
  RooFitResult *fitres = w->pdf("simPdf")->	fitTo(combData,Save(kTRUE),RooFit::Minos(),Extended(kTRUE));
  //RooRealVar *sf_effb=w->var("bposttagsf");
  // RooAbsReal *nll = w->pdf("simPdf")->createNLL(combData);
  // RooMinuit (*nll).migrad();
  // RooAbsReal* pll_sf_effb = nll->createProfile(*sf_effb);

  //display the result
  TCanvas *c = new TCanvas(tag,tag,1000,1000); 
  c->Divide(2,2);
  const string keys[]={"pretag","posttag","veto"};
  const string flavs[]={"b","udsg","c"};
  for(int i=0; i<3; i++)
    {
      std::map<string,TH1 *>::iterator it=dataH.find(keys[i]);
      if(it==dataH.end()) continue;
      
      c->cd(i+1);
      RooPlot *frame = x.frame();
    
      //project data
      RooDataHist* dataslice = (RooDataHist *)combData.reduce(("categ==categ::"+it->first).c_str());
      dataslice->plotOn(frame,DataError(RooAbsData::SumW2),Name((it->first+"data").c_str()));
      
      //project PDF
      RooCategory theCategory(it->first.c_str(),it->first.c_str());
      w->pdf("simPdf")->plotOn(frame,Slice(theCategory),ProjWData(x,*dataslice));

      //project b and c components
      w->pdf("simPdf")->plotOn(frame,
			       RooFit::Components( * w->pdf( ("b"+it->first+"pdf").c_str()  ) ), 
			       ProjWData(x,*dataslice),
			       RooFit::FillStyle(1001), RooFit::FillColor(kOrange-2),
			       MoveToBack(),            DrawOption("lf")  );
      w->pdf("simPdf")->plotOn(frame,
			       RooFit::Components( RooArgSet(* w->pdf( ("c"+it->first+"pdf").c_str()  ),
							     * w->pdf( ("b"+it->first+"pdf").c_str()  ) )
						   ), 
			       ProjWData(x,*dataslice),
			       RooFit::FillStyle(1001), RooFit::FillColor(kAzure+9),
			       MoveToBack(),            DrawOption("lf")  );
      frame->SetTitleOffset(1.3,"Y");
      frame->SetYTitle("Jets");
      frame->Draw();

      TPaveText* T = new TPaveText(0.5,0.65,0.88,0.88, "NDC");
      T->SetFillColor(0);
      T->SetFillStyle(0);
      T->SetTextAlign(12);
      T->SetLineColor(0);
      T->SetBorderSize(0);
      T->SetTextFont(42);
      std::string itag=it->first;
      if(itag=="pretag")      itag="Pre-tagged";
      else if(itag=="posttag") itag="Tagged";
      else               itag="Vetoed";
      T->AddText(itag.c_str());
      if(itag!="Vetoed")
	{
	  for(size_t iflav=0; iflav<3; iflav++)
	    {
	      RooRealVar *v=w->var((flavs[iflav]+it->first+"sf").c_str());
	      if(v==0) continue;
	      char buf[100]; 
	      TString name("SF_{");  name += flavs[iflav].c_str(); name += "}";
	      if(itag=="Pre-tagged") name += "^{flavor}";
	      sprintf(buf,"%s = %3.3f #pm %3.3f",name.Data(),v->getVal(),v->getError());
	      T->AddText(buf);
	    }
	}
      else T->SetY1(0.8);
      T->Draw("same");

      if(i>0) continue;
      T = new TPaveText(0.1,0.92,0.9,0.98, "NDC");
      T->SetFillColor(0);
      T->SetFillStyle(0);  T->SetLineColor(0);
      T->SetTextAlign(22);
      char Buffer[1024]; sprintf(Buffer, "CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int} L = 5.0 fb^{-1}, %s", tagLabel.Data());
      T->AddText(Buffer);
      T->Draw("same");
    }

  //draw the likelihood as function of epsilon_b
  c->cd(4);
  /*  
  float xcen=sf_effb->getVal();
  float xerr=5*sf_effb->getError();

  RooPlot* frame = sf_effb->frame(Range(max(float(xcen-xerr),float(0.)),xcen+xerr)); //min(float(xcen+xerr),float(1.0) ) ) );
  nll->plotOn(frame,ShiftToZero(),LineColor(kGray),LineStyle(9),Name("ll")) ;  
  pll_sf_effb->plotOn(frame,Name("pll"),ShiftToZero()) ;
  frame->Draw();
  frame->SetMinimum(0); frame->SetMaximum(3);
  frame->GetYaxis()->SetTitle();
  frame->GetXaxis()->SetTitle("SF_{b}=#varepsilon_{b}/#varepsilon_{b}^{MC}");

  TLegend *leg=new TLegend(0.5,0.65,0.88,0.88);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextAlign(12);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry("ll","-log(L/L_{max})","l");
  leg->AddEntry("pll","PLR","l");
  leg->Draw();
  */
  //update and save
  c->Modified();
  c->Update();
}
