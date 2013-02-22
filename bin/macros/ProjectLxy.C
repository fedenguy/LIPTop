{
  TFile *inF=TFile::Open("~/work/top/plotter.root");
  //TH2F *h=(TH2F *)inF->Get("t#bar{t}/jetlxyreseta");
  TH2F *h=(TH2F *)inF->Get("t#bar{t}/jetlxyrespt");
  h->SetDirectory(0);
  inF->Close();

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *c=new TCanvas("c","c",600,600);
  c->cd();
  for(int ybin=1; ybin<=h->GetYaxis()->GetNbins(); ybin++)
    {
      char buf[200];
      sprintf(buf,"%3.1f-%3.1f",h->GetYaxis()->GetBinLowEdge(ybin),h->GetYaxis()->GetBinLowEdge(ybin+1));
      TString name("proj"); name+= ybin;
      TH1D *hproj=h->ProjectionX(name,ybin,ybin);
      hproj->SetMarkerStyle(ybin+20);
      hproj->SetTitle(buf);
      hproj->Draw(ybin==1?"":"same");
      hproj->SetFillStyle(0);
    }
  c->BuildLegend()->SetBorderSize(0);
}
