{
  TString files[]={"MC8TeV_TTJets_filt1.root",
		   "MC8TeV_TTJets_166v5_filt1.root",
		   "MC8TeV_TTJets_178v5_filt1.root",
		   "MC8TeV_TTJets_matchingdown.root",
		   "MC8TeV_TTJets_matchingup.root",
		   "MC8TeV_TTJets_scaleup.root",
		   "MC8TeV_TTJets_scaledown.root"
  };
  TString titles[]={"172.5 GeV","166.5 GeV","178.5 GeV", "172.5 GeV ME-PS up", "172.5 GeV ME-PS down", "172.5 GeV (2Q)^{2}", "172.5 (Q/2)^{2}"};
  TObjArray histos;

  for(int i=0; i<sizeof(files)/sizeof(TString); i++)
    {
      TFile *fIn=TFile::Open("~/work/top/2012/"+files[i]);
      h=(TH1 *) fIn->Get("all_jetlxy");
      h->SetTitle(titles[i]);
      h->SetName(titles[i]);
      if(i==0) h->SetLineWidth(2);
      h->SetLineColor(i%3+1);
      h->SetMarkerColor(i%3+1);
      h->SetMarkerStyle(1);
      h->SetDirectory(0);
      h->Scale(1./h->Integral());
      histos.Add(h);
      fIn->Close();
    }

  gStyle->SetOptStat(1111);
  gStyle->SetOptTitle(0);

  TCanvas *c= new TCanvas("c1","c1",600,600);
  c->cd();
  histos.At(0)->Draw("hist");
  histos.At(1)->Draw("histsames");
  histos.At(2)->Draw("histsames");
  TLegend *leg=c->BuildLegend();
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);

  c= new TCanvas("c2","c2",600,600);
  c->cd();
  histos.At(0)->Draw("hist");
  histos.At(3)->Draw("histsames");
  histos.At(4)->Draw("histsames");
  TLegend *leg=c->BuildLegend();
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);

  c= new TCanvas("c3","c3",600,600);
  c->cd();
  histos.At(0)->Draw("hist");
  histos.At(5)->Draw("histsames");
  histos.At(6)->Draw("histsames");
  TLegend *leg=c->BuildLegend();
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
}
