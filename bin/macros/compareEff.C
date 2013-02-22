{
  TString ch("ee");

  TFile *inF=TFile::Open("/afs/cern.ch/user/p/psilva/work/top/plotter.root");
  TH1F *ttbarMet=(TH1F *)inF->Get("t#bar{t}/"+ch+"_met");
  ttbarMet->SetDirectory(0);
  TH1F *ttbarMvaMet=(TH1F *)inF->Get("t#bar{t}/"+ch+"_mvamet");
  ttbarMvaMet->SetDirectory(0);
  TH1F *dyMet=(TH1F *)inF->Get("Z#rightarrow ll/"+ch+"_met");
  dyMet->SetDirectory(0);
  TH1F *dyMvaMet=(TH1F *)inF->Get("Z#rightarrow ll/"+ch+"_mvamet");
  dyMvaMet->SetDirectory(0);

  TGraphErrors *metEff=new TGraphErrors;
  TGraphErrors *mvaMetEff=new TGraphErrors;
  Double_t ttbarTotal=ttbarMet->Integral(0,ttbarMet->GetXaxis()->GetNbins()+1);
  Double_t dyTotal   =dyMet->Integral(0,ttbarMet->GetXaxis()->GetNbins()+1);
  for(int ibin=1; ibin<=ttbarMet->GetXaxis()->GetNbins(); ibin++)
    {
      Double_t ttbarierr,dyierr;
      Double_t ttbari=ttbarMet->IntegralAndError(ibin,ttbarMet->GetXaxis()->GetNbins()+1,ttbarierr);
      Double_t dyi   =dyMet->IntegralAndError(ibin,ttbarMet->GetXaxis()->GetNbins()+1,dyierr);
      Int_t np=metEff->GetN();
      metEff->SetPoint(np,ttbari/ttbarTotal,dyi/dyTotal);
      metEff->SetPoint(np,ttbarierr/ttbarTotal,dyierr/dyTotal);

      ttbari=ttbarMvaMet->IntegralAndError(ibin,ttbarMet->GetXaxis()->GetNbins()+1,ttbarierr);
      dyi   =dyMvaMet->IntegralAndError(ibin,ttbarMet->GetXaxis()->GetNbins()+1,dyierr);
      mvaMetEff->SetPoint(np,ttbari/ttbarTotal,dyi/dyTotal);
      mvaMetEff->SetPoint(np,ttbarierr/ttbarTotal,dyierr/dyTotal);
    }

  inF->Close();
PiS818

  TCanvas *c=new TCanvas;
  metEff->Draw("ap"); 
  metEff->SetFillStyle(0); 
  metEff->SetTitle("PF E_{T}^{miss}"); 
  metEff->GetXaxis()->SetTitle("t#{bar} efficiency");
  metEff->GetXaxis()->SetTitle("DY efficiency");
  metEff->SetMarkerStyle(20);
  metEff->SetFillStyle(0);

  mvaMetEff->Draw("p"); 
  mvaMetEff->SetMarkerStyle(24);
  mvaMetEff->SetFillStyle(0);

  c->BuildLegend()->SetTextFont(42);

}
