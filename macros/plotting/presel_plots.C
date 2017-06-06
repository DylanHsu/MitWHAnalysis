void presel_plots() {
  gStyle->SetOptStat(0);
  TFile *f = new TFile("MitWHAnalysis/plots/histo_wh_nice.root", "READ");
  TH1D *histo_mT_top_presel   = (TH1D*) f->Get("histo_mT_presel_nj_l_1");
  TH1D *histo_mT_EWKW_presel  = (TH1D*) f->Get("histo_mT_presel_nj_l_2");
  TH1D *histo_mT_Wjets_presel = (TH1D*) f->Get("histo_mT_presel_nj_l_3");
  TH1D *histo_mT_Zjets_presel = (TH1D*) f->Get("histo_mT_presel_nj_l_4");
  TH1D *histo_mT_WZ_presel    = (TH1D*) f->Get("histo_mT_presel_nj_l_5");
  TH1D *histo_mT_ZZ_presel    = (TH1D*) f->Get("histo_mT_presel_nj_l_6");
  histo_mT_top_presel   ->Rebin(5);
  histo_mT_EWKW_presel  ->Rebin(5);
  histo_mT_Wjets_presel ->Rebin(5);
  histo_mT_Zjets_presel ->Rebin(5);
  histo_mT_WZ_presel    ->Rebin(5);
  histo_mT_ZZ_presel    ->Rebin(5);
  THStack *hs_mT_presel = new THStack("hs_mT_presel", "Transverse mass at preselection level");
  histo_mT_top_presel   ->SetFillColor(kPink+7   );   histo_mT_top_presel   ->SetLineColor(kPink+7   );   
  histo_mT_EWKW_presel  ->SetFillColor(kViolet-9 );   histo_mT_EWKW_presel  ->SetLineColor(kViolet-9 );  
  histo_mT_Wjets_presel ->SetFillColor(kViolet+8 );   histo_mT_Wjets_presel ->SetLineColor(kViolet+8 ); 
  histo_mT_Zjets_presel ->SetFillColor(901       );   histo_mT_Zjets_presel ->SetLineColor(901       ); 
  histo_mT_WZ_presel    ->SetFillColor(856       );   histo_mT_WZ_presel    ->SetLineColor(856       );     
  histo_mT_ZZ_presel    ->SetFillColor(842       );   histo_mT_ZZ_presel    ->SetLineColor(842       );     
  hs_mT_presel->Add( histo_mT_ZZ_presel    ); 
  hs_mT_presel->Add( histo_mT_EWKW_presel  ); 
  hs_mT_presel->Add( histo_mT_Zjets_presel ); 
  hs_mT_presel->Add( histo_mT_WZ_presel    ); 
  hs_mT_presel->Add( histo_mT_Wjets_presel ); 
  hs_mT_presel->Add( histo_mT_top_presel   ); 
  TCanvas *c_mT_presel = new TCanvas("c_mT_presel", "c_mT_presel");
  c_mT_presel->SetLogy();
  hs_mT_presel->Draw("HIST");
  hs_mT_presel->GetXaxis()->SetRangeUser(125.,600.);
  hs_mT_presel->GetYaxis()->SetRangeUser(1.,400000.);
  hs_mT_presel->GetXaxis()->SetTitle("m_{T} [GeV]");
  TLegend *l_mT_presel = new TLegend(0.5,0.5,0.88,0.88);
  l_mT_presel->SetFillColor(0);
  l_mT_presel->SetBorderSize(0);
  l_mT_presel->AddEntry(histo_mT_top_presel   , "Top"               , "f");
  l_mT_presel->AddEntry(histo_mT_EWKW_presel  , "EWK W+jets"        , "f");
  l_mT_presel->AddEntry(histo_mT_Wjets_presel , "Cont. W+jets"      , "f");
  l_mT_presel->AddEntry(histo_mT_Zjets_presel , "Cont. Z+jets"      , "f");
  l_mT_presel->AddEntry(histo_mT_WZ_presel    , "WZ"                , "f");  
  l_mT_presel->AddEntry(histo_mT_ZZ_presel    , "ZZ"                , "f");  
  l_mT_presel->Draw("SAME");


}
