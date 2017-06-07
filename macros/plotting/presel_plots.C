void presel_plots(
  TString input_file="MitWHAnalysis/plots/njets_test/histo_wh_nice.root",
  TString selection="presel",
  TString var="mT",
  TString jetstr="nj",
  TString xlabel="m_{T}",
  TString units="GeV",
  TString output = "histo_nice",
  int rebin=5,
  bool logScale=true,
  double xmin=0,
  double xmax=0,
  double ymin=0,
  double ymax=0
) {
  gStyle->SetOptStat(0);
  TFile *f = new TFile(input_file, "READ");
  TH1D *histo_top   = (TH1D*) f->Get(Form("histo_%s_%s_%s_l_1", var.Data(), selection.Data(), jetstr.Data()));
  TH1D *histo_EWKW  = (TH1D*) f->Get(Form("histo_%s_%s_%s_l_2", var.Data(), selection.Data(), jetstr.Data()));
  TH1D *histo_Wjets = (TH1D*) f->Get(Form("histo_%s_%s_%s_l_3", var.Data(), selection.Data(), jetstr.Data()));
  TH1D *histo_Zjets = (TH1D*) f->Get(Form("histo_%s_%s_%s_l_4", var.Data(), selection.Data(), jetstr.Data()));
  TH1D *histo_WZ    = (TH1D*) f->Get(Form("histo_%s_%s_%s_l_5", var.Data(), selection.Data(), jetstr.Data()));
  TH1D *histo_ZZ    = (TH1D*) f->Get(Form("histo_%s_%s_%s_l_6", var.Data(), selection.Data(), jetstr.Data()));
  TH1D *histo_WW    = (TH1D*) f->Get(Form("histo_%s_%s_%s_l_7", var.Data(), selection.Data(), jetstr.Data()));
  histo_top   ->Rebin(rebin);
  histo_EWKW  ->Rebin(rebin);
  histo_Wjets ->Rebin(rebin);
  histo_Zjets ->Rebin(rebin);
  histo_WZ    ->Rebin(rebin);
  histo_ZZ    ->Rebin(rebin);
  histo_WW    ->Rebin(rebin);
  THStack *hs = new THStack("hs", "");
  histo_top   ->SetFillColor(kPink+7   );   histo_top   ->SetLineColor( 1 );   
  histo_EWKW  ->SetFillColor(kViolet-9 );   histo_EWKW  ->SetLineColor( 1 );  
  histo_Wjets ->SetFillColor(kViolet+8 );   histo_Wjets ->SetLineColor( 1 ); 
  histo_Zjets ->SetFillColor(901       );   histo_Zjets ->SetLineColor( 1 ); 
  histo_WZ    ->SetFillColor(856       );   histo_WZ    ->SetLineColor( 1 );     
  histo_ZZ    ->SetFillColor(842       );   histo_ZZ    ->SetLineColor( 1 );     
  histo_WW    ->SetFillColor(kAzure-9  );   histo_WW    ->SetLineColor( 1 );     
  //histo_EWKW->SetMaximum(400000);
  hs->Add( histo_EWKW  ); 
  hs->Add( histo_ZZ    ); 
  hs->Add( histo_Zjets ); 
  hs->Add( histo_WW ); 
  hs->Add( histo_WZ    ); 
  hs->Add( histo_Wjets ); 
  hs->Add( histo_top   );
  TString canvas_title=Form("c_%s_%s", var.Data(), jetstr.Data());
  TCanvas *canvas = new TCanvas(canvas_title, canvas_title);
  if(logScale) canvas->SetLogy();
  hs->Draw("HIST");
  if(xmin!=0 || xmax!=0 )  hs->GetXaxis()->SetRangeUser(xmin, xmax);
  if((ymin!=0 || ymax!=0) && ymin<ymax && (!logScale || ymin>0) ) { 
    hs->SetMaximum(ymax);
    hs->SetMinimum(ymin);
  }
  TString title;
  if(units!="") title=Form("%s [%s]", xlabel.Data(), units.Data());
  else          title=xlabel;
  hs->GetXaxis()->SetTitle(title);
  TLegend *l = new TLegend(0.6,0.6,0.88,0.88);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(histo_top   , "Top"               , "f");
  l->AddEntry(histo_Wjets , "Cont. W+jets"      , "f");
  l->AddEntry(histo_WZ    , "WZ"                , "f");  
  l->AddEntry(histo_WW    , "WW"                , "f");  
  l->AddEntry(histo_Zjets , "Cont. Z+jets"      , "f");
  l->AddEntry(histo_ZZ    , "ZZ"                , "f");  
  l->AddEntry(histo_EWKW  , "EWK W+jets"        , "f");
  l->Draw("SAME");
  canvas->Print(Form("%s.png",output.Data()));
}
