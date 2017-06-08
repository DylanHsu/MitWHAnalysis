void finalPlot2017(
  TString input_file="MitWHAnalysis/plots/test_newplots_noprescale/histo_wh_nice.root",
  TString plotName="presel mT",
  TString jetString="nj",
  TString lepString="l",
  TString xlabel="m_{T}",
  TString units="GeV",
  TString output = "histo_nice",
  int rebin=5,
  bool logScale=true,
  double xmin=0,
  double xmax=0,
  double ymin=0,
  double ymax=0,
  bool ratioPad=true,
  bool isBlind=true
) {
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open(input_file,"READ");
  TH1D *histo_data  = (TH1D*) f->Get(Form("%s/%s/%s/Data"        , plotName.Data(), jetString.Data(), lepString.Data()));
  TH1D *histo_top   = (TH1D*) f->Get(Form("%s/%s/%s/Top"         , plotName.Data(), jetString.Data(), lepString.Data()));
  TH1D *histo_EWKW  = (TH1D*) f->Get(Form("%s/%s/%s/EWK W+jets"  , plotName.Data(), jetString.Data(), lepString.Data()));
  TH1D *histo_Wjets = (TH1D*) f->Get(Form("%s/%s/%s/Cont. W+jets", plotName.Data(), jetString.Data(), lepString.Data()));
  TH1D *histo_Zjets = (TH1D*) f->Get(Form("%s/%s/%s/Z+jets", plotName.Data(), jetString.Data(), lepString.Data()));
  TH1D *histo_WZ    = (TH1D*) f->Get(Form("%s/%s/%s/WZ"          , plotName.Data(), jetString.Data(), lepString.Data()));
  TH1D *histo_ZZ    = (TH1D*) f->Get(Form("%s/%s/%s/ZZ"          , plotName.Data(), jetString.Data(), lepString.Data()));
  TH1D *histo_WW    = (TH1D*) f->Get(Form("%s/%s/%s/WW"          , plotName.Data(), jetString.Data(), lepString.Data()));
  TH1D *histo_qcd   = (TH1D*) f->Get(Form("%s/%s/%s/QCD, y+jets"  , plotName.Data(), jetString.Data(), lepString.Data()));
  TH1D *histo_WHinv = (TH1D*) f->Get(Form("%s/%s/%s/WH(125)"     , plotName.Data(), jetString.Data(), lepString.Data()));
  histo_data  ->Rebin(rebin);printf("hi histo_data  \n");
  histo_top   ->Rebin(rebin);printf("hi histo_top   \n");
  histo_EWKW  ->Rebin(rebin);printf("hi histo_EWKW  \n");
  histo_Wjets ->Rebin(rebin);printf("hi histo_Wjets \n");
  histo_Zjets ->Rebin(rebin);printf("hi histo_Zjets \n");
  histo_WZ    ->Rebin(rebin);printf("hi histo_WZ    \n");
  histo_ZZ    ->Rebin(rebin);printf("hi histo_ZZ    \n");
  histo_WW    ->Rebin(rebin);printf("hi histo_WW    \n");
  histo_qcd   ->Rebin(rebin);printf("hi histo_qcd   \n");
  histo_WHinv ->Rebin(rebin);printf("hi histo_WHinv \n");
  THStack *hs = new THStack("hs", "");
  histo_qcd   ->SetFillColor(kGray     );   histo_qcd   ->SetLineColor( 1 );   
  histo_top   ->SetFillColor(kPink+7   );   histo_top   ->SetLineColor( 1 );   
  histo_EWKW  ->SetFillColor(kViolet-9 );   histo_EWKW  ->SetLineColor( 1 );  
  histo_Wjets ->SetFillColor(kViolet+8 );   histo_Wjets ->SetLineColor( 1 ); 
  histo_Zjets ->SetFillColor(901       );   histo_Zjets ->SetLineColor( 1 ); 
  histo_WZ    ->SetFillColor(856       );   histo_WZ    ->SetLineColor( 1 );     
  histo_ZZ    ->SetFillColor(842       );   histo_ZZ    ->SetLineColor( 1 );     
  histo_WW    ->SetFillColor(kAzure-9  );   histo_WW    ->SetLineColor( 1 );     
  histo_data  ->SetMarkerStyle(20      );   histo_data  ->SetLineColor( 1 ); 
  histo_WHinv ->SetLineWidth(4         );   histo_WHinv ->SetLineColor( kRed+3);     
  hs->Add( histo_EWKW  ); 
  hs->Add( histo_ZZ    ); 
  hs->Add( histo_Zjets ); 
  hs->Add( histo_WW ); 
  hs->Add( histo_WZ    ); 
  hs->Add( histo_Wjets ); 
  hs->Add( histo_top   );
  hs->Add( histo_qcd   );
  TString canvas_title=Form("c_%s_%s_%s", (plotName.ReplaceAll(" ","_")).Data(), jetString.Data(), lepString.Data());
  TCanvas *canvas = new TCanvas(canvas_title, canvas_title, 800, 800);
  canvas->SetBottomMargin(0.1);
  canvas->SetRightMargin(0.05);
  canvas->SetLeftMargin(0.2);
  TPad *pad1, *pad2;
  if(ratioPad) {
    pad1 = new TPad("pad1", "pad1",0.00,0.30,1.00,1.00);
    pad2 = new TPad("pad2", "pad2",0.00,0.00,1.00,0.30);
    pad1->SetBottomMargin(0.01);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.3);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    if(logScale) pad1->SetLogy();
  } else {
    if(logScale) canvas->SetLogy();
  }
  hs->Draw("HIST");
  if(xmin!=0 || xmax!=0 )  hs->GetXaxis()->SetRangeUser(xmin, xmax);
  if((ymin!=0 || ymax!=0) && ymin<ymax && (!logScale || ymin>0) ) { 
    hs->SetMaximum(ymax);
    hs->SetMinimum(ymin);
  }
  histo_WHinv->Draw("HIST SAME");
  if(!isBlind) histo_data->Draw("E X0 P SAME");
  TString title;
  if(units!="") title=Form("%s [%s]", xlabel.Data(), units.Data());
  else          title=xlabel;
  if(units!="") hs->GetYaxis()->SetTitle(Form("Events / %d %s", int(histo_Wjets->GetBinWidth(1)), units.Data()));
  else          hs->GetYaxis()->SetTitle(Form("Events / %d", int(histo_Wjets->GetBinWidth(1))));
  
  if(ratioPad) {
    hs->GetYaxis()->SetLabelSize(0.045);
    hs->GetYaxis()->SetTitleSize(0.05);
    hs->GetYaxis()->SetTitleOffset(0.8);
  } else {
    hs->GetXaxis()->SetTitle(title);
    hs->GetYaxis()->SetTitleOffset(2.1);
  }
  
  TLegend *l1, *l2;
  if(ratioPad) { l1 = new TLegend(0.57,0.4,0.75,0.88); l2 = new TLegend(0.77,0.51,0.95,0.88); }
  else         { l1 = new TLegend(0.57,0.6,0.75,0.88); l2 = new TLegend(0.77,0.655,0.95,0.88); }
  l1->SetFillColor(0);
  l1->SetBorderSize(0);
  l2->SetFillColor(0);
  l2->SetBorderSize(0);
  l1->AddEntry(histo_qcd   , "QCD/#gamma+jets"   , "f");
  l1->AddEntry(histo_top   , "Top"               , "f");
  l1->AddEntry(histo_Wjets , "Cont. W+jets"      , "f");
  l1->AddEntry(histo_WZ    , "WZ"                , "f");  
  l1->AddEntry(histo_WW    , "WW"                , "f");  
  l2->AddEntry(histo_Zjets , "Z+jets"            , "f");
  l2->AddEntry(histo_ZZ    , "ZZ"                , "f");  
  l2->AddEntry(histo_EWKW  , "EWK W+jets"        , "f");
  l2->AddEntry(histo_WHinv , "WH(125)"           , "l");
  l1->Draw("SAME");
  l2->Draw("SAME");
  
  if(ratioPad) {
    pad2->cd(); 
    TH1D *histo_sum   = (TH1D*) histo_top->Clone("histo_sum");
    histo_sum->Add(histo_EWKW ); 
    histo_sum->Add(histo_Wjets); 
    histo_sum->Add(histo_Zjets); 
    histo_sum->Add(histo_WZ   ); 
    histo_sum->Add(histo_ZZ   ); 
    histo_sum->Add(histo_WW   ); 
    histo_sum->Add(histo_qcd  ); 
    
    TH1D *histo_ratio = (TH1D*) histo_data->Clone("histo_ratio");
    histo_ratio->Divide(histo_sum);
    if(xmin!=0 || xmax!=0 )  histo_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
    histo_ratio->SetTitle("");
    histo_ratio->Draw("E X0 P");
    histo_ratio->GetXaxis()->SetTitle(title);
    histo_ratio->GetXaxis()->SetLabelSize(0.1);
    histo_ratio->GetXaxis()->SetTitleSize(0.12);
    histo_ratio->GetYaxis()->SetTitle("Data/MC");
    histo_ratio->GetYaxis()->SetRangeUser(0.2,1.8);
    histo_ratio->GetYaxis()->SetTitleOffset(0.3);
    histo_ratio->GetYaxis()->SetTitleSize(0.115);
    histo_ratio->GetYaxis()->SetLabelSize(0.09);
    histo_ratio->GetYaxis()->SetNdivisions(4);
    histo_ratio->GetYaxis()->CenterTitle();
    TLine* baseline;
    if(xmin!=0 || xmax!=0 )  baseline=new TLine(xmin,1.,xmax,1.);
    else                     baseline=new TLine(histo_ratio->GetXaxis()->GetXmin(), 1., histo_ratio->GetXaxis()->GetXmax(), 1.); 
    baseline->SetLineStyle(kDashed);
    baseline->Draw();
 
    pad2->Update();
  }
  canvas->Print(Form("%s.png",output.Data()));
}
