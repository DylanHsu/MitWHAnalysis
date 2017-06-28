void betterCorrelationPlot(TString theClass="Signal", TString extra_string="") {
  gStyle->SetPalette(kDeepSea);
  TCanvas *canvas = (TCanvas*) gROOT->FindObject(Form("//root/Canvases/CorrelationMatrix%s",theClass.Data()));
  canvas->cd();
  TH2F *matrix = (TH2F*)canvas->GetPrimitive(Form("CorrelationMatrix%s",theClass.Data()));
  matrix->LabelsOption("v","x");
  matrix->GetXaxis()->SetLabelSize(0.028);
  matrix->GetYaxis()->SetLabelSize(0.028);
  canvas->SetWindowSize(1024,1024);
  matrix->SetMarkerSize(0.7);
  canvas->SetLeftMargin(0.25);
  canvas->SetBottomMargin(0.25);
  matrix->Draw("TEXT COLZ");
  canvas->Print(Form("plots/CorrelationMatrix%s",theClass.Data()) + (extra_string == "" ? "" : "_"+extra_string) + ".pdf");
  canvas->Print(Form("plots/CorrelationMatrix%s",theClass.Data()) + (extra_string == "" ? "" : "_"+extra_string) + ".png");


}
