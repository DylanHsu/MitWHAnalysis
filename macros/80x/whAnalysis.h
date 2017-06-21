#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <Compression.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TSystem.h>
#include <TThread.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
//const int MVAVarType = 1; const int nBinMVA = 8; Double_t xbins[nBinMVA+1] = {100, 125, 150, 175, 200, 250, 300, 400, 500};
//const int MVAVarType = 2; const int nBinMVA = 10; Double_t xbins[nBinMVA+1] = {-1, 0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1};
//const int MVAVarType = 2; const int nBinMVA = 9; Double_t xbins[nBinMVA+1] = {0, 0.2, 0.4, 0.5, 0.54, 0.58, 0.62, 0.66, 0.7, 0.74};
const int MVAVarType = 2; const int nBinMVA = 9; Double_t xbins[nBinMVA+1] = {0.2, 0.4, 0.5, 0.6, 0.65, 0.70, 0.75, 0.8, 0.85, 0.9};
const unsigned int jet_cats=4; // 0-jet, 1-jet, 2-jet, and inclusive
const unsigned int process_types=10; 
const unsigned int multiclassSignal=0;
const unsigned int allPlots=58;
TH1D* makeHisto(unsigned int thePlot, TString &plotName) {
  TH1D *theHisto;
  UInt_t nBinPlot; Float_t xminPlot, xmaxPlot;
  if     (thePlot ==  0) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="presel lep pT"             ;}
  else if(thePlot ==  1) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="presel mT"                 ;}
  else if(thePlot ==  2) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="presel MET"                ;}
  else if(thePlot ==  3) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel lep eta"            ;}
  else if(thePlot ==  4) {nBinPlot = 185; xminPlot = 30.; xmaxPlot = 400.0;       plotName="presel leading jet pT"     ;}
  else if(thePlot ==  5) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;       plotName="presel max CSV2"           ;}
  else if(thePlot ==  6) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =   4.0;       plotName="presel pf calo balance"    ;}
  else if(thePlot ==  7) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel lepton MET dPhi"    ;}
  else if(thePlot ==  8) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel jet MET dPhi"       ;}
  else if(thePlot ==  9) {nBinPlot =2000; xminPlot = 0.0; xmaxPlot =  10.0;       plotName="presel jet lep dR"         ;}
  else if(thePlot == 10) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   100;       plotName="presel NPV "               ;}
  else if(thePlot == 11) {nBinPlot =1000; xminPlot =-500; xmaxPlot =   500;       plotName="presel Uperp"              ;}
  else if(thePlot == 12) {nBinPlot =1000; xminPlot =-500; xmaxPlot =   500;       plotName="presel Upara"              ;}
  else if(thePlot == 13) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel nAwayPFCH"          ;}
  else if(thePlot == 14) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel nTransversePFCH"    ;}
  else if(thePlot == 15) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel nTowardPFCH"        ;}
  else if(thePlot == 16) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel sumPtAwayPFCH"      ;}
  else if(thePlot == 17) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel sumPtTransversePFCH";}
  else if(thePlot == 18) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel sumPtTowardPFCH"    ;}
  else if(thePlot == 19) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="presel PFCH1Pt"            ;}
  else if(thePlot == 20) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="presel PFCH1Eta"           ;}
  else if(thePlot == 21) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel PFCH1Phi"           ;}
  else if(thePlot == 22) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="presel PFCH2Pt"            ;}
  else if(thePlot == 23) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="presel PFCH2Eta"           ;}
  else if(thePlot == 24) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel PFCH2Phi"           ;}
  else if(thePlot == 25) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="presel PFCH3Pt"            ;}
  else if(thePlot == 26) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="presel PFCH3Eta"           ;}
  else if(thePlot == 27) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel PFCH3Phi"           ;}
  else if(thePlot == 28) {nBinPlot =2000; xminPlot = xbins[0]; xmaxPlot = xbins[nBinMVA-1]; plotName="presel MVAVar"   ;}
  else if(thePlot == 29) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="sigsel lep pT"             ;}
  else if(thePlot == 30) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="sigsel mT"                 ;}
  else if(thePlot == 31) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="sigsel lep eta"            ;}
  else if(thePlot == 32) {nBinPlot = 185; xminPlot = 30.; xmaxPlot = 400.0;       plotName="sigsel leading jet pT"     ;}
  else if(thePlot == 33) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =   4.0;       plotName="sigsel pf calo balance"    ;}
  else if(thePlot == 34) {nBinPlot =2000; xminPlot = 0.0; xmaxPlot =  10.0;       plotName="sigsel jet lep dR"         ;}
  else if(thePlot == 35) {nBinPlot =1000; xminPlot =-500; xmaxPlot =   500;       plotName="sigsel Uperp"              ;}
  else if(thePlot == 36) {nBinPlot =1000; xminPlot =-500; xmaxPlot =   500;       plotName="sigsel Upara"              ;}
  else if(thePlot == 37) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel nAwayPFCH"          ;}
  else if(thePlot == 38) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel nTransversePFCH"    ;}
  else if(thePlot == 39) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel nTowardPFCH"        ;}
  else if(thePlot == 40) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel sumPtAwayPFCH"      ;}
  else if(thePlot == 41) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel sumPtTransversePFCH";}
  else if(thePlot == 42) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel sumPtTowardPFCH"    ;}
  else if(thePlot == 43) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="sigsel PFCH1Pt"            ;}
  else if(thePlot == 44) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="sigsel PFCH1Eta"           ;}
  else if(thePlot == 45) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="sigsel PFCH1Phi"           ;}
  else if(thePlot == 46) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="sigsel PFCH2Pt"            ;}
  else if(thePlot == 47) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="sigsel PFCH2Eta"           ;}
  else if(thePlot == 48) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="sigsel PFCH2Phi"           ;}
  else if(thePlot == 49) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="sigsel PFCH3Pt"            ;}
  else if(thePlot == 50) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="sigsel PFCH3Eta"           ;}
  else if(thePlot == 51) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="sigsel PFCH3Phi"           ;}
  else if(thePlot == 52) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =   4.0;       plotName="N-1 lepton MET balance"    ;}
  else if(thePlot == 53) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;       plotName="N-1 max CSV2"              ;}
  else if(thePlot == 54) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="N-1 MET"                   ;}
  else if(thePlot == 55) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="N-1 lepton MET dPhi"       ;}
  else if(thePlot == 56) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="N-1 jet MET dPhi"          ;}
  else if(thePlot == allPlots-1) plotName="Shape analysis";
  else { printf("error with makeHisto: thePlot out of bounds (%d). check allPlots\n", thePlot); return theHisto;}
  
  if(thePlot != allPlots-1) theHisto = new TH1D("theHisto", "theHisto", nBinPlot, xminPlot, xmaxPlot);
  else                      theHisto = new TH1D("theHisto", "theHisto", nBinMVA, xbins);
  theHisto->Sumw2();
  return theHisto;
}
void categoryNames(vector<TString> &categoryName_) {
  categoryName_.push_back("Data");
  categoryName_.push_back("Top");
  categoryName_.push_back("EWK W+jets");
  categoryName_.push_back("Cont. W+jets");
  categoryName_.push_back("Z+jets");
  categoryName_.push_back("WZ");
  categoryName_.push_back("ZZ");
  categoryName_.push_back("WW");
  categoryName_.push_back("QCD, y+jets");
  categoryName_.push_back("WH(125)");
  return;
}
void jetStrings(vector<TString> &jetString_) {
  jetString_.push_back("0j");
  jetString_.push_back("1j");
  jetString_.push_back("2j");
  jetString_.push_back("nj");
  return;
}
void lepStrings(vector<TString> &lepString_) {
  lepString_.push_back("l");
  lepString_.push_back("e");
  lepString_.push_back("m");
  return;
}

