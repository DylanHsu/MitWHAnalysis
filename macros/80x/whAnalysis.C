#include "TMVA/Reader.h"
#include "MitWHAnalysis/macros/80x/whAnalysis.h"


void whAnalysis(
 string subdirectory="",
 TString bdtWeights="",
 bool isBatch=false,
 TString batchInput="",
 Int_t batchCategory=0
) {
  // Hardcoded settings
  double mcPrescale = 1.; 
  Double_t lumi = 35.9;
  TString filesPathDA   = "/data/t3home000/dhsu/panda/merged_skims/";
  TString filesPathMC   = "/data/t3home000/dhsu/panda/merged_skims/";
  TString weightsPath   = (isBatch? "" : "/home/dhsu/CMSSW_8_0_26_patch1/src/weights/");

  TString the_BDT_weights=weightsPath+bdtWeights;

  // Set up output dirs
  if(subdirectory!="" && subdirectory.c_str()[0]!='/') subdirectory = "/"+subdirectory;
  system(("mkdir -p MitWHAnalysis/datacards"+subdirectory).c_str());
  system(("mkdir -p MitWHAnalysis/plots"+subdirectory).c_str());

  // Process types
  vector<TString> categoryName_, jetString_, lepString_;
  categoryNames(categoryName_);
  jetStrings(jetString_);
  lepStrings(lepString_);
  
  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infileName_; vector<Int_t> infileCat_;  
  if(isBatch) {
    infileName_.push_back(batchInput);
    infileCat_.push_back(batchCategory);
  } else {
    // Data files
    infileName_.push_back(Form("%sSingleElectron.root"  ,  filesPathDA.Data()));  infileCat_.push_back(0);
    infileName_.push_back(Form("%sSingleMuon.root"      ,  filesPathDA.Data()));  infileCat_.push_back(0);
    // MC files
    // Top backgrounds
    infileName_.push_back(Form("%sTTbar_Powheg.root"    ,  filesPathMC.Data()));  infileCat_.push_back(1);
    infileName_.push_back(Form("%sSingleTop_tG.root"    ,  filesPathMC.Data()));  infileCat_.push_back(1);
    infileName_.push_back(Form("%sSingleTop_tT.root"    ,  filesPathMC.Data()));  infileCat_.push_back(1);
    infileName_.push_back(Form("%sSingleTop_tTbar.root" ,  filesPathMC.Data()));  infileCat_.push_back(1);
    infileName_.push_back(Form("%sSingleTop_tW.root"    ,  filesPathMC.Data()));  infileCat_.push_back(1);
    infileName_.push_back(Form("%sSingleTop_tbarW.root" ,  filesPathMC.Data()));  infileCat_.push_back(1);
    // Single boson production
    infileName_.push_back(Form("%sWJets_EWKWMinus.root"     , filesPathMC.Data()));  infileCat_.push_back(2);
    infileName_.push_back(Form("%sWJets_EWKWPlus.root"      , filesPathMC.Data()));  infileCat_.push_back(2);
    infileName_.push_back(Form("%sWJets_Wpt0to50.root"      , filesPathMC.Data()));  infileCat_.push_back(3);
    infileName_.push_back(Form("%sWJets_Wpt50to100.root"    , filesPathMC.Data()));  infileCat_.push_back(3);
    infileName_.push_back(Form("%sWJets_pt100to250.root"    , filesPathMC.Data()));  infileCat_.push_back(3);
    infileName_.push_back(Form("%sWJets_pt250to400.root"    , filesPathMC.Data()));  infileCat_.push_back(3);
    infileName_.push_back(Form("%sWJets_pt400to600.root"    , filesPathMC.Data()));  infileCat_.push_back(3);
    infileName_.push_back(Form("%sWJets_pt600toinf.root"    , filesPathMC.Data()));  infileCat_.push_back(3);
    infileName_.push_back(Form("%sZJets_pt50to100"          , filesPathMC.Data()));  infileCat_.push_back(4);
    infileName_.push_back(Form("%sZJets_pt100to250"         , filesPathMC.Data()));  infileCat_.push_back(4);
    infileName_.push_back(Form("%sZJets_pt250to400"         , filesPathMC.Data()));  infileCat_.push_back(4);
    infileName_.push_back(Form("%sZJets_pt400to650"         , filesPathMC.Data()));  infileCat_.push_back(4);
    infileName_.push_back(Form("%sZJets_pt650toinf"         , filesPathMC.Data()));  infileCat_.push_back(4);
    infileName_.push_back(Form("%sZJets_EWK.root"           , filesPathMC.Data()));  infileCat_.push_back(4);
    infileName_.push_back(Form("%sZtoNuNu_EWK.root"         , filesPathMC.Data()));  infileCat_.push_back(4);
    infileName_.push_back(Form("%sZtoNuNu_Zpt50to100.root"  , filesPathMC.Data()));  infileCat_.push_back(4);
    infileName_.push_back(Form("%sZtoNuNu_Zpt100to250.root" , filesPathMC.Data()));  infileCat_.push_back(4);
    infileName_.push_back(Form("%sZtoNuNu_Zpt250to400.root" , filesPathMC.Data()));  infileCat_.push_back(4);
    infileName_.push_back(Form("%sZtoNuNu_Zpt400to650.root" , filesPathMC.Data()));  infileCat_.push_back(4);
    infileName_.push_back(Form("%sZtoNuNu_Zpt650toinf.root" , filesPathMC.Data()));  infileCat_.push_back(4);
    // Diboson production
    infileName_.push_back(Form("%sDiboson_wz.root"      ,  filesPathMC.Data()));  infileCat_.push_back(5);
    infileName_.push_back(Form("%sDiboson_zz.root"      ,  filesPathMC.Data()));  infileCat_.push_back(6);
    infileName_.push_back(Form("%sDiboson_ww.root"      ,  filesPathMC.Data()));  infileCat_.push_back(7);
    // Random crap samples
    infileName_.push_back(Form("%sQCD.root"             ,  filesPathMC.Data()));  infileCat_.push_back(8);
    infileName_.push_back(Form("%sGJets.root"           ,  filesPathMC.Data()));  infileCat_.push_back(8);
    // Signal samples
    infileName_.push_back(Form("%sWminusH_HToInvisible_WToLNu_M125_13TeV_powheg_pythia8.root"           ,  filesPathMC.Data()));  infileCat_.push_back(9);
    infileName_.push_back(Form("%sWplusH_HToInvisible_WToLNu_M125_13TeV_powheg_pythia8.root"            ,  filesPathMC.Data()));  infileCat_.push_back(9);
  }

  //Declaration of leaves types
  Int_t nJet, nJot, nLooseLep, nLooseElectron, nLooseMuon, nTightLep, nTightElectron, nTightMuon, looseLep1PdgId, looseLep2PdgId, looseLep1IsTight, looseLep2IsTight, nTau, jot1VBFID, jetNMBtags, nB, isGS, looseLep1IsHLTSafe, looseLep2IsHLTSafe, runNumber, lumiNumber, npv, pu, trigger, metFilter, egmFilter, genTopIsHad, genAntiTopIsHad, nIsoJet, jet1Flav, jet1IsTight, jet2Flav, isojet1Flav, isojet2Flav, jetNBtags, isojetNBtags, nFatjet, nHF, nLoosePhoton, nTightPhoton, loosePho1IsTight, nAwayPFCH, nTransversePFCH, nTowardPFCH;
  Float_t jot1Phi, jot1Pt, jot1GenPt, jot1Eta, jot2Phi, jot2Pt, jot2GenPt, jot2Eta, jot12DPhi, jot12Mass, jot12DEta, pfmetUp, pfmetDown, pfUWmagUp, pfUZmagUp, pfUAmagUp, pfUmagUp, pfUWmagDown, pfUZmagDown, pfUAmagDown, pfUmagDown, jot1EtaUp, jot1EtaDown, jot1PtUp, jot1PtDown, jot2PtUp, jot2PtDown, jot12MassUp, jot12DEtaUp, jot12DPhiUp, jot12MassDown, jot12DEtaDown, jot12DPhiDown, jot2EtaUp, jot2EtaDown, scale[6], sf_btag0, sf_btag1, sf_btag2, sf_btagGT0, sf_sjbtag0, sf_sjbtag1, sf_sjbtag2, sf_sjbtagGT0, sf_btag0BUp, sf_btag1BUp, sf_btag2BUp, sf_btagGT0BUp, sf_sjbtag0BUp, sf_sjbtag1BUp, sf_sjbtag2BUp, sf_sjbtagGT0BUp, sf_btag0BDown, sf_btag1BDown, sf_btag2BDown, sf_btagGT0BDown, sf_sjbtag0BDown, sf_sjbtag1BDown, sf_sjbtag2BDown, sf_sjbtagGT0BDown, sf_btag0MUp, sf_btag1MUp, sf_btag2MUp, sf_btagGT0MUp, sf_sjbtag0MUp, sf_sjbtag1MUp, sf_sjbtag2MUp, sf_sjbtagGT0MUp, sf_btag0MDown, sf_btag1MDown, sf_btag2MDown, sf_btagGT0MDown, sf_sjbtag0MDown, sf_sjbtag1MDown, sf_sjbtag2MDown, sf_sjbtagGT0MDown, sf_metTrigZmm, sf_qcdV_VBF, pfmetRaw, mcWeight, normalizedWeight, filter_maxRecoil, filter_whichRecoil, sf_ewkV, sf_qcdV, sf_ewkV2j, sf_qcdV2j, sf_qcdTT, sf_lepID, sf_lepIso, sf_lepTrack, sf_pho, sf_eleTrig, sf_phoTrig, sf_metTrig, sf_pu, sf_npv, sf_tt, sf_tt_ext, sf_tt_bound, sf_tt8TeV, sf_tt8TeV_ext, sf_tt8TeV_bound, sf_phoPurity, pfmet, pfmetphi, pfmetnomu, puppimet, puppimetphi, calomet, calometphi, pfcalobalance, sumET, trkmet, puppiUWmag, puppiUWphi, puppiUZmag, puppiUZphi, puppiUAmag, puppiUAphi, puppiUperp, puppiUpara, puppiUmag, puppiUphi, pfUWmag, pfUWphi, pfUZmag, pfUZphi, pfUAmag, pfUAphi, pfUperp, pfUpara, pfUmag, pfUphi, dphipfmet, dphipuppimet, dphipuppiUW, dphipuppiUZ, dphipuppiUA, dphipfUW, dphipfUZ, dphipfUA, dphipuppiU, dphipfU, trueGenBosonPt, genBosonPt, genBosonEta, genBosonMass, genBosonPhi, genWPlusPt, genWMinusPt, genWPlusEta, genWMinusEta, genTopPt, genTopEta, genAntiTopPt, genAntiTopEta, genTTPt, genTTEta, jet1Phi, jet1Pt, jet1GenPt, jet1Eta, jet1CSV, jet2Phi, jet2Pt, jet2GenPt, jet2Eta, jet2CSV, isojet1Pt, isojet1CSV, isojet2Pt, isojet2CSV, loosePho1Pt, loosePho1Eta, loosePho1Phi, looseLep1Pt, looseLep1Eta, looseLep1Phi, looseLep2Pt, looseLep2Eta, looseLep2Phi, diLepMass, mT, scaleUp, scaleDown, pdfUp, pdfDown, sumPtAwayPFCH, sumPtTransversePFCH, sumPtTowardPFCH, PFCH1Pt, PFCH1Eta, PFCH1Phi, PFCH2Pt, PFCH2Eta, PFCH2Phi, PFCH3Pt, PFCH3Eta, PFCH3Phi;
  ULong64_t eventNumber;
  
  //*******************************************************
  // Set up BDT (if applicable)
  //*******************************************************
  TMVA::Reader *reader; // =new TMVA::Reader();
  float *mva_pfmet, *mva_pfmetphi, *mva_mT, *mva_looseLep1Pt, *mva_looseLep1Eta, *mva_looseLep1Phi, *mva_pfUmag, *mva_pfUphi, *mva_sumPtAwayPFCH, *mva_sumPtTransversePFCH, *mva_sumPtTowardPFCH, *mva_PFCH1Pt, *mva_PFCH2Pt, *mva_PFCH3Pt, *mva_PFCH1Eta, *mva_PFCH2Eta, *mva_PFCH3Eta, *mva_PFCH1Phi, *mva_PFCH2Phi, *mva_PFCH3Phi;
  float mva_nAwayPFCH, mva_nTransversePFCH, mva_nTowardPFCH;
  
  bool useBDT = (MVAVarType==2);
  if(useBDT) {
    mva_pfmet               = &pfmet               ;
    mva_pfmetphi            = &pfmetphi            ;
    mva_mT                  = &mT                  ;
    mva_looseLep1Pt         = &looseLep1Pt         ;
    mva_looseLep1Eta        = &looseLep1Eta        ;
    mva_looseLep1Phi        = &looseLep1Phi        ;
    mva_pfUmag              = &pfUmag              ;
    mva_pfUphi              = &pfUphi              ;
    mva_sumPtAwayPFCH       = &sumPtAwayPFCH       ;
    mva_sumPtTransversePFCH = &sumPtTransversePFCH ;
    mva_sumPtTowardPFCH     = &sumPtTowardPFCH     ;
    mva_PFCH1Pt             = &PFCH1Pt             ;
    mva_PFCH2Pt             = &PFCH2Pt             ;
    mva_PFCH3Pt             = &PFCH3Pt             ;
    mva_PFCH1Eta            = &PFCH1Eta            ;
    mva_PFCH2Eta            = &PFCH2Eta            ;
    mva_PFCH3Eta            = &PFCH3Eta            ;
    mva_PFCH1Phi            = &PFCH1Phi            ;
    mva_PFCH2Phi            = &PFCH2Phi            ;
    mva_PFCH3Phi            = &PFCH3Phi            ;
    reader=new TMVA::Reader();
    if(MVAVarType==2) {
      reader->AddVariable( "pfmet"                  , mva_pfmet               );
      reader->AddVariable( "pfmetphi"               , mva_pfmetphi            );
      reader->AddVariable( "mT"                     , mva_mT                  );
      reader->AddVariable( "looseLep1Pt"            , mva_looseLep1Pt         );
      reader->AddVariable( "looseLep1Eta"           , mva_looseLep1Eta        );
      reader->AddVariable( "looseLep1Phi"           , mva_looseLep1Phi        );
      reader->AddVariable( "pfUmag"                 , mva_pfUmag              );
      reader->AddVariable( "pfUphi"                 , mva_pfUphi              );
      reader->AddVariable( "nAwayPFCH"              , &mva_nAwayPFCH          );
      reader->AddVariable( "nTransversePFCH"        , &mva_nTransversePFCH    );
      reader->AddVariable( "nTowardPFCH"            , &mva_nTowardPFCH        );
      reader->AddVariable( "sumPtAwayPFCH"          , mva_sumPtAwayPFCH       );
      reader->AddVariable( "sumPtTransversePFCH"    , mva_sumPtTransversePFCH );
      reader->AddVariable( "sumPtTowardPFCH"        , mva_sumPtTowardPFCH     );
      reader->AddVariable( "PFCH1Pt"                , mva_PFCH1Pt             );
      reader->AddVariable( "PFCH2Pt"                , mva_PFCH2Pt             );
      reader->AddVariable( "PFCH3Pt"                , mva_PFCH3Pt             );
      reader->AddVariable( "PFCH1Eta"               , mva_PFCH1Eta            );
      reader->AddVariable( "PFCH2Eta"               , mva_PFCH2Eta            );
      reader->AddVariable( "PFCH3Eta"               , mva_PFCH3Eta            );
      reader->AddVariable( "PFCH1Phi"               , mva_PFCH1Phi            );
      reader->AddVariable( "PFCH2Phi"               , mva_PFCH2Phi            );
      reader->AddVariable( "PFCH3Phi"               , mva_PFCH3Phi            );
    }
    reader->BookMVA("BDT", the_BDT_weights);
  }
  
  //*******************************************************
  // Set up trees (if batch mode)
  //*******************************************************
  char output[400];
  if(isBatch) sprintf(output,"histo_wh_nice_%s",batchInput.Data());
  else        sprintf(output,"MitWHAnalysis/plots%s/histo_wh_nice.root",subdirectory.c_str());
  TFile *output_plots = new TFile(output,"RECREATE","",0);
  UInt_t nBinPlot; Float_t xminPlot, xmaxPlot;
  UChar_t batch_nJot, batch_flavor, batch_category;
  UInt_t thePlot;
  Float_t theVar, totalWeight;
  char batch_plotName[64];
  TTree *tree_batchPlotIndex, *tree_batchPlots;
  if(isBatch) {
    tree_batchPlotIndex=new TTree("tree_batchPlotIndex", "tree_batchPlotIndex");
    tree_batchPlotIndex->Branch("thePlot" , &thePlot       );
    tree_batchPlotIndex->Branch("nBinPlot", &nBinPlot      );
    tree_batchPlotIndex->Branch("xminPlot", &xminPlot      );
    tree_batchPlotIndex->Branch("xmaxPlot", &xmaxPlot      );
    tree_batchPlotIndex->Branch("plotName", batch_plotName, "plotName/C", 64);
    tree_batchPlots=new TTree("tree_batchPlots", "tree_batchPlots");
    tree_batchPlots->Branch("thePlot"     , &thePlot        );
    tree_batchPlots->Branch("nJets"       , &batch_nJot     );
    tree_batchPlots->Branch("flavor"      , &batch_flavor   );
    tree_batchPlots->Branch("category"    , &batch_category );
    tree_batchPlots->Branch("theVar"      , &theVar         );
    tree_batchPlots->Branch("totalWeight" , &totalWeight    );
  }

  //*******************************************************
  // Set up histograms
  //*******************************************************
  TH1D *histo[allPlots][process_types][jet_cats][3];
  vector<TString> plotName_; TString plotName;
  for(thePlot=0; thePlot<allPlots; thePlot++){
    //if     (thePlot ==  0) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="presel lep pT"             ;}
    //else if(thePlot ==  1) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="presel mT"                 ;}
    //else if(thePlot ==  2) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="presel MET"                ;}
    //else if(thePlot ==  3) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel lep eta"            ;}
    //else if(thePlot ==  4) {nBinPlot = 185; xminPlot = 30.; xmaxPlot = 400.0;       plotName="presel leading jet pT"     ;}
    //else if(thePlot ==  5) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;       plotName="presel max CSV2"           ;}
    //else if(thePlot ==  6) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =   4.0;       plotName="presel pf calo balance"    ;}
    //else if(thePlot ==  7) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel lepton MET dPhi"    ;}
    //else if(thePlot ==  8) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel jet MET dPhi"       ;}
    //else if(thePlot ==  9) {nBinPlot =2000; xminPlot = 0.0; xmaxPlot =  10.0;       plotName="presel jet lep dR"         ;}
    //else if(thePlot == 10) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   100;       plotName="presel NPV "               ;}
    //else if(thePlot == 11) {nBinPlot =1000; xminPlot =-500; xmaxPlot =   500;       plotName="presel Uperp"              ;}
    //else if(thePlot == 12) {nBinPlot =1000; xminPlot =-500; xmaxPlot =   500;       plotName="presel Upara"              ;}
    //else if(thePlot == 13) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel nAwayPFCH"          ;}
    //else if(thePlot == 14) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel nTransversePFCH"    ;}
    //else if(thePlot == 15) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel nTowardPFCH"        ;}
    //else if(thePlot == 16) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel sumPtAwayPFCH"      ;}
    //else if(thePlot == 17) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel sumPtTransversePFCH";}
    //else if(thePlot == 18) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="presel sumPtTowardPFCH"    ;}
    //else if(thePlot == 19) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="presel PFCH1Pt"            ;}
    //else if(thePlot == 20) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="presel PFCH1Eta"           ;}
    //else if(thePlot == 21) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel PFCH1Phi"           ;}
    //else if(thePlot == 22) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="presel PFCH2Pt"            ;}
    //else if(thePlot == 23) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="presel PFCH2Eta"           ;}
    //else if(thePlot == 24) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel PFCH2Phi"           ;}
    //else if(thePlot == 25) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="presel PFCH3Pt"            ;}
    //else if(thePlot == 26) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="presel PFCH3Eta"           ;}
    //else if(thePlot == 27) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="presel PFCH3Phi"           ;}
    //else if(thePlot == 28) {nBinPlot =2000; xminPlot = xbins[0]; xmaxPlot = xbins[nBinMVA-1]; plotName="presel MVAVar"   ;}
    //else if(thePlot == 29) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="sigsel lep pT"             ;}
    //else if(thePlot == 30) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="sigsel mT"                 ;}
    //else if(thePlot == 31) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="sigsel lep eta"            ;}
    //else if(thePlot == 32) {nBinPlot = 185; xminPlot = 30.; xmaxPlot = 400.0;       plotName="sigsel leading jet pT"     ;}
    //else if(thePlot == 33) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =   4.0;       plotName="sigsel pf calo balance"    ;}
    //else if(thePlot == 34) {nBinPlot =2000; xminPlot = 0.0; xmaxPlot =  10.0;       plotName="sigsel jet lep dR"         ;}
    //else if(thePlot == 35) {nBinPlot =1000; xminPlot =-500; xmaxPlot =   500;       plotName="sigsel Uperp"              ;}
    //else if(thePlot == 36) {nBinPlot =1000; xminPlot =-500; xmaxPlot =   500;       plotName="sigsel Upara"              ;}
    //else if(thePlot == 37) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel nAwayPFCH"          ;}
    //else if(thePlot == 38) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel nTransversePFCH"    ;}
    //else if(thePlot == 39) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel nTowardPFCH"        ;}
    //else if(thePlot == 40) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel sumPtAwayPFCH"      ;}
    //else if(thePlot == 41) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel sumPtTransversePFCH";}
    //else if(thePlot == 42) {nBinPlot =1000; xminPlot = 0.0; xmaxPlot =   500;       plotName="sigsel sumPtTowardPFCH"    ;}
    //else if(thePlot == 43) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="sigsel PFCH1Pt"            ;}
    //else if(thePlot == 44) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="sigsel PFCH1Eta"           ;}
    //else if(thePlot == 45) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="sigsel PFCH1Phi"           ;}
    //else if(thePlot == 46) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="sigsel PFCH2Pt"            ;}
    //else if(thePlot == 47) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="sigsel PFCH2Eta"           ;}
    //else if(thePlot == 48) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="sigsel PFCH2Phi"           ;}
    //else if(thePlot == 49) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   100;       plotName="sigsel PFCH3Pt"            ;}
    //else if(thePlot == 50) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =   2.5;       plotName="sigsel PFCH3Eta"           ;}
    //else if(thePlot == 51) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="sigsel PFCH3Phi"           ;}
    //else if(thePlot == 52) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =   4.0;       plotName="N-1 lepton MET balance"    ;}
    //else if(thePlot == 53) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;       plotName="N-1 max CSV2"              ;}
    //else if(thePlot == 54) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =1000.0;       plotName="N-1 MET"                   ;}
    //else if(thePlot == 55) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="N-1 lepton MET dPhi"       ;}
    //else if(thePlot == 56) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi(); plotName="N-1 jet MET dPhi"          ;}
    //else if(thePlot == allPlots-1) plotName="Shape analysis";

    plotName_.push_back(plotName);
    TH1D* histos = makeHisto(thePlot, plotName);
    if(isBatch) {
      sprintf(batch_plotName,"%s",plotName_[thePlot].Data());
      tree_batchPlotIndex->Fill();
    } else {
      //if(thePlot != allPlots-1) histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
      //else                      histos = new TH1D("histos", "histos", nBinMVA, xbins);
      
      //histos->Sumw2();
      for(unsigned int i_jet=0; i_jet<jet_cats; i_jet++) { for(unsigned int i_type=0; i_type<process_types; i_type++) { for(unsigned int i_flav=0; i_flav<3; i_flav++) {
        histo[thePlot][i_type][i_jet][i_flav] = (TH1D*) histos->Clone(Form("%s %s (%s,%s)",plotName_[thePlot].Data(), categoryName_[i_type].Data(), jetString_[i_jet].Data(), lepString_[i_flav].Data()));
        histo[thePlot][i_type][i_jet][i_flav]->SetDirectory(0);
      }}}
      histos->Reset();histos->Clear();
    }
  }

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infileName_.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infileName_[ifile].Data());

    TFile *the_input_file = TFile::Open(infileName_[ifile].Data());
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
    UChar_t theCategory = infileCat_[ifile];
    // Set branch addresses.
    the_input_tree->SetBranchAddress("nJet",&nJet);
    the_input_tree->SetBranchAddress("nJot",&nJot);
    the_input_tree->SetBranchAddress("jot1Phi",&jot1Phi);
    the_input_tree->SetBranchAddress("jot1Pt",&jot1Pt);
    the_input_tree->SetBranchAddress("jot1GenPt",&jot1GenPt);
    the_input_tree->SetBranchAddress("jot1Eta",&jot1Eta);
    the_input_tree->SetBranchAddress("jot2Phi",&jot2Phi);
    the_input_tree->SetBranchAddress("jot2Pt",&jot2Pt);
    the_input_tree->SetBranchAddress("jot2GenPt",&jot2GenPt);
    the_input_tree->SetBranchAddress("jot2Eta",&jot2Eta);
    the_input_tree->SetBranchAddress("jot12DPhi",&jot12DPhi);
    the_input_tree->SetBranchAddress("jot12Mass",&jot12Mass);
    the_input_tree->SetBranchAddress("jot12DEta",&jot12DEta);
    the_input_tree->SetBranchAddress("pfmetUp",&pfmetUp);
    the_input_tree->SetBranchAddress("pfmetDown",&pfmetDown);
    //the_input_tree->SetBranchAddress("pfUWmagUp",&pfUWmagUp);
    //the_input_tree->SetBranchAddress("pfUZmagUp",&pfUZmagUp);
    //the_input_tree->SetBranchAddress("pfUAmagUp",&pfUAmagUp);
    //the_input_tree->SetBranchAddress("pfUmagUp",&pfUmagUp);
    //the_input_tree->SetBranchAddress("pfUWmagDown",&pfUWmagDown);
    //the_input_tree->SetBranchAddress("pfUZmagDown",&pfUZmagDown);
    //the_input_tree->SetBranchAddress("pfUAmagDown",&pfUAmagDown);
    //the_input_tree->SetBranchAddress("pfUmagDown",&pfUmagDown);
    //the_input_tree->SetBranchAddress("jot1EtaUp",&jot1EtaUp);
    //the_input_tree->SetBranchAddress("jot1EtaDown",&jot1EtaDown);
    //the_input_tree->SetBranchAddress("jot1PtUp",&jot1PtUp);
    //the_input_tree->SetBranchAddress("jot1PtDown",&jot1PtDown);
    //the_input_tree->SetBranchAddress("jot2PtUp",&jot2PtUp);
    //the_input_tree->SetBranchAddress("jot2PtDown",&jot2PtDown);
    //the_input_tree->SetBranchAddress("jot12MassUp",&jot12MassUp);
    //the_input_tree->SetBranchAddress("jot12DEtaUp",&jot12DEtaUp);
    //the_input_tree->SetBranchAddress("jot12DPhiUp",&jot12DPhiUp);
    //the_input_tree->SetBranchAddress("jot12MassDown",&jot12MassDown);
    //the_input_tree->SetBranchAddress("jot12DEtaDown",&jot12DEtaDown);
    //the_input_tree->SetBranchAddress("jot12DPhiDown",&jot12DPhiDown);
    //the_input_tree->SetBranchAddress("jot2EtaUp",&jot2EtaUp);
    //the_input_tree->SetBranchAddress("jot2EtaDown",&jot2EtaDown);
    //the_input_tree->SetBranchAddress("jot1VBFID",&jot1VBFID);
    the_input_tree->SetBranchAddress("jetNMBtags",&jetNMBtags);
    //the_input_tree->SetBranchAddress("pfmetRaw",&pfmetRaw);
    the_input_tree->SetBranchAddress("nB",&nB);
    the_input_tree->SetBranchAddress("isGS",&isGS);
    //the_input_tree->SetBranchAddress("looseLep1IsHLTSafe",&looseLep1IsHLTSafe);
    //the_input_tree->SetBranchAddress("looseLep2IsHLTSafe",&looseLep2IsHLTSafe);
    the_input_tree->SetBranchAddress("runNumber",&runNumber);
    the_input_tree->SetBranchAddress("lumiNumber",&lumiNumber);
    the_input_tree->SetBranchAddress("eventNumber",&eventNumber);
    the_input_tree->SetBranchAddress("npv",&npv);
    the_input_tree->SetBranchAddress("pu",&pu);
    //the_input_tree->SetBranchAddress("trigger",&trigger);
    //the_input_tree->SetBranchAddress("metFilter",&metFilter);
    //the_input_tree->SetBranchAddress("egmFilter",&egmFilter);
    //the_input_tree->SetBranchAddress("filter_maxRecoil",&filter_maxRecoil);
    //the_input_tree->SetBranchAddress("filter_whichRecoil",&filter_whichRecoil);
    the_input_tree->SetBranchAddress("pfmet",&pfmet);
    the_input_tree->SetBranchAddress("pfmetphi",&pfmetphi);
    //the_input_tree->SetBranchAddress("pfmetnomu",&pfmetnomu);
    //the_input_tree->SetBranchAddress("puppimet",&puppimet);
    //the_input_tree->SetBranchAddress("puppimetphi",&puppimetphi);
    the_input_tree->SetBranchAddress("calomet",&calomet);
    the_input_tree->SetBranchAddress("calometphi",&calometphi);
    the_input_tree->SetBranchAddress("pfcalobalance",&pfcalobalance);
    //the_input_tree->SetBranchAddress("sumET",&sumET);
    //the_input_tree->SetBranchAddress("trkmet",&trkmet);
    //the_input_tree->SetBranchAddress("puppiUWmag",&puppiUWmag);
    //the_input_tree->SetBranchAddress("puppiUWphi",&puppiUWphi);
    //the_input_tree->SetBranchAddress("puppiUZmag",&puppiUZmag);
    //the_input_tree->SetBranchAddress("puppiUZphi",&puppiUZphi);
    //the_input_tree->SetBranchAddress("puppiUAmag",&puppiUAmag);
    //the_input_tree->SetBranchAddress("puppiUAphi",&puppiUAphi);
    //the_input_tree->SetBranchAddress("puppiUperp",&puppiUperp);
    //the_input_tree->SetBranchAddress("puppiUpara",&puppiUpara);
    //the_input_tree->SetBranchAddress("puppiUmag",&puppiUmag);
    //the_input_tree->SetBranchAddress("puppiUphi",&puppiUphi);
    //the_input_tree->SetBranchAddress("pfUWmag",&pfUWmag);
    //the_input_tree->SetBranchAddress("pfUWphi",&pfUWphi);
    //the_input_tree->SetBranchAddress("pfUZmag",&pfUZmag);
    //the_input_tree->SetBranchAddress("pfUZphi",&pfUZphi);
    //the_input_tree->SetBranchAddress("pfUAmag",&pfUAmag);
    //the_input_tree->SetBranchAddress("pfUAphi",&pfUAphi);
    //the_input_tree->SetBranchAddress("pfUperp",&pfUperp);
    //the_input_tree->SetBranchAddress("pfUpara",&pfUpara);
    the_input_tree->SetBranchAddress("pfUmag",&pfUmag);
    the_input_tree->SetBranchAddress("pfUphi",&pfUphi);
    //the_input_tree->SetBranchAddress("dphipfmet",&dphipfmet);
    //the_input_tree->SetBranchAddress("dphipuppimet",&dphipuppimet);
    //the_input_tree->SetBranchAddress("dphipuppiUW",&dphipuppiUW);
    //the_input_tree->SetBranchAddress("dphipuppiUZ",&dphipuppiUZ);
    //the_input_tree->SetBranchAddress("dphipuppiUA",&dphipuppiUA);
    //the_input_tree->SetBranchAddress("dphipfUW",&dphipfUW);
    //the_input_tree->SetBranchAddress("dphipfUZ",&dphipfUZ);
    //the_input_tree->SetBranchAddress("dphipfUA",&dphipfUA);
    //the_input_tree->SetBranchAddress("dphipuppiU",&dphipuppiU);
    //the_input_tree->SetBranchAddress("dphipfU",&dphipfU);
    //the_input_tree->SetBranchAddress("trueGenBosonPt",&trueGenBosonPt);
    //the_input_tree->SetBranchAddress("genBosonPt",&genBosonPt);
    //the_input_tree->SetBranchAddress("genBosonEta",&genBosonEta);
    //the_input_tree->SetBranchAddress("genBosonMass",&genBosonMass);
    //the_input_tree->SetBranchAddress("genBosonPhi",&genBosonPhi);
    //the_input_tree->SetBranchAddress("genWPlusPt",&genWPlusPt);
    //the_input_tree->SetBranchAddress("genWMinusPt",&genWMinusPt);
    //the_input_tree->SetBranchAddress("genWPlusEta",&genWPlusEta);
    //the_input_tree->SetBranchAddress("genWMinusEta",&genWMinusEta);
    //the_input_tree->SetBranchAddress("genTopPt",&genTopPt);
    //the_input_tree->SetBranchAddress("genTopIsHad",&genTopIsHad);
    //the_input_tree->SetBranchAddress("genTopEta",&genTopEta);
    //the_input_tree->SetBranchAddress("genAntiTopPt",&genAntiTopPt);
    //the_input_tree->SetBranchAddress("genAntiTopIsHad",&genAntiTopIsHad);
    //the_input_tree->SetBranchAddress("genAntiTopEta",&genAntiTopEta);
    //the_input_tree->SetBranchAddress("genTTPt",&genTTPt);
    //the_input_tree->SetBranchAddress("genTTEta",&genTTEta);
    the_input_tree->SetBranchAddress("nIsoJet",&nIsoJet);
    the_input_tree->SetBranchAddress("jet1Flav",&jet1Flav);
    the_input_tree->SetBranchAddress("jet1Phi",&jet1Phi);
    the_input_tree->SetBranchAddress("jet1Pt",&jet1Pt);
    the_input_tree->SetBranchAddress("jet1GenPt",&jet1GenPt);
    the_input_tree->SetBranchAddress("jet1Eta",&jet1Eta);
    the_input_tree->SetBranchAddress("jet1CSV",&jet1CSV);
    the_input_tree->SetBranchAddress("jet1IsTight",&jet1IsTight);
    the_input_tree->SetBranchAddress("jet2Flav",&jet2Flav);
    the_input_tree->SetBranchAddress("jet2Phi",&jet2Phi);
    the_input_tree->SetBranchAddress("jet2Pt",&jet2Pt);
    the_input_tree->SetBranchAddress("jet2GenPt",&jet2GenPt);
    the_input_tree->SetBranchAddress("jet2Eta",&jet2Eta);
    the_input_tree->SetBranchAddress("jet2CSV",&jet2CSV);
    //the_input_tree->SetBranchAddress("isojet1Pt",&isojet1Pt);
    //the_input_tree->SetBranchAddress("isojet1CSV",&isojet1CSV);
    //the_input_tree->SetBranchAddress("isojet1Flav",&isojet1Flav);
    //the_input_tree->SetBranchAddress("isojet2Pt",&isojet2Pt);
    //the_input_tree->SetBranchAddress("isojet2CSV",&isojet2CSV);
    //the_input_tree->SetBranchAddress("isojet2Flav",&isojet2Flav);
    the_input_tree->SetBranchAddress("jetNBtags",&jetNBtags);
    //the_input_tree->SetBranchAddress("isojetNBtags",&isojetNBtags);
    //the_input_tree->SetBranchAddress("nFatjet",&nFatjet);
    //the_input_tree->SetBranchAddress("nHF",&nHF);
    //the_input_tree->SetBranchAddress("nLoosePhoton",&nLoosePhoton);
    //the_input_tree->SetBranchAddress("nTightPhoton",&nTightPhoton);
    //the_input_tree->SetBranchAddress("loosePho1IsTight",&loosePho1IsTight);
    //the_input_tree->SetBranchAddress("loosePho1Pt",&loosePho1Pt);
    //the_input_tree->SetBranchAddress("loosePho1Eta",&loosePho1Eta);
    //the_input_tree->SetBranchAddress("loosePho1Phi",&loosePho1Phi);
    the_input_tree->SetBranchAddress("nLooseLep",&nLooseLep);
    the_input_tree->SetBranchAddress("nLooseElectron",&nLooseElectron);
    the_input_tree->SetBranchAddress("nLooseMuon",&nLooseMuon);
    //the_input_tree->SetBranchAddress("nTightLep",&nTightLep);
    //the_input_tree->SetBranchAddress("nTightElectron",&nTightElectron);
    //the_input_tree->SetBranchAddress("nTightMuon",&nTightMuon);
    the_input_tree->SetBranchAddress("looseLep1PdgId",&looseLep1PdgId);
    //the_input_tree->SetBranchAddress("looseLep2PdgId",&looseLep2PdgId);
    the_input_tree->SetBranchAddress("looseLep1IsTight",&looseLep1IsTight);
    //the_input_tree->SetBranchAddress("looseLep2IsTight",&looseLep2IsTight);
    the_input_tree->SetBranchAddress("looseLep1Pt",&looseLep1Pt);
    the_input_tree->SetBranchAddress("looseLep1Eta",&looseLep1Eta);
    the_input_tree->SetBranchAddress("looseLep1Phi",&looseLep1Phi);
    the_input_tree->SetBranchAddress("looseLep2Pt",&looseLep2Pt);
    the_input_tree->SetBranchAddress("looseLep2Eta",&looseLep2Eta);
    the_input_tree->SetBranchAddress("looseLep2Phi",&looseLep2Phi);
    //the_input_tree->SetBranchAddress("diLepMass",&diLepMass);
    //the_input_tree->SetBranchAddress("nTau",&nTau);
    the_input_tree->SetBranchAddress("mT",&mT);
    the_input_tree->SetBranchAddress("nAwayPFCH",&nAwayPFCH);
    the_input_tree->SetBranchAddress("nTransversePFCH",&nTransversePFCH);
    the_input_tree->SetBranchAddress("nTowardPFCH",&nTowardPFCH);
    the_input_tree->SetBranchAddress("sumPtAwayPFCH",&sumPtAwayPFCH);
    the_input_tree->SetBranchAddress("sumPtTransversePFCH",&sumPtTransversePFCH);
    the_input_tree->SetBranchAddress("sumPtTowardPFCH",&sumPtTowardPFCH);
    the_input_tree->SetBranchAddress("PFCH1Pt",&PFCH1Pt);
    the_input_tree->SetBranchAddress("PFCH1Eta",&PFCH1Eta);
    the_input_tree->SetBranchAddress("PFCH1Phi",&PFCH1Phi);
    the_input_tree->SetBranchAddress("PFCH2Pt",&PFCH2Pt);
    the_input_tree->SetBranchAddress("PFCH2Eta",&PFCH2Eta);
    the_input_tree->SetBranchAddress("PFCH2Phi",&PFCH2Phi);
    the_input_tree->SetBranchAddress("PFCH3Pt",&PFCH3Pt);
    the_input_tree->SetBranchAddress("PFCH3Eta",&PFCH3Eta);
    the_input_tree->SetBranchAddress("PFCH3Phi",&PFCH3Phi);
    if(theCategory != 0) {
      the_input_tree->SetBranchAddress("mcWeight",&mcWeight);
      the_input_tree->SetBranchAddress("normalizedWeight",&normalizedWeight);
      //the_input_tree->SetBranchAddress("sf_ewkV",&sf_ewkV);
      //the_input_tree->SetBranchAddress("sf_qcdV",&sf_qcdV);
      //the_input_tree->SetBranchAddress("sf_ewkV2j",&sf_ewkV2j);
      //the_input_tree->SetBranchAddress("sf_qcdV2j",&sf_qcdV2j);
      //the_input_tree->SetBranchAddress("sf_qcdTT",&sf_qcdTT);
      the_input_tree->SetBranchAddress("sf_lepID",&sf_lepID);
      the_input_tree->SetBranchAddress("sf_lepIso",&sf_lepIso);
      the_input_tree->SetBranchAddress("sf_lepTrack",&sf_lepTrack);
      //the_input_tree->SetBranchAddress("sf_pho",&sf_pho);
      //the_input_tree->SetBranchAddress("sf_eleTrig",&sf_eleTrig);
      //the_input_tree->SetBranchAddress("sf_phoTrig",&sf_phoTrig);
      //the_input_tree->SetBranchAddress("sf_metTrig",&sf_metTrig);
      the_input_tree->SetBranchAddress("sf_pu",&sf_pu);
      the_input_tree->SetBranchAddress("sf_npv",&sf_npv);
      //the_input_tree->SetBranchAddress("sf_tt",&sf_tt);
      //the_input_tree->SetBranchAddress("sf_tt_ext",&sf_tt_ext);
      //the_input_tree->SetBranchAddress("sf_tt_bound",&sf_tt_bound);
      //the_input_tree->SetBranchAddress("sf_tt8TeV",&sf_tt8TeV);
      //the_input_tree->SetBranchAddress("sf_tt8TeV_ext",&sf_tt8TeV_ext);
      //the_input_tree->SetBranchAddress("sf_tt8TeV_bound",&sf_tt8TeV_bound);
      //the_input_tree->SetBranchAddress("sf_phoPurity",&sf_phoPurity);
      the_input_tree->SetBranchAddress("scale",scale);
      //the_input_tree->SetBranchAddress("sf_btag0",&sf_btag0);
      //the_input_tree->SetBranchAddress("sf_btag1",&sf_btag1);
      //the_input_tree->SetBranchAddress("sf_btag2",&sf_btag2);
      //the_input_tree->SetBranchAddress("sf_btagGT0",&sf_btagGT0);
      //the_input_tree->SetBranchAddress("sf_sjbtag0",&sf_sjbtag0);
      //the_input_tree->SetBranchAddress("sf_sjbtag1",&sf_sjbtag1);
      //the_input_tree->SetBranchAddress("sf_sjbtag2",&sf_sjbtag2);
      //the_input_tree->SetBranchAddress("sf_sjbtagGT0",&sf_sjbtagGT0);
      //the_input_tree->SetBranchAddress("sf_btag0BUp",&sf_btag0BUp);
      //the_input_tree->SetBranchAddress("sf_btag1BUp",&sf_btag1BUp);
      //the_input_tree->SetBranchAddress("sf_btag2BUp",&sf_btag2BUp);
      //the_input_tree->SetBranchAddress("sf_btagGT0BUp",&sf_btagGT0BUp);
      //the_input_tree->SetBranchAddress("sf_sjbtag0BUp",&sf_sjbtag0BUp);
      //the_input_tree->SetBranchAddress("sf_sjbtag1BUp",&sf_sjbtag1BUp);
      //the_input_tree->SetBranchAddress("sf_sjbtag2BUp",&sf_sjbtag2BUp);
      //the_input_tree->SetBranchAddress("sf_sjbtagGT0BUp",&sf_sjbtagGT0BUp);
      //the_input_tree->SetBranchAddress("sf_btag0BDown",&sf_btag0BDown);
      //the_input_tree->SetBranchAddress("sf_btag1BDown",&sf_btag1BDown);
      //the_input_tree->SetBranchAddress("sf_btag2BDown",&sf_btag2BDown);
      //the_input_tree->SetBranchAddress("sf_btagGT0BDown",&sf_btagGT0BDown);
      //the_input_tree->SetBranchAddress("sf_sjbtag0BDown",&sf_sjbtag0BDown);
      //the_input_tree->SetBranchAddress("sf_sjbtag1BDown",&sf_sjbtag1BDown);
      //the_input_tree->SetBranchAddress("sf_sjbtag2BDown",&sf_sjbtag2BDown);
      //the_input_tree->SetBranchAddress("sf_sjbtagGT0BDown",&sf_sjbtagGT0BDown);
      //the_input_tree->SetBranchAddress("sf_btag0MUp",&sf_btag0MUp);
      //the_input_tree->SetBranchAddress("sf_btag1MUp",&sf_btag1MUp);
      //the_input_tree->SetBranchAddress("sf_btag2MUp",&sf_btag2MUp);
      //the_input_tree->SetBranchAddress("sf_btagGT0MUp",&sf_btagGT0MUp);
      //the_input_tree->SetBranchAddress("sf_sjbtag0MUp",&sf_sjbtag0MUp);
      //the_input_tree->SetBranchAddress("sf_sjbtag1MUp",&sf_sjbtag1MUp);
      //the_input_tree->SetBranchAddress("sf_sjbtag2MUp",&sf_sjbtag2MUp);
      //the_input_tree->SetBranchAddress("sf_sjbtagGT0MUp",&sf_sjbtagGT0MUp);
      //the_input_tree->SetBranchAddress("sf_btag0MDown",&sf_btag0MDown);
      //the_input_tree->SetBranchAddress("sf_btag1MDown",&sf_btag1MDown);
      //the_input_tree->SetBranchAddress("sf_btag2MDown",&sf_btag2MDown);
      //the_input_tree->SetBranchAddress("sf_btagGT0MDown",&sf_btagGT0MDown);
      //the_input_tree->SetBranchAddress("sf_sjbtag0MDown",&sf_sjbtag0MDown);
      //the_input_tree->SetBranchAddress("sf_sjbtag1MDown",&sf_sjbtag1MDown);
      //the_input_tree->SetBranchAddress("sf_sjbtag2MDown",&sf_sjbtag2MDown);
      //the_input_tree->SetBranchAddress("sf_sjbtagGT0MDown",&sf_sjbtagGT0MDown);
      //the_input_tree->SetBranchAddress("sf_metTrigZmm",&sf_metTrigZmm);
      //the_input_tree->SetBranchAddress("sf_qcdV_VBF",&sf_qcdV_VBF);
      the_input_tree->SetBranchAddress("scaleUp",&scaleUp);
      the_input_tree->SetBranchAddress("scaleDown",&scaleDown);
      the_input_tree->SetBranchAddress("pdfUp",&pdfUp);
      the_input_tree->SetBranchAddress("pdfDown",&pdfDown);
    }

    Long64_t nentries = the_input_tree->GetEntries();

    Long64_t nbytes = 0;
    double theMCPrescale = mcPrescale;
    if(theCategory == 0) theMCPrescale = 1.0;
    // Loop over events in the file
    for (int i=0; i<int(nentries/theMCPrescale); ++i) {
      if(i%100000==0 || i+1==int(nentries/theMCPrescale) ) printf("event %d out of %d\n",i+1,int(nentries/theMCPrescale));
      the_input_tree->GetEntry(i);
      if(nLooseLep<1) continue;
      
      // Set MVA input variables
      mva_nAwayPFCH           = float(nAwayPFCH           );
      mva_nTransversePFCH     = float(nTransversePFCH     );
      mva_nTowardPFCH         = float(nTowardPFCH         );

      // Analysis calculations
      unsigned int flavor=0; {
        unsigned int absPdgId = TMath::Abs(looseLep1PdgId);
        if     (absPdgId==11) flavor=1;
        else if(absPdgId==13) flavor=2;
      }
      double ptFrac = looseLep1Pt/pfmet;
      double bDiscrMax = TMath::Max(Float_t(0.), Float_t(TMath::Max(jet1CSV, jet2CSV)));
      double pfCaloBalance = TMath::Abs(pfmet-calomet)/pfmet;
      double deltaPhiLepMET, deltaPhiJetMET=-1, deltaRJetLep=-1; {
        TVector2 vLep1, vPFMET;
        vLep1.SetMagPhi(looseLep1Pt, looseLep1Phi);
        vPFMET.SetMagPhi(pfmet, pfmetphi);
        deltaPhiLepMET = TMath::Abs(vLep1.DeltaPhi(vPFMET));
        if(nJot>0) { 
          TVector3 vLep1, vJet1;
          vLep1.SetPtEtaPhi(looseLep1Pt, looseLep1Eta, looseLep1Phi);
          vJet1.SetPtEtaPhi(jot1Pt, jot1Eta, jot1Phi);
          deltaPhiJetMET=TMath::Abs(vPFMET.DeltaPhi(vJet1.XYvector()));
          deltaRJetLep=vLep1.DeltaR(vJet1);
        }
      }
      double Upara = pfUmag*TMath::Cos(pfUphi - looseLep1Phi);
      double Uperp = pfUmag*TMath::Sin(pfUphi - looseLep1Phi);
      double bdt_value=-1;
      if(useBDT) {
        if(MVAVarType==2) bdt_value = (reader->EvaluateMulticlass(multiclassSignal, "BDT"));
        //printf("bdt_value=%f\n", bdt_value);
      }

      double MVAVar=-1;
      if     (MVAVarType==1) MVAVar=pfmet;
      else if(MVAVarType==2) MVAVar=bdt_value;

      // Analysis booleans
      bool pass1LepSel, passNjets, passLepPt, passMT, passMet, passMetTight, passPtFrac, passDPhiLepMet, passDPhiJetMet, passBveto;
      switch(MVAVarType) {
        case 1:
          pass1LepSel    = (nLooseLep==1 && looseLep1IsTight && (flavor>0));
          passNjets      = (unsigned(nJot) < (jet_cats-1));
          passLepPt      = (looseLep1Pt>=50);
          passMT         = (mT>=50);
          passMet        = (pfmet>=50);
          passMetTight   = (pfmet>=100);
          passPtFrac     = (ptFrac > 0.4 && ptFrac < 1.5);
          passDPhiLepMet = (deltaPhiLepMET > 2.5);
          passDPhiJetMet = (deltaPhiJetMET > 0.5) || (deltaPhiJetMET == -1);
          passBveto      = jetNBtags==0 && jet1CSV<0.5426 && jet2CSV<0.5426;
          break;
        case 2:
          // testing 0-jet category for multiclass BDT
          pass1LepSel    = (nLooseLep==1 && looseLep1IsTight && (flavor>0));
          passNjets      = (unsigned(nJot) < (jet_cats-1));
          passLepPt      = (looseLep1Pt>=50);
          passMT         = (mT>=50);
          passMet        = (pfmet>=50);
          passMetTight   = (pfmet>=100);
          passPtFrac     = (ptFrac > 0.4 && ptFrac < 1.5);
          passDPhiLepMet = (deltaPhiLepMET > 2);
          passDPhiJetMet = true;
          passBveto      = true;
          break;
      } 
      
      // Selection booleans
      std::map<TString, bool> passAllCuts, passNMinusOne;
      passAllCuts["presel"] = pass1LepSel && passNjets && passLepPt && passMT && passMet && !passMetTight;
      passAllCuts["sigsel"] = pass1LepSel && passNjets && passLepPt && passMT &&             passMetTight  && passPtFrac && passDPhiLepMet && passDPhiJetMet && passBveto;

      passNMinusOne["pass1LepSel"   ] =                passNjets && passLepPt && passMT && passMetTight && passPtFrac && passDPhiLepMet && passDPhiJetMet && passBveto; 
      passNMinusOne["passNjets"     ] = pass1LepSel &&              passLepPt && passMT && passMetTight && passPtFrac && passDPhiLepMet && passDPhiJetMet && passBveto; 
      passNMinusOne["passLepPt"     ] = pass1LepSel && passNjets &&              passMT && passMetTight && passPtFrac && passDPhiLepMet && passDPhiJetMet && passBveto; 
      passNMinusOne["passMT"        ] = pass1LepSel && passNjets && passLepPt &&           passMetTight && passPtFrac && passDPhiLepMet && passDPhiJetMet && passBveto; 
      passNMinusOne["passMetTight"  ] = pass1LepSel && passNjets && passLepPt && passMT &&                 passPtFrac && passDPhiLepMet && passDPhiJetMet && passBveto; 
      passNMinusOne["passPtFrac"    ] = pass1LepSel && passNjets && passLepPt && passMT && passMetTight &&               passDPhiLepMet && passDPhiJetMet && passBveto; 
      passNMinusOne["passDPhiLepMet"] = pass1LepSel && passNjets && passLepPt && passMT && passMetTight && passPtFrac &&                   passDPhiJetMet && passBveto; 
      passNMinusOne["passDPhiJetMet"] = pass1LepSel && passNjets && passLepPt && passMT && passMetTight && passPtFrac && passDPhiLepMet &&                   passBveto; 
      passNMinusOne["passBveto"     ] = pass1LepSel && passNjets && passLepPt && passMT && passMetTight && passPtFrac && passDPhiLepMet && passDPhiJetMet             ; 
      
      // debug
      //if(passAllCuts["presel"]) printf("passed all presel cuts\n");
      //else printf("%d %d %d %d %d %d\n",int(pass1LepSel),int(passNjets),int(passLepPt),int(passMT),int(passMet),int(!passMetTight));

      //begin event weighting
      totalWeight = 1;
      if(theCategory != 0) {
        totalWeight *= normalizedWeight;
        totalWeight *= 1000. * lumi;
        totalWeight *= theMCPrescale;
        totalWeight *= sf_npv;
        //totalWeight *= sf_pu;
        totalWeight *= sf_lepID * sf_lepIso * sf_lepTrack;
      }
      // fill the plots
      for(thePlot=0; thePlot<allPlots; thePlot++){
        bool makePlot = false;
        if     (thePlot==  0 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=looseLep1Pt             ;} 
        else if(thePlot==  1 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=mT                      ;} 
        else if(thePlot==  2 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=pfmet                   ;} 
        else if(thePlot==  3 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=TMath::Abs(looseLep1Eta);} 
        else if(thePlot==  4 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=jot1Pt                  ;} 
        else if(thePlot==  5 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=bDiscrMax               ;} 
        else if(thePlot==  6 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=pfCaloBalance           ;} 
        else if(thePlot==  7 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=deltaPhiLepMET          ;} 
        else if(thePlot==  8 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=deltaPhiJetMET          ;} 
        else if(thePlot==  9 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=deltaRJetLep            ;} 
        else if(thePlot== 10 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=npv                     ;} 
        else if(thePlot== 11 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=Uperp                   ;} 
        else if(2 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=Upara                   ;} 
        else if(thePlot== 13 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=nAwayPFCH               ;} 
        else if(thePlot== 14 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=nTransversePFCH         ;} 
        else if(thePlot== 15 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=nTowardPFCH             ;} 
        else if(thePlot== 16 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=sumPtAwayPFCH           ;} 
        else if(thePlot== 17 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=sumPtTransversePFCH     ;} 
        else if(thePlot== 18 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=sumPtTowardPFCH         ;} 
        else if(thePlot== 19 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=PFCH1Pt                 ;} 
        else if(thePlot== 20 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=PFCH1Eta                ;} 
        else if(thePlot== 21 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=PFCH1Phi                ;} 
        else if(thePlot== 22 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=PFCH2Pt                 ;} 
        else if(thePlot== 23 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=PFCH2Eta                ;} 
        else if(thePlot== 24 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=PFCH2Phi                ;} 
        else if(thePlot== 25 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=PFCH3Pt                 ;} 
        else if(thePlot== 26 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=PFCH3Eta                ;} 
        else if(thePlot== 27 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=PFCH3Phi                ;} 
        else if(thePlot== 28 && passAllCuts  ["presel"        ]) {makePlot=true;theVar=MVAVar                  ;} 
        else if(thePlot== 29 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=looseLep1Pt             ;} 
        else if(thePlot== 30 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=mT                      ;} 
        else if(thePlot== 31 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=TMath::Abs(looseLep1Eta);} 
        else if(thePlot== 32 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=jot1Pt                  ;} 
        else if(thePlot== 33 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=pfCaloBalance           ;} 
        else if(thePlot== 34 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=deltaRJetLep            ;} 
        else if(thePlot== 35 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=Uperp                   ;} 
        else if(thePlot== 36 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=Upara                   ;} 
        else if(thePlot== 37 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=nAwayPFCH               ;} 
        else if(thePlot== 38 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=nTransversePFCH         ;} 
        else if(thePlot== 39 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=nTowardPFCH             ;} 
        else if(thePlot== 40 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=sumPtAwayPFCH           ;} 
        else if(thePlot== 41 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=sumPtTransversePFCH     ;} 
        else if(thePlot== 42 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=sumPtTowardPFCH         ;} 
        else if(thePlot== 43 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=PFCH1Pt                 ;} 
        else if(thePlot== 44 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=PFCH1Eta                ;} 
        else if(thePlot== 45 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=PFCH1Phi                ;} 
        else if(thePlot== 46 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=PFCH2Pt                 ;} 
        else if(thePlot== 47 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=PFCH2Eta                ;} 
        else if(thePlot== 48 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=PFCH2Phi                ;} 
        else if(thePlot== 49 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=PFCH3Pt                 ;} 
        else if(thePlot== 50 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=PFCH3Eta                ;} 
        else if(thePlot== 51 && passAllCuts  ["sigsel"        ]) {makePlot=true;theVar=PFCH3Phi                ;} 
        else if(thePlot== 52 && passNMinusOne["passPtFrac"    ]) {makePlot=true;theVar=ptFrac                  ;} 
        else if(thePlot== 53 && passNMinusOne["passBveto"     ]) {makePlot=true;theVar=bDiscrMax               ;} 
        else if(thePlot== 54 && passNMinusOne["passMetTight"  ]) {makePlot=true;theVar=pfmet                   ;} 
        else if(thePlot== 55 && passNMinusOne["passDPhiLepMet"]) {makePlot=true;theVar=deltaPhiLepMET          ;} 
        else if(thePlot== 56 && passNMinusOne["passDPhiJetMet"]) {makePlot=true;theVar=deltaPhiJetMET          ;} 
        else if(thePlot==allPlots-1 && passAllCuts["sigsel"   ]) {makePlot=true;theVar=MVAVar                  ;} 
        if(makePlot) {
          if(!isBatch) {
            nBinPlot = histo[thePlot][theCategory][0][0]->GetNbinsX();
            xminPlot = histo[thePlot][theCategory][0][0]->GetBinLowEdge(1);
            xmaxPlot = histo[thePlot][theCategory][0][0]->GetBinLowEdge(nBinPlot+1);
            theVar = TMath::Min(theVar, Float_t(xmaxPlot-0.001)); //overflow binning
            histo[thePlot][theCategory][nJot      ][flavor]->Fill( theVar, totalWeight);
            histo[thePlot][theCategory][nJot      ][0     ]->Fill( theVar, totalWeight);
            if(jet_cats==4) {
              histo[thePlot][theCategory][jet_cats-1][flavor]->Fill( theVar, totalWeight);
            }
          } else {
            batch_nJot=UChar_t(nJot);
            batch_flavor=UChar_t(flavor);
            batch_category=theCategory;
            tree_batchPlots->Fill();
          }
        }
      } // done looping over plots for filling
    } // done looping over entries
    the_input_file->Close();
  } // done with this flat tree file
  output_plots->cd();
  if(isBatch) {
    tree_batchPlotIndex->Write();
    tree_batchPlots->Write();
  } else {
    for(thePlot=0; thePlot<allPlots; thePlot++) {
      TDirectory *plotDir = output_plots->mkdir(plotName_[thePlot]);
      plotDir->cd();
      for(unsigned int i_jet=0; i_jet<jet_cats; i_jet++) {
        TDirectory *jetDir = plotDir->mkdir(jetString_[i_jet]);
        jetDir->cd();
        for(unsigned int i_flav=0; i_flav<3; i_flav++) {
          TDirectory *flavDir = jetDir->mkdir(lepString_[i_flav]);
          flavDir->cd();
          for(unsigned int i_type=0; i_type<process_types; i_type++) {
            histo[thePlot][i_type][i_jet][i_flav] -> Write(categoryName_[i_type]);
          }
        }
      }
    }
  }
  output_plots->Close();

}
void plotsFromBatchTree(
 TString batchTree,
 string subdirectory=""
) {
  printf("initializing...\n");
  if(subdirectory!="" && subdirectory.c_str()[0]!='/') subdirectory = "/"+subdirectory;
  system(("mkdir -p MitWHAnalysis/plots"+subdirectory).c_str());
  
  vector<TString> categoryName_, jetString_, lepString_;
  categoryNames(categoryName_);
  jetStrings(jetString_);
  lepStrings(lepString_);
  UInt_t nBinPlot; Float_t xminPlot, xmaxPlot;
  UChar_t batch_nJot, batch_flavor, batch_category;
  UInt_t thePlot;
  Float_t theVar, totalWeight;
  
  char output[400];
  sprintf(output,"MitWHAnalysis/plots%s/histo_wh_nice.root",subdirectory.c_str());
  TFile *output_plots = new TFile(output,"RECREATE","",0);
  TH1D *histo[allPlots][process_types][jet_cats][3];
  TString plotName; vector<TString> plotName_;
  for(thePlot=0; thePlot<allPlots; thePlot++){
    TH1D* histos = makeHisto(thePlot, plotName);
    plotName_.push_back(plotName);
    for(unsigned int i_jet=0; i_jet<jet_cats; i_jet++) { for(unsigned int i_type=0; i_type<process_types; i_type++) { for(unsigned int i_flav=0; i_flav<3; i_flav++) {
      histo[thePlot][i_type][i_jet][i_flav] = (TH1D*) histos->Clone(Form("%s %s (%s,%s)",plotName.Data(), categoryName_[i_type].Data(), jetString_[i_jet].Data(), lepString_[i_flav].Data()));
    }}}
    histos->Reset(); histos->Clear();
  }

  // set up tree
  TFile *batchFile = TFile::Open(batchTree,"READ");
  TTree *tree_batchPlots     = (TTree*) batchFile->Get("tree_batchPlots"    ),
        *tree_batchPlotIndex = (TTree*) batchFile->Get("tree_batchPlotIndex");
  tree_batchPlots->SetBranchAddress("thePlot"     , &thePlot        );
  tree_batchPlots->SetBranchAddress("nJets"       , &batch_nJot     );
  tree_batchPlots->SetBranchAddress("flavor"      , &batch_flavor   );
  tree_batchPlots->SetBranchAddress("category"    , &batch_category );
  tree_batchPlots->SetBranchAddress("theVar"      , &theVar         );
  tree_batchPlots->SetBranchAddress("totalWeight" , &totalWeight    );

  printf("scanning batch tree...\n");
  Long64_t nentries=tree_batchPlots->GetEntries();
  for(Long64_t i=0; i<nentries; i++) {
    if(i%1000000==0 || i+1==nentries ) printf("event %lld out of %lld\n",i+1,nentries);
    tree_batchPlots->GetEntry(i);
    histo[thePlot][batch_category][batch_nJot][batch_flavor]->Fill(theVar, totalWeight);
    histo[thePlot][batch_category][batch_nJot][0           ]->Fill(theVar, totalWeight);
    histo[thePlot][batch_category][jet_cats-1][batch_flavor]->Fill(theVar, totalWeight);
    histo[thePlot][batch_category][jet_cats-1][0           ]->Fill(theVar, totalWeight);
  }
  batchFile->Close();
  output_plots->cd();
  for(thePlot=0; thePlot<allPlots; thePlot++) {
    TDirectory *plotDir = output_plots->mkdir(plotName_[thePlot]);
    plotDir->cd();
    for(unsigned int i_jet=0; i_jet<jet_cats; i_jet++) {
      TDirectory *jetDir = plotDir->mkdir(jetString_[i_jet]);
      jetDir->cd();
      for(unsigned int i_flav=0; i_flav<3; i_flav++) {
        TDirectory *flavDir = jetDir->mkdir(lepString_[i_flav]);
        flavDir->cd();
        for(unsigned int i_type=0; i_type<process_types; i_type++) {
          histo[thePlot][i_type][i_jet][i_flav] -> Write(categoryName_[i_type]);
        }
      }
    }
  }
  output_plots->Close();
  printf("done!\n");
}


