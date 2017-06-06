#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector2.h>
#include <iostream>
#include <fstream>

//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Mon Jun  5 17:00:40 2017 by ROOT version6.06/01)
//   from TTree events/events
//   found on file: Diboson_wz_0_0.root
//////////////////////////////////////////////////////////
void whAnalysis(
 string subdirectory=""
) {
  // Hardcoded settings
  double mcPrescale = 1; 
  Double_t lumi = 35.9;
  const unsigned int jet_categories=4; // 0-jet, 1-jet, 2-jet, and inclusive
  const unsigned int process_types=8; 
  TString filesPathDA   = "/data/t3home000/dhsu/panda/v_004_0/";
  TString filesPathMC   = "/data/t3home000/dhsu/panda/v_004_0/";
  
  // Set up output dirs
  if(subdirectory!="" && subdirectory.c_str()[0]!='/') subdirectory = "/"+subdirectory;
  system(("mkdir -p MitWHAnalysis/datacards"+subdirectory).c_str());
  system(("mkdir -p MitWHAnalysis/plots"+subdirectory).c_str());

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infileName_, signalName_, categoryName_, jetString_, lepString_;
  vector<Int_t> infileCat_, signalIndex_;  

  
  // Data files
  
  // MC files
  infileName_.push_back(Form("%sTTbar_Powheg.root"   ,  filesPathMC.Data()));  infileCat_.push_back(1);
  infileName_.push_back(Form("%sWJets_EWKWPlus.root" ,  filesPathMC.Data()));  infileCat_.push_back(2);
  infileName_.push_back(Form("%sWJets_EWKWMinus.root",  filesPathMC.Data()));  infileCat_.push_back(2);
  infileName_.push_back(Form("%sWJets_nlo.root"      ,  filesPathMC.Data()));  infileCat_.push_back(3);
  infileName_.push_back(Form("%sZJets_nlo.root"      ,  filesPathMC.Data()));  infileCat_.push_back(4);
  infileName_.push_back(Form("%sDiboson_wz.root"     ,  filesPathMC.Data()));  infileCat_.push_back(5);
  infileName_.push_back(Form("%sDiboson_zz.root"     ,  filesPathMC.Data()));  infileCat_.push_back(6);
  //infileName_.push_back(Form("%sDiboson_ww.root"     ,  filesPathMC.Data()));  infileCat_.push_back(7);

  // Process types
  categoryName_.push_back("Data");
  categoryName_.push_back("Top");
  categoryName_.push_back("EWK W+jets");
  categoryName_.push_back("Cont. W+jets");
  categoryName_.push_back("Cont. Z+jets");
  categoryName_.push_back("WZ");
  categoryName_.push_back("ZZ");
  categoryName_.push_back("WW");

  jetString_.push_back("0j");
  jetString_.push_back("1j");
  jetString_.push_back("2j");
  jetString_.push_back("nj");

  lepString_.push_back("e");
  lepString_.push_back("m");
  lepString_.push_back("l");


  
  //Declaration of leaves types
  Int_t nJet, nJot, nLooseLep, nLooseElectron, nLooseMuon, nTightLep, nTightElectron, nTightMuon, looseLep1PdgId, looseLep2PdgId, looseLep1IsTight, looseLep2IsTight, nTau, jot1VBFID, jetNMBtags, nB, isGS, looseLep1IsHLTSafe, looseLep2IsHLTSafe, runNumber, lumiNumber, npv, pu, trigger, metFilter, egmFilter, genTopIsHad, genAntiTopIsHad, nIsoJet, jet1Flav, jet1IsTight, jet2Flav, isojet1Flav, isojet2Flav, jetNBtags, isojetNBtags, nFatjet, nHF, nLoosePhoton, nTightPhoton, loosePho1IsTight;
  Float_t jot1Phi, jot1Pt, jot1GenPt, jot1Eta, jot2Phi, jot2Pt, jot2GenPt, jot2Eta, jot12DPhi, jot12Mass, jot12DEta, pfmetUp, pfmetDown, pfUWmagUp, pfUZmagUp, pfUAmagUp, pfUmagUp, pfUWmagDown, pfUZmagDown, pfUAmagDown, pfUmagDown, jot1EtaUp, jot1EtaDown, jot1PtUp, jot1PtDown, jot2PtUp, jot2PtDown, jot12MassUp, jot12DEtaUp, jot12DPhiUp, jot12MassDown, jot12DEtaDown, jot12DPhiDown, jot2EtaUp, jot2EtaDown, scale[6], sf_btag0, sf_btag1, sf_btag2, sf_btagGT0, sf_sjbtag0, sf_sjbtag1, sf_sjbtag2, sf_sjbtagGT0, sf_btag0BUp, sf_btag1BUp, sf_btag2BUp, sf_btagGT0BUp, sf_sjbtag0BUp, sf_sjbtag1BUp, sf_sjbtag2BUp, sf_sjbtagGT0BUp, sf_btag0BDown, sf_btag1BDown, sf_btag2BDown, sf_btagGT0BDown, sf_sjbtag0BDown, sf_sjbtag1BDown, sf_sjbtag2BDown, sf_sjbtagGT0BDown, sf_btag0MUp, sf_btag1MUp, sf_btag2MUp, sf_btagGT0MUp, sf_sjbtag0MUp, sf_sjbtag1MUp, sf_sjbtag2MUp, sf_sjbtagGT0MUp, sf_btag0MDown, sf_btag1MDown, sf_btag2MDown, sf_btagGT0MDown, sf_sjbtag0MDown, sf_sjbtag1MDown, sf_sjbtag2MDown, sf_sjbtagGT0MDown, sf_metTrigZmm, sf_qcdV_VBF, pfmetRaw, mcWeight, normalizedWeight, filter_maxRecoil, filter_whichRecoil, sf_ewkV, sf_qcdV, sf_ewkV2j, sf_qcdV2j, sf_qcdTT, sf_lepID, sf_lepIso, sf_lepTrack, sf_pho, sf_eleTrig, sf_phoTrig, sf_metTrig, sf_pu, sf_npv, sf_tt, sf_tt_ext, sf_tt_bound, sf_tt8TeV, sf_tt8TeV_ext, sf_tt8TeV_bound, sf_phoPurity, pfmet, pfmetphi, pfmetnomu, puppimet, puppimetphi, calomet, calometphi, pfcalobalance, sumET, trkmet, puppiUWmag, puppiUWphi, puppiUZmag, puppiUZphi, puppiUAmag, puppiUAphi, puppiUperp, puppiUpara, puppiUmag, puppiUphi, pfUWmag, pfUWphi, pfUZmag, pfUZphi, pfUAmag, pfUAphi, pfUperp, pfUpara, pfUmag, pfUphi, dphipfmet, dphipuppimet, dphipuppiUW, dphipuppiUZ, dphipuppiUA, dphipfUW, dphipfUZ, dphipfUA, dphipuppiU, dphipfU, trueGenBosonPt, genBosonPt, genBosonEta, genBosonMass, genBosonPhi, genWPlusPt, genWMinusPt, genWPlusEta, genWMinusEta, genTopPt, genTopEta, genAntiTopPt, genAntiTopEta, genTTPt, genTTEta, jet1Phi, jet1Pt, jet1GenPt, jet1Eta, jet1CSV, jet2Phi, jet2Pt, jet2GenPt, jet2Eta, jet2CSV, isojet1Pt, isojet1CSV, isojet2Pt, isojet2CSV, loosePho1Pt, loosePho1Eta, loosePho1Phi, looseLep1Pt, looseLep1Eta, looseLep1Phi, looseLep2Pt, looseLep2Eta, looseLep2Phi, diLepMass, mT, scaleUp, scaleDown, pdfUp, pdfDown;
  ULong64_t eventNumber;
  
  // Set up histograms
  char output[200];
  sprintf(output,"MitWHAnalysis/plots%s/histo_wh_nice.root",subdirectory.c_str());
  TFile *output_plots = new TFile(output,"RECREATE");
  TH1D *histo_mT_presel[jet_categories][process_types][3];
  for(unsigned int i_jet=0; i_jet<jet_categories; i_jet++) { for(unsigned int i_type=0; i_type<process_types; i_type++) { for(unsigned int i_flav=0; i_flav<3; i_flav++) {
      histo_mT_presel[i_jet][i_type][i_flav] = new TH1D(
        Form("histo_mT_presel_%s_%s_%d", jetString_[i_jet].Data(), lepString_[i_flav].Data(), i_type),
        Form("preselection m_{T} for %s (%s channel, %s)", categoryName_[i_type].Data(), lepString_[i_flav].Data(), jetString_[i_jet].Data()),
      200, 0, 1000);

  }}}
  
  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infileName_.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infileName_[ifile].Data());

    TFile *the_input_file = TFile::Open(infileName_[ifile].Data());
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
   
    // Set branch addresses.
    the_input_tree->SetBranchAddress("nJet",&nJet);
    //the_input_tree->SetBranchAddress("nJot",&nJot);
    //the_input_tree->SetBranchAddress("jot1Phi",&jot1Phi);
    //the_input_tree->SetBranchAddress("jot1Pt",&jot1Pt);
    //the_input_tree->SetBranchAddress("jot1GenPt",&jot1GenPt);
    //the_input_tree->SetBranchAddress("jot1Eta",&jot1Eta);
    //the_input_tree->SetBranchAddress("jot2Phi",&jot2Phi);
    //the_input_tree->SetBranchAddress("jot2Pt",&jot2Pt);
    //the_input_tree->SetBranchAddress("jot2GenPt",&jot2GenPt);
    //the_input_tree->SetBranchAddress("jot2Eta",&jot2Eta);
    //the_input_tree->SetBranchAddress("jot12DPhi",&jot12DPhi);
    //the_input_tree->SetBranchAddress("jot12Mass",&jot12Mass);
    //the_input_tree->SetBranchAddress("jot12DEta",&jot12DEta);
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
    the_input_tree->SetBranchAddress("mcWeight",&mcWeight);
    the_input_tree->SetBranchAddress("normalizedWeight",&normalizedWeight);
    //the_input_tree->SetBranchAddress("trigger",&trigger);
    //the_input_tree->SetBranchAddress("metFilter",&metFilter);
    //the_input_tree->SetBranchAddress("egmFilter",&egmFilter);
    //the_input_tree->SetBranchAddress("filter_maxRecoil",&filter_maxRecoil);
    //the_input_tree->SetBranchAddress("filter_whichRecoil",&filter_whichRecoil);
    //the_input_tree->SetBranchAddress("sf_ewkV",&sf_ewkV);
    //the_input_tree->SetBranchAddress("sf_qcdV",&sf_qcdV);
    //the_input_tree->SetBranchAddress("sf_ewkV2j",&sf_ewkV2j);
    //the_input_tree->SetBranchAddress("sf_qcdV2j",&sf_qcdV2j);
    //the_input_tree->SetBranchAddress("sf_qcdTT",&sf_qcdTT);
    //the_input_tree->SetBranchAddress("sf_lepID",&sf_lepID);
    //the_input_tree->SetBranchAddress("sf_lepIso",&sf_lepIso);
    //the_input_tree->SetBranchAddress("sf_lepTrack",&sf_lepTrack);
    //the_input_tree->SetBranchAddress("sf_pho",&sf_pho);
    //the_input_tree->SetBranchAddress("sf_eleTrig",&sf_eleTrig);
    //the_input_tree->SetBranchAddress("sf_phoTrig",&sf_phoTrig);
    //the_input_tree->SetBranchAddress("sf_metTrig",&sf_metTrig);
    //the_input_tree->SetBranchAddress("sf_pu",&sf_pu);
    //the_input_tree->SetBranchAddress("sf_npv",&sf_npv);
    //the_input_tree->SetBranchAddress("sf_tt",&sf_tt);
    //the_input_tree->SetBranchAddress("sf_tt_ext",&sf_tt_ext);
    //the_input_tree->SetBranchAddress("sf_tt_bound",&sf_tt_bound);
    //the_input_tree->SetBranchAddress("sf_tt8TeV",&sf_tt8TeV);
    //the_input_tree->SetBranchAddress("sf_tt8TeV_ext",&sf_tt8TeV_ext);
    //the_input_tree->SetBranchAddress("sf_tt8TeV_bound",&sf_tt8TeV_bound);
    //the_input_tree->SetBranchAddress("sf_phoPurity",&sf_phoPurity);
    the_input_tree->SetBranchAddress("pfmet",&pfmet);
    the_input_tree->SetBranchAddress("pfmetphi",&pfmetphi);
    //the_input_tree->SetBranchAddress("pfmetnomu",&pfmetnomu);
    //the_input_tree->SetBranchAddress("puppimet",&puppimet);
    //the_input_tree->SetBranchAddress("puppimetphi",&puppimetphi);
    the_input_tree->SetBranchAddress("calomet",&calomet);
    the_input_tree->SetBranchAddress("calometphi",&calometphi);
    //the_input_tree->SetBranchAddress("pfcalobalance",&pfcalobalance);
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
    //the_input_tree->SetBranchAddress("pfUmag",&pfUmag);
    //the_input_tree->SetBranchAddress("pfUphi",&pfUphi);
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
    //the_input_tree->SetBranchAddress("jetNBtags",&jetNBtags);
    //the_input_tree->SetBranchAddress("isojetNBtags",&isojetNBtags);
    //the_input_tree->SetBranchAddress("nFatjet",&nFatjet);
    the_input_tree->SetBranchAddress("nHF",&nHF);
    //the_input_tree->SetBranchAddress("nLoosePhoton",&nLoosePhoton);
    //the_input_tree->SetBranchAddress("nTightPhoton",&nTightPhoton);
    the_input_tree->SetBranchAddress("loosePho1IsTight",&loosePho1IsTight);
    the_input_tree->SetBranchAddress("loosePho1Pt",&loosePho1Pt);
    the_input_tree->SetBranchAddress("loosePho1Eta",&loosePho1Eta);
    the_input_tree->SetBranchAddress("loosePho1Phi",&loosePho1Phi);
    the_input_tree->SetBranchAddress("nLooseLep",&nLooseLep);
    the_input_tree->SetBranchAddress("nLooseElectron",&nLooseElectron);
    the_input_tree->SetBranchAddress("nLooseMuon",&nLooseMuon);
    the_input_tree->SetBranchAddress("nTightLep",&nTightLep);
    the_input_tree->SetBranchAddress("nTightElectron",&nTightElectron);
    the_input_tree->SetBranchAddress("nTightMuon",&nTightMuon);
    the_input_tree->SetBranchAddress("looseLep1PdgId",&looseLep1PdgId);
    the_input_tree->SetBranchAddress("looseLep2PdgId",&looseLep2PdgId);
    the_input_tree->SetBranchAddress("looseLep1IsTight",&looseLep1IsTight);
    the_input_tree->SetBranchAddress("looseLep2IsTight",&looseLep2IsTight);
    the_input_tree->SetBranchAddress("looseLep1Pt",&looseLep1Pt);
    the_input_tree->SetBranchAddress("looseLep1Eta",&looseLep1Eta);
    the_input_tree->SetBranchAddress("looseLep1Phi",&looseLep1Phi);
    the_input_tree->SetBranchAddress("looseLep2Pt",&looseLep2Pt);
    the_input_tree->SetBranchAddress("looseLep2Eta",&looseLep2Eta);
    the_input_tree->SetBranchAddress("looseLep2Phi",&looseLep2Phi);
    //the_input_tree->SetBranchAddress("diLepMass",&diLepMass);
    //the_input_tree->SetBranchAddress("nTau",&nTau);
    the_input_tree->SetBranchAddress("mT",&mT);
    the_input_tree->SetBranchAddress("scaleUp",&scaleUp);
    the_input_tree->SetBranchAddress("scaleDown",&scaleDown);
    the_input_tree->SetBranchAddress("pdfUp",&pdfUp);
    the_input_tree->SetBranchAddress("pdfDown",&pdfDown);

    //     This is the loop skeleton
    //       To read only selected branches, Insert statements like:
    // events->SetBranchStatus("*",0);  // disable all branches
    // TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

    Long64_t nentries = the_input_tree->GetEntries();

    Long64_t nbytes = 0;
    double theMCPrescale = mcPrescale;
    if(infileCat_[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(nentries/theMCPrescale); ++i) {
      if(i%100000==0 || i+1==int(nentries/theMCPrescale) ) printf("event %d out of %d\n",i+1,int(nentries/theMCPrescale));
      the_input_tree->GetEntry(i);
      bool passPresel=false;
      int flavor=-1;
      double ptFrac = looseLep1Pt/pfmet;
      
      if(nLooseLep<1) continue;
      if     (nTightElectron ==1) flavor=0;
      else if(nTightMuon     ==1) flavor=1;
      TVector2 vLep1, vPFMET; 
      vLep1.SetMagPhi(looseLep1Pt, looseLep1Phi);
      vPFMET.SetMagPhi(pfmet, pfmetphi);
      double deltaPhiLepMET = vLep1.DeltaPhi(vPFMET);
      
      unsigned int nNiceJets=0;
      if(jet1Pt > 30) nNiceJets++;
      if(jet2Pt > 30) nNiceJets++;

      if(
        nLooseLep==1 && looseLep1IsTight && (flavor>=0) &&
        looseLep1Pt>50 &&
        pfmet>100 &&
        ptFrac > 0.4 && ptFrac < 1.5 &&
        TMath::Abs(deltaPhiLepMET) > 2 &&
        jet1CSV<0.5426 && jet2CSV<0.5426
      ) passPresel=true;
      
      //begin event weighting
      double totalWeight = normalizedWeight * 1000. * lumi * theMCPrescale;
      if(passPresel && nNiceJets < (jet_categories-1) ) histo_mT_presel[nNiceJets][infileCat_[ifile]][flavor]->Fill(mT, mcWeight);
      if(passPresel) {
        histo_mT_presel[jet_categories-1][infileCat_[ifile]][flavor]->Fill(mT, totalWeight);
        histo_mT_presel[jet_categories-1][infileCat_[ifile]][2]     ->Fill(mT, totalWeight);
      }

    }
    the_input_file->Close();
  }
  output_plots->cd();
  for(unsigned int i_jet=0; i_jet<jet_categories; i_jet++) { for(unsigned int i_type=0; i_type<process_types; i_type++) { for(unsigned int i_flav=0; i_flav<3; i_flav++) {
      if(i_flav<2) continue;
      histo_mT_presel[i_jet][i_type][i_flav] -> Write();
  }}}
  output_plots->Close();

}
