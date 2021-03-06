#include <TMVA/Config.h>
#include <TMVA/Factory.h>
#include <TMVA/Types.h>
#include <TFile.h>
#include <TCut.h>
#include <TTree.h>
#include <TString.h>
#include <TVector2.h>

void whMVA(unsigned int mode=1, string extra_string="", unsigned int njets=0) {
  // Use TMVA to analyze/train on the trees generated by PandaAnalyzer
  
  TString filesPathDA   = "/data/t3home000/dhsu/panda/merged_skims/";
  TString filesPathMC   = "/data/t3home000/dhsu/panda/merged_skims/";
  
  TFile *output_file;
  TMVA::Factory *factory;
  TMVA::gConfig().GetVariablePlotting().fMaxNumOfAllowedVariablesForScatterPlots=50;
  
  // Initialize the factory
  system("mkdir -p MitWHAnalysis/mva");
  if(mode==1) {
    output_file=TFile::Open(("MitWHAnalysis/mva/training_result_BDT_WHinv"+(extra_string == "" ? "" : "_"+extra_string)+".root").c_str(), "RECREATE");
    factory = new TMVA::Factory("bdt", output_file, "V:!Silent:DrawProgressBar:AnalysisType=Multiclass");
  } else if(mode==2) {
    output_file=TFile::Open(("MitWHAnalysis/mva/training_result_BDT_WHinv"+(extra_string == "" ? "" : "_"+extra_string)+".root").c_str(), "RECREATE");
    factory = new TMVA::Factory("bdt", output_file, "V:!Silent:DrawProgressBar:AnalysisType=Multiclass");
  } else return;

  //if(mode==1) {
  //  output_file=TFile::Open(("MitZHAnalysis/mva/training_result_RCO_massIndependent_"+mediator+".root").c_str(), "RECREATE");
  //  factory = new TMVA::Factory("rco", output_file, "AnalysisType=Classification");
  //} else if(mode==2) {
  //  output_file=TFile::Open(("MitZHAnalysis/mva/training_result_BDT_massIndependent_"+mediator+".root").c_str(), "RECREATE");
  //  factory = new TMVA::Factory("bdt", output_file, "!V:!Silent:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
  //} else if(mode==3) {
  //  output_file=TFile::Open(("MitZHAnalysis/mva/training_result_BDT_massDependent_"+signal_model+".root").c_str(), "RECREATE");
  //  factory = new TMVA::Factory("bdt", output_file, "!V:!Silent:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
  //} else if(mode==4) {
  //  output_file=TFile::Open(("MitZHAnalysis/mva/training_result_BDT_multiClass_"+signal_model+(extra_string == "" ? "" : "_"+extra_string)+".root").c_str(), "RECREATE");
  //  factory = new TMVA::Factory("bdt", output_file, "!V:!Silent:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Multiclass");
  //} else { assert(0); return; }

  // Determine the input trees
  vector<TString> infileName_, categoryName_;
  vector<Int_t> infileCat_;  
  //categoryName_.push_back("Signal");
  //categoryName_.push_back("Top");
  //categoryName_.push_back("EWK W+jets");
  //categoryName_.push_back("Cont. W+jets");
  //categoryName_.push_back("Z+jets");
  //categoryName_.push_back("WZ");
  //categoryName_.push_back("ZZ");
  //categoryName_.push_back("WW");
  //categoryName_.push_back("QCD, y+jets");
  categoryName_.push_back("Signal");
  categoryName_.push_back("W+jets");
  categoryName_.push_back("Z+jets");
  categoryName_.push_back("VV");
  categoryName_.push_back("Top");
  
  infileName_.push_back(Form("%sWminusH_HToInvisible_WToLNu_M125_13TeV_powheg_pythia8.root", filesPathMC.Data())); infileCat_.push_back(0);
  infileName_.push_back(Form("%sWplusH_HToInvisible_WToLNu_M125_13TeV_powheg_pythia8.root" , filesPathMC.Data())); infileCat_.push_back(0);
  
  infileName_.push_back(Form("%sWJets_EWKWMinus.root"                                      , filesPathMC.Data())); infileCat_.push_back(1);
  infileName_.push_back(Form("%sWJets_EWKWPlus.root"                                       , filesPathMC.Data())); infileCat_.push_back(1);
  
  infileName_.push_back(Form("%sWJets_Wpt0to50.root"                                       , filesPathMC.Data())); infileCat_.push_back(1);
  infileName_.push_back(Form("%sWJets_Wpt50to100.root"                                     , filesPathMC.Data())); infileCat_.push_back(1);
  infileName_.push_back(Form("%sWJets_pt100to250.root"                                     , filesPathMC.Data())); infileCat_.push_back(1);
  infileName_.push_back(Form("%sWJets_pt250to400.root"                                     , filesPathMC.Data())); infileCat_.push_back(1);
  infileName_.push_back(Form("%sWJets_pt400to600.root"                                     , filesPathMC.Data())); infileCat_.push_back(1);
  infileName_.push_back(Form("%sWJets_pt600toinf.root"                                     , filesPathMC.Data())); infileCat_.push_back(1);

  //infileName_.push_back(Form("%sZJets_pt50to100.root"                                      , filesPathMC.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sZJets_pt100to250.root"                                     , filesPathMC.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sZJets_pt250to400.root"                                     , filesPathMC.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sZJets_pt400to650.root"                                     , filesPathMC.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sZJets_pt650toinf.root"                                     , filesPathMC.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sZJets_EWK.root"                                            , filesPathMC.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sZtoNuNu_EWK.root"                                          , filesPathMC.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sZtoNuNu_Zpt50to100.root"                                   , filesPathMC.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sZtoNuNu_Zpt100to250.root"                                  , filesPathMC.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sZtoNuNu_Zpt250to400.root"                                  , filesPathMC.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sZtoNuNu_Zpt400to650.root"                                  , filesPathMC.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sZtoNuNu_Zpt650toinf.root"                                  , filesPathMC.Data())); infileCat_.push_back(2);

  infileName_.push_back(Form("%sDiboson_wz.root"                                           , filesPathMC.Data())); infileCat_.push_back(3);
  infileName_.push_back(Form("%sDiboson_zz.root"                                           , filesPathMC.Data())); infileCat_.push_back(3);
  infileName_.push_back(Form("%sDiboson_ww.root"                                           , filesPathMC.Data())); infileCat_.push_back(3);
  
  infileName_.push_back(Form("%sTTbar_Powheg.root"                                         , filesPathMC.Data())); infileCat_.push_back(4);
  infileName_.push_back(Form("%sSingleTop_tG.root"                                         , filesPathMC.Data())); infileCat_.push_back(4);
  infileName_.push_back(Form("%sSingleTop_tT.root"                                         , filesPathMC.Data())); infileCat_.push_back(4);
  infileName_.push_back(Form("%sSingleTop_tTbar.root"                                      , filesPathMC.Data())); infileCat_.push_back(4);
  infileName_.push_back(Form("%sSingleTop_tW.root"                                         , filesPathMC.Data())); infileCat_.push_back(4);
  infileName_.push_back(Form("%sSingleTop_tbarW.root"                                      , filesPathMC.Data())); infileCat_.push_back(4);
  
  //infileName_.push_back(Form("%sQCD.root"                                                  , filesPathMC.Data())); infileCat_.push_back(8);
  //infileName_.push_back(Form("%sGJets.root"                                                , filesPathMC.Data())); infileCat_.push_back(8);

  TTree *the_trees[64]; TFile *the_files[64]; Int_t iCat=-1;
  for(unsigned int ifile=0; ifile<infileName_.size(); ifile++) {
    the_files[ifile]=TFile::Open(infileName_[ifile], "READ");
    if(!the_files[ifile]) { printf("problem with file %s\n", infileName_[ifile].Data()); assert(0); return; }
    the_trees[ifile] = (TTree*) the_files[ifile]->Get("events");
    if(!the_trees[ifile]) { printf("problem with file %s\n", infileName_[ifile].Data()); assert(0); return; }
    //printf("calling factory->AddTree(the_trees[%d], \"%s\")\n", ifile, categoryName_[infileCat_[ifile]].Data());
    factory->AddTree(the_trees[ifile], categoryName_[infileCat_[ifile]].Data());
    // if(infileCat_[ifile] > iCat) {
    //   iCat=infileCat_[ifile];
    //   factory->SetWeightExpression("normalizedWeight * 35900. * sf_npv * sf_lepID * sf_lepIso * sf_lepTrack", categoryName_[infileCat_[ifile]].Data());
    // }
    factory->SetWeightExpression("normalizedWeight * 35900. * sf_npv * sf_lepID * sf_lepIso * sf_lepTrack", categoryName_[infileCat_[ifile]].Data());
  }
  //factory->SetWeightExpression("normalizedWeight * 35900. * sf_npv * sf_lepID * sf_lepIso * sf_lepTrack");
  if(mode==1) {
    factory->AddVariable( "pfmet"                  , "|E_{T}^{miss}|"                 , ""  , 'F');
    factory->AddVariable( "pfmetphi"               , "E_{T}^{miss} #phi"              , ""  , 'F');
    factory->AddVariable( "mT"                     , "m_{T}(l1,E_{T}^{miss})"         , ""  , 'F');
    factory->AddVariable( "looseLep1Pt"            , "p_{T} l1"                       , ""  , 'F');
    factory->AddVariable( "looseLep1Eta"           , "#eta l1"                        , ""  , 'F');
    factory->AddVariable( "looseLep1Phi"           , "#phi l1"                        , ""  , 'F');
    factory->AddVariable( "pfUmag"                 , "|U_{T}|"                        , ""  , 'F');
    factory->AddVariable( "pfUphi"                 , "U_{T} #phi"                     , ""  , 'F');
    factory->AddVariable( "nAwayPFCH"              , "N PFCH away"                    , ""  , 'I');
    factory->AddVariable( "nTransversePFCH"        , "N PFCH transverse"              , ""  , 'I');
    factory->AddVariable( "nTowardPFCH"            , "N PFCH toward"                  , ""  , 'I');
    factory->AddVariable( "sumPtAwayPFCH"          , "sum p_{T} PFCH away"            , ""  , 'F');
    factory->AddVariable( "sumPtTransversePFCH"    , "sum p_{T} PFCH transverse"      , ""  , 'F');
    factory->AddVariable( "sumPtTowardPFCH"        , "sum p_{T} PFCH toward"          , ""  , 'F');
    factory->AddVariable( "PFCH1Pt"                , "p_{T} PFCH1"                    , ""  , 'F');
    factory->AddVariable( "PFCH2Pt"                , "p_{T} PFCH2"                    , ""  , 'F');
    factory->AddVariable( "PFCH3Pt"                , "p_{T} PFCH3"                    , ""  , 'F');
    factory->AddVariable( "PFCH1Eta"               , "#eta PFCH1"                     , ""  , 'F');
    factory->AddVariable( "PFCH2Eta"               , "#eta PFCH2"                     , ""  , 'F');
    factory->AddVariable( "PFCH3Eta"               , "#eta PFCH3"                     , ""  , 'F');
    factory->AddVariable( "PFCH1Phi"               , "#phi PFCH1"                     , ""  , 'F');
    factory->AddVariable( "PFCH2Phi"               , "#phi PFCH2"                     , ""  , 'F');
    factory->AddVariable( "PFCH3Phi"               , "#phi PFCH3"                     , ""  , 'F');
    //TCut preselectionCut = "nJot==0";
    //TCut preselectionCut = "nJot==0 && looseLep1Pt/pfmet>0.4 && looseLep1Pt/pfmet<1.5 && TMath::Abs(TVector2::Phi_mpi_pi(pfmetphi-looseLep1Phi))>2"; // presel2
    TCut preselectionCut = "nJot==0 && looseLep1Pt/pfmet>0.4 && looseLep1Pt/pfmet<1.5 && TMath::Abs(TVector2::Phi_mpi_pi(pfmetphi-looseLep1Phi))>2.5 && pfmet>100"; // sigsel 1
    //TCut preselectionCut = "nJot==0 && looseLep1Pt/pfmet>0.7 && looseLep1Pt/pfmet<1.2 && TMath::Abs(TVector2::Phi_mpi_pi(pfmetphi-looseLep1Phi))>2.8 && pfmet>100"; // sigsel 2
    //TCut preselectionCut = "nJot==0 && looseLep1Pt/pfmet>0.7 && looseLep1Pt/pfmet<1.2 && TMath::Abs(TVector2::Phi_mpi_pi(pfmetphi-looseLep1Phi))>2.8 && pfmet>200"; // sigsel 3
    factory->PrepareTrainingAndTestTree(preselectionCut, "");
    //factory->BookMethod( TMVA::Types::kBDT, "BDT_WHinv"+(extra_string == "" ? "" : "_"+extra_string), "!H:!V:NTrees=500:MinNodeSize=5%:MaxDepth=5:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=10000:PruneMethod=NoPruning");
    factory->BookMethod( TMVA::Types::kBDT, "BDT_WHinv"+(extra_string == "" ? "" : "_"+extra_string), "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=4:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=10000:PruneMethod=NoPruning");
  } else if(mode==2) {
    factory->AddVariable( "pfmet"                                                    , "|E_{T}^{miss}|"                      , ""  , 'F');
    factory->AddVariable( "mT"                                                       , "m_{T}(l1,E_{T}^{miss})"              , ""  , 'F');
    factory->AddVariable( "looseLep1Pt"                                              , "p_{T} l1"                            , ""  , 'F');
    factory->AddVariable( "looseLep1Eta"                                             , "#eta l1"                             , ""  , 'F');
    factory->AddVariable( "pfUmag"                                                   , "|U_{T}|"                             , ""  , 'F');
    factory->AddVariable( "TMath::Abs(TVector2::Phi_mpi_pi(pfmetphi-looseLep1Phi))"  , "#Delta#phi(E_{T}^{miss}, p_{T}^{l1})", ""  , 'F');
    factory->AddVariable( "TMath::Abs(TVector2::Phi_mpi_pi(pfUphi-looseLep1Phi))"    , "#Delta#phi(U_{T}, p_{T}^{l1})"       , ""  , 'F');
    factory->AddVariable( "nAwayPFCH"                                                , "N PFCH away"                         , ""  , 'I');
    factory->AddVariable( "nTransversePFCH"                                          , "N PFCH transverse"                   , ""  , 'I');
    factory->AddVariable( "nTowardPFCH"                                              , "N PFCH toward"                       , ""  , 'I');
    factory->AddVariable( "sumPtAwayPFCH"                                            , "sum p_{T} PFCH away"                 , ""  , 'F');
    factory->AddVariable( "sumPtTransversePFCH"                                      , "sum p_{T} PFCH transverse"           , ""  , 'F');
    factory->AddVariable( "sumPtTowardPFCH"                                          , "sum p_{T} PFCH toward"               , ""  , 'F');
    factory->AddVariable( "PFCH1Pt"                                                  , "p_{T} PFCH1"                         , ""  , 'F');
    factory->AddVariable( "PFCH2Pt"                                                  , "p_{T} PFCH2"                         , ""  , 'F');
    factory->AddVariable( "PFCH3Pt"                                                  , "p_{T} PFCH3"                         , ""  , 'F');
    factory->AddVariable( "PFCH1Eta"                                                 , "#eta PFCH1"                          , ""  , 'F');
    factory->AddVariable( "PFCH2Eta"                                                 , "#eta PFCH2"                          , ""  , 'F');
    factory->AddVariable( "PFCH3Eta"                                                 , "#eta PFCH3"                          , ""  , 'F');
    factory->AddVariable( "PFCH1Phi"                                                 , "#phi PFCH1"                          , ""  , 'F');
    factory->AddVariable( "PFCH2Phi"                                                 , "#phi PFCH2"                          , ""  , 'F');
    factory->AddVariable( "PFCH3Phi"                                                 , "#phi PFCH3"                          , ""  , 'F');
    if(njets>=1) {
      factory->AddVariable( "jot1Pt"                                                                                    , "p_{T} jet1"                          , "", 'F');
      factory->AddVariable( "jot1Eta"                                                                                   , "#eta jet1"                           , "", 'F');
      factory->AddVariable( "TMath::Abs(TVector2::Phi_mpi_pi(pfmetphi-jot1Phi))"                                        , "#Delta#phi(E_{T}^{miss}, p_{T}^{j1})", "", 'F');
      factory->AddVariable( "TMath::Sqrt(pow(TVector2::Phi_mpi_pi(looseLep1Phi-jot1Phi),2)+pow(looseLep1Eta-jot1Eta,2))", "#DeltaR(p^{l1}, p^{j1})"             , "", 'F');
    }
    if(njets==2) {
      factory->AddVariable( "jot2Pt"                                                                                    , "p_{T} jet2"                          , "", 'F');
      factory->AddVariable( "jot2Eta"                                                                                   , "#eta jet2"                           , "", 'F');
      factory->AddVariable( "TMath::Abs(TVector2::Phi_mpi_pi(pfmetphi-jot2Phi))"                                        , "#Delta#phi(E_{T}^{miss}, p_{T}^{j2})", "", 'F');
      factory->AddVariable( "TMath::Sqrt(pow(TVector2::Phi_mpi_pi(looseLep1Phi-jot2Phi),2)+pow(looseLep1Eta-jot2Eta,2))", "#DeltaR(p^{l1}, p^{j2})"             , "", 'F');
      factory->AddVariable( "TMath::Sqrt(pow(TVector2::Phi_mpi_pi(jot1Phi-jot2Phi),2)+pow(jot1Eta-jot2Eta,2))"          , "#DeltaR(p^{j1}, p^{j2})"             , "", 'F');
    }
    char cutstring[1024];
    sprintf(cutstring,"nJot==%d && looseLep1Pt/pfmet>0.4 && looseLep1Pt/pfmet<1.5 && TMath::Abs(TVector2::Phi_mpi_pi(pfmetphi-looseLep1Phi))>2.5 && pfmet>100 && jetNBtags==0", njets);
    TCut preselectionCut=cutstring;
    factory->PrepareTrainingAndTestTree(preselectionCut, "");
    //factory->BookMethod( TMVA::Types::kBDT, "BDT_WHinv"+(extra_string == "" ? "" : "_"+extra_string), "!H:!V:NTrees=500:MinNodeSize=5%:MaxDepth=5:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=10000:PruneMethod=NoPruning");
    factory->BookMethod( TMVA::Types::kBDT, "BDT_WHinv"+(extra_string == "" ? "" : "_"+extra_string), "!H:!V:NTrees=500:MinNodeSize=5%:MaxDepth=6:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=10000:PruneMethod=NoPruning");
  }

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
}
