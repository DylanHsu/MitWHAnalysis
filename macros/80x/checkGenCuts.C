#include "PandaTree/Objects/interface/Event.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include <fstream>

void checkGenCuts(TString cut="", unsigned int num_files=100, bool verbose=false) {
  ifstream ifs;
  ifs.open("pandafiles_Wpt0To50.txt");

  unsigned int i_file=0;
  unsigned int num_files_read=0;
  Long64_t nentries=0, nentries_pass=0, file_nentries, file_nentries_pass;
  for (std::string local_pandafile; std::getline(ifs, local_pandafile) && i_file<num_files; i_file++) {
    file_nentries=0; file_nentries_pass=0;
    string pandafile="root://xrootd.cmsaf.mit.edu/"+local_pandafile;
    if(verbose) printf("Opening %s ... \n", pandafile.c_str());
    TFile *the_file = TFile::Open(pandafile.c_str());
    if(!the_file || !the_file->IsOpen()) {
      if(verbose) printf("File could not be opened\n");
      continue;
    } else num_files_read++;
    TTree *tree=(TTree*)the_file->Get("events");
    panda::Event event; // create an Event object
    event.setStatus(*tree, {"!*"});
    event.setAddress(*tree, {"genParticles","pfMet","genMet",}); 
    file_nentries      = tree->GetEntries();
    for (Long64_t iEntry = 0; iEntry < file_nentries; ++iEntry) {
      event.getEntry(*tree, iEntry);

      if(event.genMet.pt <= 40) continue; //gen MET > 40

      // Find at least one energetic generator electron/muon that has delta phi of pi or greater with the MET.
      int nGoodGenLep=0;
      vector< unsigned > goodGenLeps;
      for (unsigned iG = 0; iG != event.genParticles.size(); ++iG){
        auto& genParticle = event.genParticles[iG];
        //if (genParticle.statusFlags != 1) continue;
        if (genParticle.pt() < 40) continue;
        if (TMath::Abs(genParticle.pdgid) != 11 && TMath::Abs(genParticle.pdgid) != 13) continue;
        //TVector3 genParticle3;
        //genParticle3.SetPtEtaPhi(genParticle.pt(),genParticle.eta(),genParticle.phi());
        if(TMath::Abs(TVector2::Phi_mpi_pi(event.pfMet.phi-genParticle.phi())) <= 3.14159/2.) continue;
        if(TMath::Sqrt(2*genParticle.pt()*event.genMet.pt * (1-TMath::Cos(TMath::Abs(TVector2::Phi_mpi_pi(event.pfMet.phi-genParticle.phi()))))) <= 40) continue; //gen MT > 40;

        nGoodGenLep++; goodGenLeps.push_back(iG);
        
      }
      if(nGoodGenLep==0) continue;
      file_nentries_pass++;
    }
    if(verbose) printf("%lld / %lld passed\n", file_nentries_pass, file_nentries);
    the_file->Close();
    nentries+=file_nentries;
    nentries_pass+=file_nentries_pass;
  }
  printf("%lld / %lld passed the cut \"%s\"\n", nentries_pass, nentries, cut.Data());
  printf("read %d / %d files\n", num_files_read, num_files);
  ifs.close();

}
