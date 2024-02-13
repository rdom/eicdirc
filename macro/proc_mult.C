#if defined(__ACLIC__)
#include "PrtTools.h"
#else
R__LOAD_LIBRARY(../build/libPrt.so)
#endif

void proc_mult(TString in = "../build/hits.root") {

  PrtTools t(in);
  double mom = t.run()->getMomentum();
  double theta = t.run()->getTheta();
  double phi = t.run()->getPhi();
  int pid = t.run()->getPid();

  TH1F *hmult = new TH1F("hmult", "hmult", 250, 0, 250);

  while (t.next() && t.i() < 10000) {
    hmult->Fill(t.event()->getHits().size());
  }

  hmult->Draw();
  double nph = hmult->GetMean();
  std::cout << "nph " << nph << std::endl;
    
  TFile fc(in + "_r.root", "recreate");
  TTree *tc = new TTree("reco", "reco");
  tc->Branch("mom", &mom, "mom/D");
  tc->Branch("theta", &theta, "theta/D");
  tc->Branch("phi", &phi, "phi/D");
  tc->Branch("nph", &nph, "nph/D");
  tc->Branch("pid", &pid, "nph/I");
  tc->Fill();
  tc->Write();
}
