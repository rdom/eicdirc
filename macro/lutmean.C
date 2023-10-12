#if defined(__ACLIC__)
#include "PrtTools.h"
#else
R__LOAD_LIBRARY(../build/libPrt.so)
#endif

void lutmean(TString inFile = "../data/lut.root") {
  TString outFile = inFile.Copy().ReplaceAll(".root", ".avr.root");

  PrtTools t;
  TClonesArray *fLutNew = new TClonesArray("PrtLutNode");
  TTree *fTreeNew = new TTree("prtlut", "Look-up table for DIRC. Averaged");
  fTreeNew->Branch("LUT", &fLutNew, 256000, 2);

  auto run = t.get_run(inFile);
  int nch = run->getNpmt() * run->getNpix();

  TClonesArray &fLutaNew = *fLutNew;
  for (Long64_t n = 0; n < nch; n++) {
    new ((fLutaNew)[n]) PrtLutNode(-1);
  }

  TFile *f = TFile::Open(inFile, "READ");
  TTree *tree = (TTree *)f->Get("prtlut");
  TClonesArray *fLut = new TClonesArray("PrtLutNode");
  tree->SetBranchAddress("LUT", &fLut);
  tree->GetEntry(0);

  TCanvas *c = new TCanvas("c", "c", 0, 0, 800, 1200);
  c->Divide(1, 2);
  TH1F *histNode = new TH1F("LutNode", "Node vs Multiplicity", 30000, 0, 150000);
  TH1F *hTime = new TH1F("hTime", "Time", 5000, 0, 10);
  TH1F *hDir = new TH1F("hDir", "X component", 1000, -1, 1);

  std::vector<TVector3> vArray[100];
  std::vector<Double_t> tArray[100];
  std::vector<Long_t> pArray;
  std::vector<int> rArray;

  TVector3 dir, dir2, sum;
  Double_t angle, minangle, time, sumt;
  Long_t pathid;
  int nrefl;
  PrtLutNode *node;

  for (Int_t inode = 0; inode < fLut->GetEntriesFast(); inode++) {
    if (inode % 1000 == 0) cout << "Node # " << inode << endl;
    node = (PrtLutNode *)fLut->At(inode);

    for (int i = 0; i < node->Entries(); i++) {
      dir = node->GetEntry(i).Unit();
      time = node->GetTime(i);
      pathid = node->GetPathId(i);
      nrefl = node->GetNRefl(i);

      bool newid = true;
      for (int j = 0; j < pArray.size(); j++) {
        if (pathid == pArray[j]) {
          vArray[j].push_back(dir);
          tArray[j].push_back(time);
          newid = false;
        }
      }
      if (newid) {
        vArray[pArray.size()].push_back(dir);
        tArray[pArray.size()].push_back(time);
        pArray.push_back(pathid);
        rArray.push_back(nrefl);
      }
    }

    for (int j = 0; j < pArray.size(); j++) {
      sum = TVector3(0, 0, 0);
      sumt = 0;
      hDir->Reset();
      hTime->Reset();

      for (int v = 0; v < vArray[j].size(); v++) {
        sum += vArray[j][v];
        sumt += tArray[j][v];

        hDir->Fill(vArray[j][v].Y());
        hTime->Fill(tArray[j][v]);
      }

      if (vArray[j].size() < 2) continue;
      if (rArray[j] > 15) continue;

      if (hDir->GetStdDev() > 0.02) {
        std::cout << inode << " " << pArray[j]
                  << " hDir->GetStdDev() ================  " << hDir->GetStdDev() << std::endl;

        c->cd(1);
        hTime->Draw();
        c->cd(2);
        hDir->Draw();
        c->Update();
        c->WaitPrimitive();
      }

      sum *= 1 / (Double_t)vArray[j].size();
      sumt *= 1. / (Double_t)tArray[j].size();

      ((PrtLutNode *)(fLutNew->At(inode)))
        ->AddEntry(inode, sum, pArray[j], rArray[j], sumt, node->GetDigiPos(), node->GetDigiPos());
    }

    for (int i = 0; i < 100; i++) {
      vArray[i].clear();
      tArray[i].clear();
    }
    pArray.clear();
    rArray.clear();
  }

  TFile *fFileNew = TFile::Open(outFile, "RECREATE");
  fTreeNew->Fill();
  fTreeNew->Write();
}
