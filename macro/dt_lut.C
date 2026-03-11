#if defined(__ACLIC__)
#include "../src/PrtTools.h"
#include "../src/PrtLutNode.h" 
#else
R__LOAD_LIBRARY(../build/libPrt.so)
#endif

void dt_lut(TString in="../data/lut2.root"){

  TH1F *hTime[4];
  for (Int_t i = 0; i < 4; i++) {
    hTime[i] = new TH1F(Form("id - %d", i),"", 500, 0, 70);
    hTime[i]->SetLineColor(i + 1);
    hTime[i]->GetXaxis()->SetTitle("length of the photon path in prism [cm]");
    hTime[i]->GetYaxis()->SetTitle("entries [#]");
  }

  in="~/net/lust/data/eic/ev_c5/5/1/lut/lut.avr.root";
  TFile* file = new TFile(in);
  TTree* tree= (TTree *)file->Get("prtlut") ;
  TClonesArray* lut = new TClonesArray("PrtLutNode");
  tree->SetBranchAddress("LUT",&lut); 
  tree->GetEntry(0);
 
  PrtLutNode *node;
  for(Int_t ientry=0; ientry<lut->GetEntriesFast(); ientry++){
    if(ientry%1000==0)std::cout<<"entry  "<<ientry <<std::endl;    
    node = (PrtLutNode*) lut->At(ientry);
    Int_t size = node->Entries(); 
    for(int i=0; i<size; i++){
      //TVector3 dird = node->GetEntry(i);
      Double_t time = node->GetTime(i);
      Int_t nrefl = node->GetNRefl(i);

      hTime[0]->Fill(time*19.8);
      // if(nrefl>=4) continue;
      // hTime[nrefl]->Fill(time*19.8);
    }
  }

  in = "~/net/lust/data/eic/ev_c5/5/3/lut/lut.avr.root";
  file = new TFile(in);
  tree = (TTree *)file->Get("prtlut");
  lut = new TClonesArray("PrtLutNode");
  tree->SetBranchAddress("LUT", &lut);
  tree->GetEntry(0);

  for (Int_t ientry = 0; ientry < lut->GetEntriesFast(); ientry++) {
    if (ientry % 1000 == 0) std::cout << "entry  " << ientry << std::endl;
    node = (PrtLutNode *)lut->At(ientry);
    Int_t size = node->Entries();
    for (int i = 0; i < size; i++) {
      Double_t time = node->GetTime(i);
      Int_t nrefl = node->GetNRefl(i);

      hTime[1]->Fill(time * 19.8);
    }
  }

  // TCanvas* c3 = new TCanvas("c3","c3",900,400);
  // for(Int_t i=0; i<4; i++){
  //   hTime[i]->SetTitle(0);
  //   hTime[i]->GetXaxis()->SetTitleSize(0.05);
  //   hTime[i]->GetYaxis()->SetTitleSize(0.05);
  //   hTime[i]->GetXaxis()->SetTitleOffset(0.85);
  //   hTime[i]->GetYaxis()->SetTitleOffset(0.45);
  //   if(i==0)  hTime[i]->Draw();
  //   else hTime[i]->Draw("same");   
  // }

  PrtTools t;
  gStyle->SetOptStat(0);
  t.add_canvas("lut_length",1200,600);
  hTime[0]->Draw();
  hTime[1]->Draw("same");  

  TLegend *l = new TLegend(0.56, 0.69, 0.86, 0.89);
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(hTime[0], "prism ", "lp");
  l->AddEntry(hTime[1], "rectangular prism", "lp");
  l->Draw();
  gPad->SetGrid();
  t.save_canvas("data/dt_lut", 1);
}
