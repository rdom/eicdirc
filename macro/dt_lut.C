#include "../../prttools/prttools.C"

#include "../src/PrtLutNode.h" 

void dt_lut(TString inFile="../data/lut2.root"){
  TFile* file = new TFile(inFile);
  TTree* tree= (TTree *)file->Get("prtlut") ;
  TClonesArray* lut = new TClonesArray("PrtLutNode");
  tree->SetBranchAddress("LUT",&lut); 
  tree->GetEntry(0);

  TH1F *hTime[4];
  for(Int_t i=0; i<4; i++){
    hTime[i] = new TH1F(Form("id - %d", i),Form("id%d", i),500,30,70);
    hTime[i]->SetLineColor(i+1);
  }
  PrtLutNode *node;
  std::cout<<"lut->GetEntriesFast()  "<<lut->GetEntriesFast() <<std::endl;
  
  for(Int_t ientry=0; ientry<lut->GetEntriesFast(); ientry++){
    if(ientry%1000==0)std::cout<<"entry  "<<ientry <<std::endl;    
    node = (PrtLutNode*) lut->At(ientry);
    Int_t size = node->Entries(); 
    for(int i=0; i<size; i++){
      //TVector3 dird = node->GetEntry(i);
      Double_t time = node->GetTime(i);
      Int_t nrefl = node->GetNRefl(i);
      if(nrefl>=4) continue;
      hTime[nrefl]->Fill(time*19.8);
    }
  }
  std::cout<<"End " <<std::endl;
  
  
  TCanvas* c3 = new TCanvas("c3","c3",900,400);
  for(Int_t i=0; i<4; i++){
    hTime[i]->SetTitle(0);
    hTime[i]->GetXaxis()->SetTitle("length of the photon path in prism [cm]");
    hTime[i]->GetYaxis()->SetTitle("entries [#]");
    hTime[i]->GetXaxis()->SetTitleSize(0.05);
    hTime[i]->GetYaxis()->SetTitleSize(0.05);
    hTime[i]->GetXaxis()->SetTitleOffset(0.85);
    hTime[i]->GetYaxis()->SetTitleOffset(0.45);
    if(i==0)  hTime[i]->Draw();
    else hTime[i]->Draw("same");   
  }
 
}
