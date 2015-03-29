#include "../../prttools/prttools.C"
void da_scan(TString inFile = "r_spr.root", TString outFile="c_spr.root"){

  TChain ch("dirc"); ch.Add(inFile);
  Double_t cangle,spr,trr,nph,par1,par2,par3,par4,par5,par6,par7,par8,theta,phi; 

  TGraph *gSpr = new TGraph();
  TGraph *gNph = new TGraph();
  TGraph *gTrr = new TGraph();

  ch.SetBranchAddress("spr",&spr);
  ch.SetBranchAddress("trr",&trr);
  ch.SetBranchAddress("nph",&nph);
  ch.SetBranchAddress("cangle",&cangle);
  ch.SetBranchAddress("par4",&par4);
  ch.SetBranchAddress("par5",&par5);
  ch.SetBranchAddress("par6",&par6);
  ch.SetBranchAddress("par7",&par7);
  ch.SetBranchAddress("par8",&par8);
  ch.SetBranchAddress("theta",&theta);
  ch.SetBranchAddress("phi",&phi);
  
  Int_t nent = ch.GetEntries();
  std::cout<<"# entries  "<< nent <<std::endl;
  for (Int_t i = 0; i < nent; i++) {
    ch.GetEvent(i);
    gSpr->SetPoint(i,theta,TMath::Abs(spr));
    gNph->SetPoint(i,theta,nph);
    gTrr->SetPoint(i,theta,TMath::Abs(trr));
  }
  gSpr->Sort();
  gNph->Sort();
  gTrr->Sort();

  gSpr->SetLineColor(38);
  gNph->SetLineColor(38);
  gTrr->SetLineColor(38);
  gSpr->SetMarkerStyle(20);
  gNph->SetMarkerStyle(20);
  gTrr->SetMarkerStyle(20);
  gSpr->SetMarkerSize(0.7);
  gNph->SetMarkerSize(0.7);
  gTrr->SetMarkerSize(0.7);
  gNph->GetYaxis()->SetRangeUser(0,150);
  gSpr->GetYaxis()->SetRangeUser(0,12);
  gTrr->GetYaxis()->SetRangeUser(0,2);

  gSpr->GetYaxis()->SetTitle("SPR [mrad]");
  gNph->GetYaxis()->SetTitle("multiplicity [#]");
  gTrr->GetYaxis()->SetTitle("#sigma_{#theta_{C} tr} [mrad]");
  
  gSpr->GetXaxis()->SetLabelSize(0.05);
  gSpr->GetXaxis()->SetTitleSize(0.06);
  gSpr->GetXaxis()->SetTitleOffset(0.84);

  gTrr->GetXaxis()->SetLabelSize(0.05);
  gTrr->GetXaxis()->SetTitleSize(0.06);
  gTrr->GetXaxis()->SetTitleOffset(0.84);

  gNph->GetXaxis()->SetLabelSize(0.05);
  gNph->GetXaxis()->SetTitleSize(0.06);
  gNph->GetXaxis()->SetTitleOffset(0.84);

  gSpr->GetYaxis()->SetLabelSize(0.05);
  gSpr->GetYaxis()->SetTitleSize(0.06);
  gSpr->GetYaxis()->SetTitleOffset(0.7);

  gTrr->GetYaxis()->SetLabelSize(0.05);
  gTrr->GetYaxis()->SetTitleSize(0.06);
  gTrr->GetYaxis()->SetTitleOffset(0.7);

  gNph->GetYaxis()->SetLabelSize(0.05);
  gNph->GetYaxis()->SetTitleSize(0.06);
  gNph->GetYaxis()->SetTitleOffset(0.7);


  gSpr->GetXaxis()->SetTitle("#theta_{track} [#circ]");
  gNph->GetXaxis()->SetTitle("#theta_{track} [#circ]");
  gTrr->GetXaxis()->SetTitle("#theta_{track} [#circ]");

  TFile *file = new TFile(outFile,"RECREATE");
  TCanvas* c1 = new TCanvas("c1","c1",800,500);c1->SetBottomMargin(0.12);
  gSpr->Draw("APL");
  canvasAdd(c1);
  TCanvas* c2 = new TCanvas("c2","c2",800,500);c2->SetBottomMargin(0.12);
  gNph->Draw("APL");
  canvasAdd(c2);
  TCanvas* c3 = new TCanvas("c3","c3",800,500);c3->SetBottomMargin(0.12);
  gTrr->Draw("APL");
  canvasAdd(c3);

  canvasSave(0,"addcanvases.C",1,"data/temp");
  
  file->cd();
  c1->Write();
  c2->Write();
  c3->Write();
  file->Close();
}
