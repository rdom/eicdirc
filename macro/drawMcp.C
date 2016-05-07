// #include "TROOT.h"
// #include "TSystem.h"
// #include "TStyle.h"
// #include "TCanvas.h"
// #include "TPad.h"
// #include "TH1.h"
// #include "TH2.h"
// #include "TF1.h"
// #include "TFile.h"
// #include "TTree.h"
// #include "TClonesArray.h"
// #include "TVector3.h"
// #include "TMath.h"
// #include "TChain.h"

// #include <vector>
// #include <iostream>

// #include "prt/PrtHit.h"
// #include "prt/PrtEvent.h" 

#include "../../prttools/prttools.C"

const Int_t nmcp = 15;
const Int_t nrow = 3;
const Int_t ncol = 5;

TH2D* hist[nmcp];

void drawMcp(TString path = ".", TString name = ""){
  fSavePath = "load";
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  //gStyle->SetOptStat(0);
  Int_t referenceTime, referenceChannel;
  TH2F *hist1 = new TH2F("hHist",";y [cm];x [cm]",500,-20,20,500,-2,32 );
  TH2F *hVertex = new TH2F("hVertex","",5000,-2,2, 1000,-10,10 );

  hist1->SetStats(0);
  hVertex->SetStats(0);

  TH1F *hTime  = new TH1F("hTime",";time [ns];entries [#]",500,0,150);
  Int_t angle, startInd = 0;

  for(int i=0; i<nmcp;i++){
    hist[i] = new TH2D( Form("mcp - %d", i+startInd),Form("%d", i+startInd),8,0.,8.,8,0.,8.);
    hist[i]->SetStats(0);
    hist[i]->GetXaxis()->SetNdivisions(10);
    hist[i]->GetYaxis()->SetNdivisions(10);
    hist[i]->GetXaxis()->SetLabelOffset(100);
    hist[i]->GetYaxis()->SetLabelOffset(100);
    hist[i]->GetXaxis()->SetTickLength(1);
    hist[i]->GetYaxis()->SetTickLength(1);
    hist[i]->GetXaxis()->SetAxisColor(15);
    hist[i]->GetYaxis()->SetAxisColor(15);
    // hist[i]->SetLineColor(ci);
  }

  PrtEvent* fEvent = 0;
  PrtHit hit;
  TChain *ch = new TChain("data");
  ch->Add("../build/hits.root");
  ch->SetBranchAddress("PrtEvent", &fEvent);
  std::cout<<"ch->GetEntries()  "<< ch->GetEntries()<<std::endl;
 
  for (Int_t ievent=0; ievent<ch->GetEntries(); ievent++){
    ch->GetEntry(ievent);
    if(ievent%100==0) cout<<"Event # "<<ievent<<endl;
    for(Int_t h=0; h<fEvent->GetHitSize(); h++){
      hit = fEvent->GetHit(h);
      angle =fEvent->GetAngle() + 0.01;
      
      hist[hit.GetMcpId()-startInd]->Fill((hit.GetPixelId()-startInd)/8,
					  (hit.GetPixelId()-startInd)%8);
      hist1->Fill(hit.GetGlobalPos().Y()/10.,hit.GetGlobalPos().X()/10.);
      hTime->Fill(hit.GetLeadTime());
    }
  }  


  int tmax, max=0;
  for(int i=0; i<9;i++){
    tmax = hist[i]->GetMaximum();
    if(max<tmax) max = tmax;
  }


 //  canvasAdd(Form("load_%d",angle),800,500);
 //  TPad* pads[nmcp];
 //  float tbw, tbh, bw = 0.01, bh = 0.05, ntw = 3, nth = 3;
 //  float tbw = bw, tbh = bh;
 //  int padi = 0;
 //  for(int i=0; i<nrow; i++){
 //    for(int j=0; j<ncol; j++){
 //      //  if(i==0 || i==2) {tbw = 0.001; }
 //      // else{tbw = bw;}
 //      pads[padi] =  new TPad("P","T", i/ncol+tbw , j/nrow+tbw, (i+1)/ncol-tbw, (1+j)/nrow-tbh, 21);
 //      pads[padi]->SetFillColor(kCyan-5);
 //      //pads[padi]->SetFrameFillColor(ci);
 //      pads[padi]->SetMargin(0.04,0.04,0.04,0.04);
 //      pads[padi]->Draw(); 
 //      padi++;
 //    }
 //  }
 // for(int i=0; i<nmcp;i++){
 //    pads[i]->cd();
 //    hist[i]->Draw("col");
 //    hist[i]->SetMaximum(max);
 // }

 canvasAdd(Form("loadi_%d",angle),800,500);
 hist1->Draw("colz");
 TGaxis::SetMaxDigits(3);
 gStyle->SetOptStat(111);

 canvasAdd(Form("time_%d",angle),800,500);
 hTime->SetTitle(0);
 hTime->Draw();

 canvasSave(1,0); 
}


