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

const Int_t nmcp = 15;
const Int_t nrow = 3;
const Int_t ncol = 5;

TH2D* hist[nmcp];
#include "save.C"

void drawMcp(TString path = ".", TString name = ""){

  Int_t saveflag = 1;
  TString info = "hits";
  TString path = createDir("rdata", info, saveflag); 

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  gSystem->Load("../build/libeicdirclib.so");
  //gStyle->SetOptStat(0);
  Int_t referenceTime;
  Int_t referenceChannel;
  //TH2F *hist1 = new TH2F("name","name",500,-32,19,500,-10,10 );
  // TH2F *hist1 = new TH2F("hHist","",500,-2,32,500,-10,10 );
 TH2F *hist1 = new TH2F("hHist","",500,-2,32,500,-20,20 );
  TH2F *hVertex = new TH2F("hVertex","",5000,-2,2, 1000,-10,10 );

  hist1->SetStats(0);
  hVertex->SetStats(0);

  TH1F *hTime[4];
  for(Int_t i=0; i<4; i++){
    hTime[i] = new TH1F(Form("id - %d", i),Form("id%d", i),500,0,25);
    hTime[i]->SetLineColor(i+1);
  }
  Int_t startInd = 0;

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
      //	std::cout<<"x "<<hit.GetDigiPos().X() <<std::endl;
      // hist[hit.GetMcpId()]->Fill(hit.GetLocalPos().X(),hit.GetLocalPos().Y());
      // if(hit.GetPixelId()/8==7 && hit.GetPixelId()%8 ==0)
      //if(hit.GetMcpId() ==6)
      Double_t time = hit.GetLeadTime(0);
      { 
	if(time<5000){
	  hist[hit.GetMcpId()-startInd]->Fill((hit.GetPixelId()-startInd)/8, (hit.GetPixelId()-startInd)%8);
	  hist1->Fill(hit.GetGlobalPos().X()/10.,hit.GetGlobalPos().Y()/10.);
	}
	//if(hit.GetPixelId()/8==7 && hit.GetPixelId()%8 ==0)  
	// if(hit.GetLeadTime(0)>6000)
	// if(hit.GetNreflectionsInPrizm() < 5)
	hVertex->Fill(hit.GetPosition().Y(),hit.GetPosition().Z());
	  
	for(Int_t i=0; i<4; i++){
	  hTime[i]->Fill(hit.GetLeadTime(i)/1000.); 
	}
      }
    }
  }  


  int tmax, max=0;
  for(int i=0; i<9;i++){
    tmax = hist[i]->GetMaximum();
    if(max<tmax) max = tmax;
  }

  TCanvas* c1 = new TCanvas("c1","c1",500,600);

  TPad* pads[nmcp];
  float tbw, tbh, bw = 0.01, bh = 0.05, ntw = 3, nth = 3;
  float tbw = bw, tbh = bh;
  int padi = 0;
  for(int i=0; i<nrow; i++){
    for(int j=0; j<ncol; j++){
      //  if(i==0 || i==2) {tbw = 0.001; }
      // else{tbw = bw;}
      pads[padi] =  new TPad("P","T", i/ncol+tbw , j/nrow+tbw, (i+1)/ncol-tbw, (1+j)/nrow-tbh, 21);
      pads[padi]->SetFillColor(kCyan-5);
      //pads[padi]->SetFrameFillColor(ci);
      pads[padi]->SetMargin(0.04,0.04,0.04,0.04);
      pads[padi]->Draw(); 
      padi++;
    }
  }
 for(int i=0; i<nmcp;i++){
    pads[i]->cd();
    hist[i]->Draw("col");
    hist[i]->SetMaximum(max);
 }
 Int_t nofe = hist1->GetEntries();
 std::cout<<"nofe  "<<nofe <<std::endl;

 TCanvas* c2 = new TCanvas("c2","c2",600,900);
 hist1->SetTitle(Form("%d hits",nofe));
 hist1->GetXaxis()->SetTitle("x, [cm]");
 hist1->GetYaxis()->SetTitle("y, [cm]");
 hist1->GetXaxis()->SetTitleOffset(1.15);
 hist1->GetYaxis()->SetTitleOffset(1.2);
 hist1->Draw("colz");
 TGaxis::SetMaxDigits(3);
 gStyle->SetOptStat(111);
 TCanvas* c3 = new TCanvas("c3","c3",900,400);
 for(Int_t i=0; i<4; i++){
   hTime[i]->SetTitle(0);
   hTime[i]->GetXaxis()->SetTitle("time, [ns]");
   hTime[i]->GetYaxis()->SetTitle("entries, [#]");
   hTime[i]->GetXaxis()->SetTitleSize(0.05);
   hTime[i]->GetYaxis()->SetTitleSize(0.05);
   hTime[i]->GetXaxis()->SetTitleOffset(0.85);
   hTime[i]->GetYaxis()->SetTitleOffset(0.45);
   if(i==0)  hTime[i]->Draw();
   //else hTime[i]->Draw("same");   
 }

 TCanvas* c4 = new TCanvas("c4","c4",900,600);
 TPad *pad1 = new TPad("pad1", "The pad 60% of the height", 0.7,0, 1.0,1.0, 0);
 TPad *pad2 = new TPad("pad2", "The pad 40% of the height", 0,0.0, 0.7,1.0, 0);
 pad1->Draw();
 pad2->Draw();

 pad2->cd();
 hist1->Draw("colz");

 pad1->cd();
 hVertex->GetXaxis()->SetTitle("y, [mm]");
 hVertex->GetYaxis()->SetTitle("z, [mm]");
 hVertex->Draw("colz");

 c1->Print(path+"/digi"+name+".png");
 c2->Print(path+"/sim"+name+".png");
 c3->Print(path+"/time"+name+".png");
}


