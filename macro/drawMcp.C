#include "../src/PrtHit.h"
#include "../src/PrtEvent.h" 

#include "../../prttools/prttools.C"

const Int_t nmcp = 150;
const Int_t nrow = 3;
const Int_t ncol = 5;

TH2D* hist[nmcp];

void drawMcp(TString path = ".", TString name = ""){
  prt_savepath = "load";
  
  //gStyle->SetOptStat(0);
  Int_t referenceTime, referenceChannel;
  TH2F *hist1 = new TH2F("hHist1",";y [cm];x [cm]",500,-20,20,500,-2,32 );
  TH2F *hist2 = new TH2F("hHist2",";y [cm];x [cm]",500,-20,20,500,-2,32 );
  TH2F *hVertex = new TH2F("hVertex","",5000,-2,2, 1000,-10,10 );

  hist1->SetStats(0);
  hVertex->SetStats(0);
  prt_setRootPalette(1);
  
  TH1F *hTime  = new TH1F("hTime",";time [ns];entries [#]",500,0,150);
  Int_t pdg(0), angle(0), startInd(0);

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
  ch->Add("hits.root");
  ch->SetBranchAddress("PrtEvent", &fEvent);
  std::cout<<"ch->GetEntries()  "<< ch->GetEntries()<<std::endl;
 
  for (Int_t ievent=0; ievent< ch->GetEntries()
	 ; ievent++){
    ch->GetEntry(ievent);
    if(ievent%100==0) cout<<"Event # "<<ievent<<endl;
    for(Int_t h=0; h<fEvent->GetHitSize(); h++){
      hit = fEvent->GetHit(h);
      pdg = fEvent->GetParticle();
      angle =fEvent->GetAngle() + 0.01;
      
      hist[hit.GetMcpId()-startInd]->Fill((hit.GetPixelId()-startInd)/8,
					  (hit.GetPixelId()-startInd)%8);
      if(pdg==211) hist1->Fill(hit.GetGlobalPos().Y()/10.,hit.GetGlobalPos().X()/10.);
      if(pdg==321) hist2->Fill(hit.GetGlobalPos().Y()/10.,hit.GetGlobalPos().X()/10.);
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

 prt_canvasAdd(Form("loadi_%d",angle),800,500);
 hist1->SetMarkerColor(4);
 hist2->SetMarkerColor(2);
 hist1->Draw();
 hist2->Draw("same");
 TGaxis::SetMaxDigits(3);
 gStyle->SetOptStat(111);

 prt_canvasAdd(Form("time_%d",angle),800,500);
 hTime->SetTitle(0);
 hTime->Draw();

 prt_canvasSave(1,0); 
}


