#define prt__beam
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/prttools.C"

TSpectrum *spect = new TSpectrum(2);
TF1 * fitpdf(TH1F *h){
  TF1 *gaust = NULL;
  double rmin(0),rmax(60);
  double integral = h->Integral(h->GetXaxis()->FindBin(rmin),h->GetXaxis()->FindBin(rmax));
  if(h->GetEntries()>50){
    int nfound = spect->Search(h,2,"",0.1);
    std::cout<<"nfound  "<<nfound <<std::endl;
    if(nfound==1){
      gaust =new TF1("gaust","gaus(0)",rmin,rmax);
      gaust->SetNpx(500);
      gaust->FixParameter(1,spect->GetPositionX()[0]);
      gaust->SetParameter(2,0.3);
      gaust->SetParLimits(2,0.2,1);
    }
    if(nfound==2){
      gaust =new TF1("gaust","gaus(0)+gaus(3)",rmin,rmax);
      gaust->SetNpx(500);
      gaust->FixParameter(1,spect->GetPositionX()[0]);
      gaust->FixParameter(4,spect->GetPositionX()[1]);
      gaust->SetParameter(2,0.3);
      gaust->SetParameter(5,0.3);
      gaust->SetParLimits(2,0.2,1);
      gaust->SetParLimits(5,0.2,1);
    
      std::cout<<spect->GetPositionX()[0]<< " "<<spect->GetPositionX()[1] <<std::endl;
    
    }
    h->Fit("gaust","","MQN",0,60);
  }else{
    gaust =new TF1("gaust","pol0",rmin,rmax);
    gaust->FixParameter(0,0);
  }
	 
  return gaust;
}

void createPdf(TString in="hits.root",int end=0, int pid=321){
  if(!prt_init(in,1,"data/createPdf")) return;
   
  const int nch(24*256);
  TH1F *hlef[nch], *hles[nch];

  for(int i=0; i<nch; i++){
    hlef[i] = new TH1F(Form("lef_%d",i),";LE time [ns]; entries/N [#]", 2000,0,100);
    hles[i] = new TH1F(Form("les_%d",i),";LE time [ns]; entries/N [#]", 2000,0,100);
  }
  
  double time;
  PrtHit hit;
  if(end==0) end=prt_entries;
  int pdg(0), totalf(0),totals(0), ch;
  for (int e=10000; e<prt_entries; e++){ //prt_entries
    prt_nextEvent(e, 1000);
    pdg = prt_event->GetParticle();

    if(prt_event->GetHitSize()<5) continue;
    for(int i=0; i<prt_event->GetHitSize(); i++){
      hit = prt_event->GetHit(i);
      ch = hit.GetMcpId()*256+hit.GetPixelId();      
      time = hit.GetLeadTime()+gRandom->Gaus(0,0.1);
      
      //if(prt_event->GetParticle()==pid) time +=0.01; // remove TOF difference;
      if(pdg==pid){
	totalf++;
	hlef[ch]->Fill(time);
      }
      if(pdg==211){
	totals++;
	hles[ch]->Fill(time);
      }
    }
    // if(pdg==pid) totalf++;
    // if(pdg==211 ) totals++;
    // if(totalf>end || totals>end) break;
  }

  std::cout<<"#1 "<< totalf <<"  #2 "<<totals <<std::endl;
  
  if(totalf>=0 && totals>0) {
    in.ReplaceAll(".root",".pdf.root");
    TFile efile(in,"RECREATE");
    
    for(int i=0; i<nch; i++){
      hlef[i]->Scale(1/(double)totalf);
      hles[i]->Scale(1/(double)totals);
      
      hlef[i]->SetName(Form("hf_%d",i));
      hles[i]->SetName(Form("hs_%d",i));
      hlef[i]->Write();
      hles[i]->Write();
      
      if(0){
	int nrow=4, ncol=6, p = i/256;
	int np =p%ncol*nrow + p/ncol;
	hles[i]->SetName(Form("mcp%dpix%d",np,i%256));
	prt_canvasAdd(Form("pdf_mcp%dpix%d",np,i%256),800,400);
      	prt_normalize(hlef[i],hles[i]);
      	hlef[i]->SetLineColor(2);
	hles[i]->SetLineColor(4);
	hles[i]->GetXaxis()->SetRangeUser(5, 20);
      	hles[i]->Draw("hist");
      	hlef[i]->Draw("hist same");
	prt_canvasSave("data/pdfs_7",0,0,1);
      }
    }
    
    efile.Write();
    efile.Close();
  }

  // prt_drawDigi("m,p,v\n",prt_geometry,0,0);
  // prt_cdigi->SetName("hits");
  // prt_canvasAdd(prt_cdigi);
  // prt_canvasSave(1,0);
}
