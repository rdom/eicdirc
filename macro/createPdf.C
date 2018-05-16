#define prt__beam
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/prttools.C"

TSpectrum *spect = new TSpectrum(2);
TF1 * fitpdf(TH1F *h){
  TF1 *gaust = NULL;
  Double_t rmin(0),rmax(60);
  Double_t integral = h->Integral(h->GetXaxis()->FindBin(rmin),h->GetXaxis()->FindBin(rmax));
  if(h->GetEntries()>50){
    Int_t nfound = spect->Search(h,2,"",0.1);
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

void createPdf(TString in="hits.root",int end=0){
  if(!prt_init(in,1,"data/createPdf")) return;
   
  const Int_t nch(15000);
  TH1F *hlef[nch], *hles[nch];

  for(Int_t i=0; i<nch; i++){
    hlef[i] = new TH1F(Form("lef_%d",i),"pdf;LE time [ns]; entries [#]", 2000,0,100);
    hles[i] = new TH1F(Form("les_%d",i),"pdf;LE time [ns]; entries [#]", 2000,0,100);
  }
  
  Double_t time;
  PrtHit hit;
  if(end==0) end=prt_entries;
  Int_t pdg(0), totalf(0),totals(0), ch;
  for (Int_t e=0; e<prt_entries; e++){ //prt_entries
    prt_nextEvent(e,1000);
    pdg =prt_event->GetParticle();
    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      hit = prt_event->GetHit(i);
      ch=300*hit.GetMcpId()+hit.GetPixelId();      
      time = hit.GetLeadTime();//+gRandom->Gaus(0,0.1);
      if(pdg==321){
	//totalf++;
	hlef[ch]->Fill(time);
      }
      if(pdg==211){
	//totals++;
	hles[ch]->Fill(time);
      }
    }
    if(pdg==321) totalf++;
    if(pdg==211 ) totals++;
    if(totalf>end || totals>end) break;
  }

  std::cout<<"#1 "<< totalf <<"  #2 "<<totals <<std::endl;
  //TCanvas *cExport = new TCanvas("cExport","cExport",0,0,800,400);
  
  if(totalf>0 && totals>0) {
    in.ReplaceAll(".root",".pdf.root");
    TFile efile(in,"RECREATE");
    
    for(Int_t i=0; i<nch; i++){
      hlef[i]->Scale(1/(Double_t)totalf);
      hles[i]->Scale(1/(Double_t)totals);
      
      hlef[i]->SetName(Form("hf_%d",i));
      hles[i]->SetName(Form("hs_%d",i));
      hlef[i]->Write();
      hles[i]->Write();
      
      // if(false){
      // 	cExport->cd();
      // 	//	canvasAdd(Form("pdf_%d",i),800,500);
      // 	cExport->SetName(Form("pdf_%d",i));
      // 	//canvasAdd(cExport);
      // 	//hlef[i]->GetYaxis()->SetRangeUser(0,1.5);
      // 	prt_normalize(hlef[i],hles[i]);
      // 	prt_axisTime800x500(hlef[i]);
      // 	prt_axisTime800x500(hles[i]);
      // 	hlef[i]->SetLineColor(2);
      // 	hlef[i]->Draw();
      // 	hles[i]->SetLineColor(4);
      // 	hles[i]->Draw("same");
      // 	// // f->Draw();
      // 	// s->SetLineColor(4);
      // 	// s->Draw("same");
      // 	cExport->Print(prt_savepath+Form("/pdf_mcp%d_pix_%d.png",map_mcp[i],map_pix[i]));
      // 	//canvasSave(1,0);
      // 	//canvasDel(cExport->GetName());
      // }
    }
    
    efile.Write();
    efile.Close();
  }

  // prt_drawDigi("m,p,v\n",prt_geometry,0,0);
  // prt_cdigi->SetName("hits");
  // prt_canvasAdd(prt_cdigi);
  // prt_canvasSave(1,0);
}
