#define prt__beam
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/prttools.C"
#include <TVirtualFitter.h>
#include <TKey.h>
#include <TRandom.h>

void recoPdf(TString in="hits.root", TString pdf="hits.pdf.root", Double_t timeres=0.1, int pid =321, TString nameid="", Double_t r1=0, Double_t r2=0){

  if(!prt_init(in,1,"data/recoPdf_"+nameid)) return;
  TGaxis::SetMaxDigits(4);
  
  TCanvas *cc = new TCanvas("cc","cc");
  TH1F *hllf= new TH1F("hllf","hllf;ln L(K) - ln L(#pi); entries [#]",180,-100,100);
  TH1F *hlls= new TH1F("hlls","hlls;ln L(K) - ln L(#pi); entries [#]",180,-100,100);  
  TH1F *hl1 = new TH1F("hl1","pdf;LE time [ns]; entries [#]", 2000,0,100);
  TH1F *hl2 = new TH1F("hl2","pdf;LE time [ns]; entries [#]", 2000,0,100);
  TH1F *hl3 = new TH1F("hl3","pdf;LE time [ns]; entries [#]", 2000,0,100);

  TH1F *hnph = new TH1F("hnph",";multiplicity [#]; entries [#]", 200,0,200);
  
  TRandom rand;
  const Int_t nch(8000);
  TF1 *pdff[nch],*pdfs[nch];
  TH1F *hpdff[nch],*hpdfs[nch];
  TFile f(pdf);

  Int_t rebin=timeres/(100/2000.);
  std::cout<<"rebin "<<rebin <<std::endl;
  
  if(rebin >0) hl3->Rebin(rebin);
  Int_t integ1(0), integ2(0);
  for(Int_t i=0; i<nch; i++){
    hpdff[i] = (TH1F*)f.Get(Form("hf_%d",i));
    hpdfs[i] = (TH1F*)f.Get(Form("hs_%d",i));
    if(rebin >0) hpdff[i]->Rebin(rebin);
    if(rebin >0) hpdfs[i]->Rebin(rebin);
    integ1+= hpdff[i]->Integral();
    integ2+= hpdfs[i]->Integral();
    hpdff[i]->Smooth(1);
    hpdfs[i]->Smooth(1);
    hl3->Add(hpdff[i]);
    hl3->Add(hpdfs[i]);
  }

  TVirtualFitter *fitter;
  Double_t time;
  PrtHit hit;
  Int_t tnph(0),totalf(0),totals(0), ch(0), pdg(0);
  
  for (Int_t ievent=0; ievent<4000; ievent++){
    prt_nextEvent(ievent,100);
    //timeres = prt_event->GetTimeRes();
    pdg = prt_event->GetParticle();
    Double_t aminf,amins, sum(0),sumf(0),sums(0);
    Int_t nHits =prt_event->GetHitSize();
    tnph = 0;
    
    if(prt_event->GetParticle()==pid && hllf->GetEntries()>1800)continue;
    if(prt_event->GetParticle()==211 && hlls->GetEntries()>1800) continue;  

    for(Int_t i=0; i<nHits; i++){
      hit = prt_event->GetHit(i);
      ch=300*hit.GetMcpId()+hit.GetPixelId();
      
      time = hit.GetLeadTime() + rand.Gaus(0,timeres);
      //if(prt_event->GetParticle()==pid) time +=0.01; // remove TOF difference;
      
      //if(time<5 || time>100) continue;
      
      aminf = hpdff[ch]->GetBinContent(hpdff[ch]->FindBin(time)); 
      amins = hpdfs[ch]->GetBinContent(hpdfs[ch]->FindBin(time));
      tnph++;
      
      Double_t noise = 1e-4; //1e-7;
      //Double_t noise = 1e-5; //1e-7;
      sumf+=TMath::Log((aminf+noise));
      sums+=TMath::Log((amins+noise));    

      // std::cout<< "K="<<aminf<< " pi="<<amins << " pid="<< pdg <<std::endl;
      // //  if(aminf>amins)
      // //	if(aminf==0 || amins==0)
      // {
      // 	prt_normalize(hpdff[ch],hpdfs[ch]);
      // 	hpdff[ch]->SetLineColor(2);
      // 	hpdff[ch]->Draw("hist");
      // 	hpdfs[ch]->SetLineColor(4);
      // 	hpdfs[ch]->Draw("hist same");
      // 	cc->Update();
      // 	TLine *line = new TLine(time,0,time,cc->GetUymax());
      // 	line->Draw("same");
      // 	cc->Update();
      // 	cc->WaitPrimitive();
      // }

      if(prt_event->GetParticle()==pid) hl1->Fill(time);
      if(prt_event->GetParticle()==211) hl2->Fill(time);

    }
    hnph->Fill(tnph);
    sum = sumf-sums;
    if(fabs(sum)<0.1) continue;
    
    //    std::cout<<"sum ===========  "<<sum  << "  "<< sumf<< "  "<< sums<<std::endl;
    
    if(prt_event->GetParticle()==pid) hllf->Fill(sum);
    if(prt_event->GetParticle()==211) hlls->Fill(sum);
     
  }

  gStyle->SetOptStat(0);

  TString name = Form("%1.2f_%1.2f",prt_theta,timeres);

  prt_canvasAdd("nph_"+name,800,400);
  hnph->Draw();  
  Int_t nph = prt_fit(hnph,50,20,50).X();
  
  prt_canvasAdd("sep_"+name,600,500);

  prt_normalize(hllf,hlls);
  
  TF1 *ff;
  Double_t m1(0),m2(0),s1(100),s2(100); 
  if(hllf->GetEntries()>10){
    hllf->Fit("gaus","Sq");
    ff = hllf->GetFunction("gaus");
    ff->SetLineColor(1);
    m1=ff->GetParameter(1);
    s1=ff->GetParameter(2);
  }
  if(hlls->GetEntries()>10){
    hlls->Fit("gaus","Sq");
    ff = hlls->GetFunction("gaus");
    ff->SetLineColor(1);
    m2=ff->GetParameter(1);
    s2=ff->GetParameter(2);
  }
  Double_t sep = (fabs(m1-m2))/(0.5*(s1+s2));
  std::cout<<in<<" separation "<< sep <<std::endl;
  hllf->SetTitle(Form("#theta = %1.2f       #sigma = %1.2f",prt_theta, sep));

  
  hllf->SetLineColor(2);
  hllf->Draw();
  hlls->SetLineColor(4);
  hlls->Draw("same");

  hl1->Scale(1/hl1->GetMaximum());
  hl2->Scale(1/hl2->GetMaximum());
  hl3->Scale(1/hl3->GetMaximum());

  prt_normalize(hl1,hl2);
  prt_canvasAdd("tim_"+name,800,500);
  hl1->Draw();
  hl2->SetLineColor(4);
  hl2->Draw("same");
  hl3->SetLineColor(2);
  hl3->Draw("same");
  prt_canvasSave(2,0);


  TFile fc(prt_savepath+"/res_"+name+".root","recreate");
  TTree *tc = new TTree("reco","reco");
  tc->Branch("theta",&prt_theta,"prt_theta/D");
  tc->Branch("sep",&sep,"sep/D");
  tc->Branch("timeres",&timeres,"timeres/D");
  tc->Branch("nph",&nph,"nph/I");
  tc->Branch("r1",&r1,"r1/D");
  tc->Branch("r2",&r2,"r2/D");
  tc->Fill();
  tc->Write();
}
