#define prt__beam
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/prttools.C"
#include <TVirtualFitter.h>
#include <TKey.h>
#include <TRandom.h>

void recoPdf(TString in="hits.root", TString pdf="hits.pdf.root", double timeres=0.1, int pid =3, TString nameid="", double r1=0, double r2=0){

  if(!prt_init(in,1)) return;
  TGaxis::SetMaxDigits(4);
  
  TCanvas *cc = new TCanvas("cc","cc");

  TH1F *hl[5],*hll[5],*hnph[5];
  for(int i=0; i<5; i++){
    hl[i] = new TH1F(Form("hl_%d",i),";LE time [ns]; entries [#]", 2000,0,100);    
    hll[i]= new TH1F(Form("ll_i%d",i),";ln L("+prt_lname[pid]+") - ln L(#pi); entries [#]",140,-140,140);
    hnph[i] = new TH1F(Form("hnph_%d",i),";multiplicity [#]; entries [#]", 200,0,200);
    hnph[i]->SetLineColor(prt_color[i]);
    hll[i]->SetLineColor(prt_color[i]);
  }
  
  const int nch(8000);
  TF1 *pdff[nch],*pdfs[nch];
  TH1F *hpdff[nch],*hpdfs[nch];
  TFile f(pdf);

  int rebin=timeres/(100/2000.);
  std::cout<<"rebin "<<rebin <<std::endl;
  
  int integ1(0), integ2(0);
  for(int i=0; i<nch; i++){
    hpdff[i] = (TH1F*)f.Get(Form("hf_%d",i));
    hpdfs[i] = (TH1F*)f.Get(Form("hs_%d",i));
    if(rebin >0) hpdff[i]->Rebin(rebin);
    if(rebin >0) hpdfs[i]->Rebin(rebin);
    integ1+= hpdff[i]->Integral();
    integ2+= hpdfs[i]->Integral();
    hpdff[i]->Smooth(1);
    hpdfs[i]->Smooth(1);
  }

  double time;
  double nph[5]={0};
  int tnph(0),totalf(0),totals(0), ch(0), pdg(0);
  
  for (int ievent=0; ievent<10000 && ievent<prt_entries; ievent++){
    prt_nextEvent(ievent,100);
    //timeres = prt_event->GetTimeRes();
    int id = prt_get_pid(prt_event->GetParticle());    
    double aminf,amins, sum(0),sumf(0),sums(0);
    tnph = 0;    
    if(hll[id]->GetEntries()>4000)continue;

    for(PrtHit hit : prt_event->GetHits()){

      ch=300*hit.GetMcpId()+hit.GetPixelId();      
      time = hit.GetLeadTime() + gRandom->Gaus(0,timeres)+0.05;
      //if(prt_event->GetParticle()==pid) time +=0.01; // remove TOF difference;
      
      //if(time<5 || time>100) continue;      
      aminf = hpdff[ch]->GetBinContent(hpdff[ch]->FindBin(time)); 
      amins = hpdfs[ch]->GetBinContent(hpdfs[ch]->FindBin(time));
      tnph++;
      
      double noise = 1e-4; //1e-7;
      // double noise = 1e-5; //1e-7;
      sumf+=TMath::Log((aminf+noise));
      sums+=TMath::Log((amins+noise));    
            
      hl[id]->Fill(time);

    }
    if(tnph>4) hnph[id]->Fill(tnph);
    sum = sumf-sums;
    if(fabs(sum)<0.1) continue;
    
    //    std::cout<<"sum ===========  "<<sum  << "  "<< sumf<< "  "<< sums<<std::endl;
    
    hll[id]->Fill(sum);     
  }

  gStyle->SetOptStat(0);

  TString name = Form("%1.2f_%1.2f_pik_6gev_hiqe",prt_theta,timeres);
  //TString name = Form("%1.2f_%1.2f",r1,r2);

  prt_canvasAdd("nph_"+name,800,400);
  for(int i=0; i<5; i++) nph[i] = prt_fit(hnph[i],50,20,50).X();
  hnph[2]->Draw();
  hnph[pid]->Draw("same");

  prt_canvasAdd("sep_"+name,800,500);  
  prt_normalize(hll,5);
  
  TF1 *ff;
  double m1(0),m2(0),s1(100),s2(100); 
  if(hll[2]->GetEntries()>10){
    hll[2]->Fit("gaus","Sq");
    ff = hll[2]->GetFunction("gaus");
    ff->SetLineColor(1);
    m1=ff->GetParameter(1);
    s1=ff->GetParameter(2);
  }
  if(hll[pid]->GetEntries()>10){
    hll[pid]->Fit("gaus","Sq");
    ff = hll[pid]->GetFunction("gaus");
    ff->SetLineColor(1);
    m2=ff->GetParameter(1);
    s2=ff->GetParameter(2);
  }
  double sep = (fabs(m1-m2))/(0.5*(s1+s2));
  std::cout<<in<<" separation "<< sep <<std::endl;
  hll[2]->SetTitle(Form("#theta = %1.2f       #sigma = %1.2f",prt_theta, sep));

  
  hll[2]->Draw();
  hll[pid]->Draw("same");

  for(int i=0; i<5; i++){
    hl[i]->Scale(1/hl[i]->GetMaximum());
  }

  prt_normalize(hl,3);
  prt_canvasAdd("tim_"+name,800,400);
  hl[2]->Draw();
  hl[pid]->SetLineColor(4);
  hl[2]->Draw("same");
  hl[pid]->SetLineColor(2);
  hl[pid]->Draw("same");
  prt_canvasSave("data/recoPdf_pik_6.0"+nameid,2);


  TFile fc(in+"_r.root","recreate");
  TTree *tc = new TTree("reco","reco");
  tc->Branch("theta",&prt_theta,"theta/D");
  tc->Branch("sep",&sep,"sep/D");
  tc->Branch("timeres",&timeres,"timeres/D");
  tc->Branch("nph",&nph,"nph[5]/D");
  tc->Branch("r1",&r1,"r1/D");
  tc->Branch("r2",&r2,"r2/D");
  tc->Fill();
  tc->Write();
}
