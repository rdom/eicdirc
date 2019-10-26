#define prt__beam
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/prttools.C"
#include <TVirtualFitter.h>
#include <TKey.h>

void reco_start_time(TString in="hits.root", TString pdf="hits.pdf.root", double timeres=0.1, int pid =321, TString nameid="", double var1=0, double var2=0){

  if(!prt_init(in,1,"data/reco_start_time_"+nameid)) return;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  
  TCanvas *cc = new TCanvas("cc","cc");
  TH1F *hllf= new TH1F("hllf","hllf;ln L(K) - ln L(#pi); entries [#]",180,-100,100);
  TH1F *hlls= new TH1F("hlls","hlls;ln L(K) - ln L(#pi); entries [#]",180,-100,100);  
  TH1F *hl1 = new TH1F("hl1","pdf;LE time [ns]; entries [#]", 2000,0,100);
  TH1F *hl2 = new TH1F("hl2","pdf;LE time [ns]; entries [#]", 2000,0,100);
  TH1F *hl3 = new TH1F("hl3","pdf;LE time [ns]; entries [#]", 2000,0,100);
  TGraph *glhf = new TGraph();
  TGraph *glhs = new TGraph();
  TH1F *hlhf= new TH1F("hlhf","hllf;ln L(K) - ln L(#pi); entries [#]",200,-5,5);
  TH1F *hlhs= new TH1F("hlhs","hlls;ln L(K) - ln L(#pi); entries [#]",200,-5,5);
  TH1F *hnph = new TH1F("hnph",";multiplicity [#]; entries [#]", 200,0,200);
  TH1F *hmeanf = new TH1F("hmeanf","hmean",400,-2,2);
  TH1F *hmeans = new TH1F("hmeans","hmean",400,-2,2);
  hlhf->GetXaxis()->SetRangeUser(-4,4);
  
  const int nch(8000);
  TF1 *pdff[nch],*pdfs[nch];
  TH1F *hpdff[nch],*hpdfs[nch];
  TFile f(pdf);

  int rebin=timeres/(100/2000.);
  std::cout<<"rebin "<<rebin <<std::endl;
  
  if(rebin >0) hl3->Rebin(rebin);
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
    hl3->Add(hpdff[i]);
    hl3->Add(hpdfs[i]);
  }

  TVirtualFitter *fitter;
  double time;
  int totalf(0),totals(0), ch(0), pdg(0);
  double noise = 1e-4;
  
  for (int ievent=0; ievent<1000; ievent++){
    prt_nextEvent(ievent,100);
    pdg = prt_event->GetParticle();
    int nHits =prt_event->GetHitSize();
    double aminf,amins,sum(0),sumf(0),sums(0);
    if(prt_event->GetParticle()==pid && hllf->GetEntries()>1800)continue;
    if(prt_event->GetParticle()==211 && hlls->GetEntries()>1800) continue;  
    hnph->Fill(prt_event->GetHitSize());
    
    int pp = 0;
    for(double s=-5; s<5; s=s+0.05){      
      sumf=0; sums=0;
      for(PrtHit hit : prt_event->GetHits()){
	ch=300*hit.GetMcpId()+hit.GetPixelId();      
	time = hit.GetLeadTime() + prt_rand.Gaus(0,timeres)+s;
      
	aminf = hpdff[ch]->GetBinContent(hpdff[ch]->FindBin(time)); 
	amins = hpdfs[ch]->GetBinContent(hpdfs[ch]->FindBin(time));
      
	sumf+=TMath::Log((aminf+noise));
	sums+=TMath::Log((amins+noise));    

	if(prt_event->GetParticle()==pid) hl1->Fill(time);
	if(prt_event->GetParticle()==211) hl2->Fill(time);
      }
      //std::cout<<s<<" sum "<<sumf<<" "<<sums<<std::endl;
      glhf->SetPoint(pp,s,sumf);
      glhs->SetPoint(pp,s,sums);
      hlhf->Fill(s,sumf);
      hlhs->Fill(s,sums);
      pp++;
    }

    TF1 *flh =new TF1("flh","gaus(0)+pol1(3)",-4,4);
    flh->SetParameters(1,0,1,1,1);
    flh->SetParameters(1,0,1,1,1);
    flh->SetParLimits(1,-1,1);
    flh->SetParLimits(2,0.5,5);
    hlhf->Fit("flh","","",-4,4);

    double mean = flh->GetParameter(1);
    hmeanf->Fill(mean);
    hlhs->Fit("flh","","",-4,4);
    mean = flh->GetParameter(1);
    hmeans->Fill(mean);

    // gPad->Update();
    // gPad->WaitPrimitive();
    hlhf->Reset();
    hlhs->Reset();
    
    
    sum = sumf-sums;
    if(fabs(sum)<0.1) continue;
    
    if(prt_event->GetParticle()==pid) hllf->Fill(sum);
    if(prt_event->GetParticle()==211) hlls->Fill(sum);     
  }

  
  TString name = Form("%1.2f_%1.2f",prt_theta,timeres);  
  // prt_canvasAdd("scan_"+name,600,500);
  // glhf->SetLineColor(kRed+1);
  // glhf->Draw("apl");
  // glhs->SetLineColor(kBlue+1);
  // glhs->Draw("pl same");


  // hlhf->GetXaxis()->SetRangeUser(-4,4);
  // hlhf->SetLineColor(kRed+1);
  // TF1 *flh =new TF1("flh","gaus(0)+pol1(3)",-4,4);
  // flh->SetParameters(1,0,1,1,1);
  // flh->SetParameters(1,0,1,1,1);
  // flh->SetParLimits(1,-1,1);
  // flh->SetParLimits(2,0.5,5);
  // hlhf->Fit("flh","","",-4,4);
  // hlhf->Draw("hist");
  // flh->Draw("same");
  // hlhs->SetLineColor(kBlue+1);
  // hlhs->Draw("hist same");

  prt_canvasAdd("nph_"+name,800,400);
  hnph->Draw();  
  double nph = prt_fit(hnph,50,20,50).X();

  prt_canvasAdd("mean_"+name,600,500);
  var1 = prt_fit(hmeanf,5,20,5).Y();
  var2 = prt_fit(hmeans,5,20,5).Y();
  double var3 = prt_fit(hmeanf,5,20,5).X();
  double var4 = prt_fit(hmeans,5,20,5).X();
  std::cout<<var1<<" var "<<var2<<std::endl;
  
  hmeanf->SetLineColor(kRed+1);
  hmeanf->Draw();
  hmeans->SetLineColor(kBlue+1);
  hmeans->Draw("same");  

  prt_canvasAdd("sep_"+name,600,500);
  prt_normalize(hllf,hlls);
  
  TF1 *ff;
  double m1(0),m2(0),s1(100),s2(100); 
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
  double sep = (fabs(m1-m2))/(0.5*(s1+s2));
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
  tc->Branch("nph",&nph,"nph/D");
  tc->Branch("var1",&var1,"var1/D");
  tc->Branch("var2",&var2,"var2/D");
  tc->Branch("var3",&var3,"var3/D");
  tc->Branch("var4",&var4,"var4/D");
  tc->Fill();
  tc->Write();
}
