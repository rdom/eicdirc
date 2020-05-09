#include "../../prttools/prttools.C"
void da_scan(TString inFile = "r_spr.root", TString outFile="c_spr.root"){
  prt_savepath ="data/da_scan";
  TChain ch("reco"); ch.Add(inFile);
  Double_t cangle,spr,trr,nph[5],par1,par2,par3,par4,par5,par6,test1,test2,theta,phi,sep; 

  TGraph *gSpr = new TGraph();
  TGraph *gNph = new TGraph();
  TGraph *gTrr = new TGraph();
  TGraph *gTrr1 = new TGraph();
  TGraph *gTrr2 = new TGraph();
  TGraph *gSep = new TGraph();

  ch.SetBranchAddress("spr",&spr);
  ch.SetBranchAddress("trr",&trr);
  ch.SetBranchAddress("sep",&sep);
  ch.SetBranchAddress("nph",&nph);
  ch.SetBranchAddress("cangle",&cangle);
  ch.SetBranchAddress("par5",&par5);
  ch.SetBranchAddress("par6",&par6);
  ch.SetBranchAddress("test1",&test1);
  ch.SetBranchAddress("test2",&test2);
  ch.SetBranchAddress("theta",&theta);
  ch.SetBranchAddress("phi",&phi);
  
  int pid=0, nent = ch.GetEntries();
  std::cout<<"# entries  "<< nent <<std::endl;  
  for (int i = 0; i < nent; i++) {
    ch.GetEvent(i);

    if(theta==80 || theta==100) continue;
    std::cout<<"theta "<<theta<<std::endl;
    
    
    gSpr->SetPoint(pid,theta,TMath::Abs(spr));
    gNph->SetPoint(pid,theta,nph[2]);
    gSep->SetPoint(pid,theta,sep);
    trr = TMath::Abs(trr);
    gTrr->SetPoint(pid,theta,trr);
    gTrr1->SetPoint(pid,theta,sqrt(trr*trr+0.5*0.5));
    gTrr2->SetPoint(pid,theta,sqrt(trr*trr+1));
    pid++;
    
  }
  gSpr->Sort();
  gSep->Sort();
  gNph->Sort();
  gTrr->Sort();
  gTrr1->Sort();
  gTrr2->Sort();

  gSep->SetLineColor(38);
  gSpr->SetLineColor(38);
  gNph->SetLineColor(38);
  gTrr->SetLineColor(38);
  gTrr1->SetLineColor(2);
  gTrr2->SetLineColor(4);

  gTrr1->SetMarkerColor(2);
  gTrr2->SetMarkerColor(4);
  
  gSep->SetMarkerStyle(20);
  gSpr->SetMarkerStyle(20);
  gNph->SetMarkerStyle(20);
  gTrr->SetMarkerStyle(20);
  gTrr1->SetMarkerStyle(20);
  gTrr2->SetMarkerStyle(20);

  gSep->SetMarkerSize(0.7);
  gSpr->SetMarkerSize(0.7);
  gNph->SetMarkerSize(0.7);
  gTrr->SetMarkerSize(0.7);
  gTrr1->SetMarkerSize(0.7);
  gTrr2->SetMarkerSize(0.7);
	  
  gNph->GetYaxis()->SetRangeUser(0,150);
  gSpr->GetYaxis()->SetRangeUser(0,10);
  gTrr->GetYaxis()->SetRangeUser(0,2);

  gSpr->GetYaxis()->SetTitle("SPR [mrad]");
  gNph->GetYaxis()->SetTitle("multiplicity [#]");
  gTrr->GetYaxis()->SetTitle("#sigma_{#theta_{C} tr} [mrad]");
  gSpr->GetYaxis()->SetTitle("separation [s.d.]");
  
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


  gSpr->GetXaxis()->SetTitle("#theta_{track} [deg]");
  gNph->GetXaxis()->SetTitle("#theta_{track} [deg]");
  gTrr->GetXaxis()->SetTitle("#theta_{track} [deg]");
  gSep->GetXaxis()->SetTitle("#theta_{track} [deg]");

  prt_canvasAdd("spr",800,400);
  gSpr->Draw("APL");

  prt_canvasAdd("nph",800,400);
  gNph->Draw("APL");
  
  prt_canvasAdd("trr",800,400);
  gTrr->Draw("APL");
  gTrr1->Draw("same PL");
  gTrr2->Draw("same PL");

  prt_canvasAdd("sep",800,400);
  gSep->Draw("APL");
  
  prt_canvasSave();  
}
