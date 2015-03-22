void da_scan(TString inFile = "r_spr.root"){

  TChain ch("dirc"); ch.Add(inFile);
  Double_t cangle,spr,trr,nph,par1,par2,par3,par4,par5,par6,par7,par8,theta,phi; 

  TGraph *gSpr = new TGraph();
  TGraph *gNph = new TGraph();
  TH2F * hResol = new TH2F("Resol2","EV inner offset vs. Theta vs. SPR",100,30,150,100,0,10);

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
    hResol->Fill(theta,spr);
    Int_t adi = (Int_t)theta;
    if(adi<30 && adi>150) continue;
    gSpr->SetPoint(i,theta,spr);
    gNph->SetPoint(i,theta,nph);
  }
  gSpr->Sort();
  gNph->Sort();
  gSpr->SetMarkerStyle(21);
  gSpr->SetMarkerSize(0.7);
  gNph->SetMarkerStyle(21);
  gNph->SetMarkerSize(0.7);
  gNph->GetYaxis()->SetRangeUser(0,150);
  gSpr->GetYaxis()->SetRangeUser(0,12);

  TCanvas* c1 = new TCanvas("c1","c1",600,900);
  gSpr->Draw("APL");
  TCanvas* c2 = new TCanvas("c2","c2",600,900);
  gNph->Draw("APL");

}
