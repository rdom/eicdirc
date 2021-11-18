#if defined(__ACLIC__)
#include "PrtTools.h"
#else
R__LOAD_LIBRARY(../build/libPrt.so)
#endif

void reco_start_time(TString in = "hits.root", TString pdf = "hits.pdf.root", double timeres = 0.1,
                     int pid = 321, TString nameid = "0", double sigma = 0, double var2 = 0) {

  PrtTools t(in);
  TCanvas *cc = new TCanvas("cc", "cc");
  TH1F *hllf = new TH1F("hllf", "hllf;ln L(K) - ln L(#pi); entries [#]", 180, -100, 100);
  TH1F *hlls = new TH1F("hlls", "hlls;ln L(K) - ln L(#pi); entries [#]", 180, -100, 100);
  TH1F *hl1 = new TH1F("hl1", "pdf;LE time [ns]; entries [#]", 2000, 0, 100);
  TH1F *hl2 = new TH1F("hl2", "pdf;LE time [ns]; entries [#]", 2000, 0, 100);
  TH1F *hl3 = new TH1F("hl3", "pdf;LE time [ns]; entries [#]", 2000, 0, 100);
  TGraph *glhf = new TGraph();
  TGraph *glhs = new TGraph();
  TH1F *hlhf = new TH1F("hlhf", "hllf;ln L(K) - ln L(#pi); entries [#]", 200, -5, 5);
  TH1F *hlhs = new TH1F("hlhs", "hlls;ln L(K) - ln L(#pi); entries [#]", 200, -5, 5);
  TH1F *hnph = new TH1F("hnph", ";multiplicity [#]; entries [#]", 200, 0, 200);
  TH1F *hmeanf = new TH1F("hmeanf", ";reconstructed time shift [ns]; entries [#]", 200, -0.7, 0.7);
  TH1F *hmeans = new TH1F("hmeans", ";reconstructed time shift [ns]; entries [#]", 200, -0.7, 0.7);
  hlhf->GetXaxis()->SetRangeUser(-4, 4);

  const int nch(8000);
  TF1 *pdff[nch], *pdfs[nch];
  TH1F *hpdff[nch], *hpdfs[nch];
  TFile f(pdf);

  int rebin = timeres / (100 / 2000.);
  int nnch = t.maxdircch();
 
  if (rebin > 0) hl3->Rebin(rebin);
  int integ1(0), integ2(0);
  for (int i = 0; i < nnch; i++) {
    hpdff[i] = (TH1F *)f.Get(Form("h_2_%d", i));
    hpdfs[i] = (TH1F *)f.Get(Form("h_3_%d", i));
    if (rebin > 0) hpdff[i]->Rebin(rebin);
    if (rebin > 0) hpdfs[i]->Rebin(rebin);
    integ1 += hpdff[i]->Integral();
    integ2 += hpdfs[i]->Integral();
    hpdff[i]->Smooth(1);
    hpdfs[i]->Smooth(1);
    hl3->Add(hpdff[i]);
    hl3->Add(hpdfs[i]);
  }

  TVirtualFitter *fitter;
  double time;
  int totalf(0), totals(0), ch(0);
  double noise = 1e-4;
  double theta = t.run()->getTheta();
  
  while (t.next() && t.i() < 5000) {
    double aminf, amins, sum(0), sumf(0), sums(0);
    int pdg = t.event()->getPid();
    if (pdg == pid && hllf->GetEntries() > 1800) continue;
    if (pdg == 211 && hlls->GetEntries() > 1800) continue;
    hnph->Fill(t.event()->getHits().size());

    int pp = 0;
    for (double s = -5; s < 5; s = s + 0.2) {
      sumf = 0;
      sums = 0;
      for (auto hit : t.event()->getHits()) {
	ch = hit.getChannel();
        time = hit.getLeadTime() + gRandom->Gaus(0, timeres) + s;

        aminf = hpdff[ch]->GetBinContent(hpdff[ch]->FindBin(time));
        amins = hpdfs[ch]->GetBinContent(hpdfs[ch]->FindBin(time));

        sumf += TMath::Log((aminf + noise));
        sums += TMath::Log((amins + noise));

        if (pdg == pid) hl1->Fill(time);
        if (pdg == 211) hl2->Fill(time);
      }
      // std::cout<<s<<" sum "<<sumf<<" "<<sums<<std::endl;
      glhf->SetPoint(pp, s, sumf);
      glhs->SetPoint(pp, s, sums);
      hlhf->Fill(s, sumf);
      hlhs->Fill(s, sums);
      pp++;
    }

    TF1 *flh = new TF1("flh", "gaus(0)+pol1(3)", -4, 4);
    flh->SetParameters(1, 0, 0.3, 1, 1);
    flh->SetParameters(1, 0, 0.3, 1, 1);
    flh->SetParLimits(1, -0.8, 0.8);
    flh->SetParLimits(2, 0.3, 5);

    hlhf->Fit("flh", "", "", -4, 4);
    double mean = flh->GetParameter(1);
    hmeanf->Fill(mean);

    hlhs->Fit("flh", "", "", -4, 4);
    mean = flh->GetParameter(1);
    hmeans->Fill(mean);
    double err[100]={1};

    // // hlhf->SetError(err);
    // hlhf->SetMarkerStyle(20);
    // hlhf->Draw("");
    // gPad->Update();
    // gPad->WaitPrimitive();

    hlhf->Reset();
    hlhs->Reset();

    sum = sumf - sums;
    if (fabs(sum) < 0.1) continue;

    if (pdg == pid) hllf->Fill(sum);
    if (pdg == 211) hlls->Fill(sum);
  }

  TString name = Form("%d_%1.2f", (int)theta, timeres);
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

  t.add_canvas("nph_" + name, 800, 400);
  hnph->Draw();
  double nph = t.fit(hnph, 50, 20, 50).X();

  t.add_canvas("mean_" + name, 800, 400);
  // sigma = t.fit(hmeanf, 5, 20, 5).Y();
  // var2 = t.fit(hmeans, 5, 20, 5).Y();
  double var3 = t.fit(hmeanf, 5, 20, 5).X();
  double var4 = t.fit(hmeans, 5, 20, 5).X();
  double rms = hmeanf->GetStdDev();
  
  TFitResultPtr rx = hmeanf->Fit("gaus", "S", "");
  sigma = rx->Parameter(2);
  rx = hmeans->Fit("gaus", "S", "");
  var2 = rx->Parameter(2);
  std::cout << sigma << " var " << var2 << std::endl;

  hmeanf->SetLineColor(kRed + 1);
  hmeanf->Draw();
  hmeans->SetLineColor(kBlue + 1);
  hmeans->Draw("same");

  t.add_canvas("sep_" + name, 600, 500);
  t.normalize(hllf, hlls);

  TF1 *ff;
  double m1(0), m2(0), s1(100), s2(100);
  if (hllf->GetEntries() > 10) {
    hllf->Fit("gaus", "Sq");
    ff = hllf->GetFunction("gaus");
    ff->SetLineColor(1);
    m1 = ff->GetParameter(1);
    s1 = ff->GetParameter(2);
  }
  if (hlls->GetEntries() > 10) {
    hlls->Fit("gaus", "Sq");
    ff = hlls->GetFunction("gaus");
    ff->SetLineColor(1);
    m2 = ff->GetParameter(1);
    s2 = ff->GetParameter(2);
  }
  double sep = (fabs(m1 - m2)) / (0.5 * (s1 + s2));
  std::cout << in << " separation " << sep << std::endl;
  hllf->SetTitle(Form("#theta = %d       #sigma = %1.2f", (int)theta, sep));

  hllf->SetLineColor(2);
  hllf->Draw();
  hlls->SetLineColor(4);
  hlls->Draw("same");

  hl1->Scale(1 / hl1->GetMaximum());
  hl2->Scale(1 / hl2->GetMaximum());
  hl3->Scale(1 / hl3->GetMaximum());

  t.normalize(hl1, hl2);
  t.add_canvas("time_" + name, 800, 500);
  hl1->Draw();
  hl2->SetLineColor(4);
  hl2->Draw("same");
  hl3->SetLineColor(2);
  hl3->Draw("same");
  t.save_canvas("data/reco_start_time", 0);
 
  TFile fc("data/reco_start_time/res_" + name + ".root", "recreate");
  TTree *tc = new TTree("reco", "reco");
  tc->Branch("theta", &theta, "theta/D");
  tc->Branch("sep", &sep, "sep/D");
  tc->Branch("timeres", &timeres, "timeres/D");
  tc->Branch("nph", &nph, "nph/D");
  tc->Branch("sigma", &sigma, "sigma/D");
  tc->Branch("rms", &rms, "rms/D");
  tc->Branch("var2", &var2, "var2/D");
  tc->Branch("var3", &var3, "var3/D");
  tc->Branch("var4", &var4, "var4/D");
  tc->Fill();
  tc->Write();
}
