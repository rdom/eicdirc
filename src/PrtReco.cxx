// -----------------------------------------
// PrtReco.cxx
// implementation of DIRC reconstruction methods 
// Created on: 27.07.2024
// Author: r.dzhygadlo at gsi.de
// -----------------------------------------

#include "PrtReco.h"

#ifdef AI
#include "cppflow/cppflow.h"
#endif

using std::cout;
using std::endl;


struct {
  double ca;
  double td;
  double tr;
} cor[28] = {0, 0, 0};

int glob_i(0);
TGraph glob_gr;

// -----   Default constructor   -------------------------------------------
PrtReco::PrtReco(TString infile, TString lutfile, TString pdffile, TString nnfile, int verbose) {
  fVerbose = verbose;  
  fCriticalAngle = asin(1.00028 / 1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
  fRingFit = false;
  
  fNNPath = nnfile;
  fChain = new TChain("data");
  fChain->Add(infile);
  fEvent = new PrtEvent();
  frun = PrtManager::Instance()->getRun();
  ft = PrtTools(frun);
  fmaxch = ft.maxdircch();
  fnpmt = frun->getNpmt();
  fnpix = frun->getNpix();
  fMethod = frun->getRunType();
  fCorrType = frun->getCorrection();
  fStudyId = frun->getStudy();
  fMomentum = frun->getMomentum();
  fRadiator = frun->getRadiator();
  fRadiatorL = frun->getRadiatorL();
  fPhysList = frun->getPhysList();
  fTimeRes = frun->getTimeSigma();  
  fTimeCut = frun->getTimeCut();

  if(fMethod == 20) {
    fRingFit = true;
    fMethod = 2;
  }

  int rpid = frun->getPid();
  fp2 = 2;                    // pi
  if (rpid == 10000) fp1 = 0; // e
  if (rpid == 10001) fp1 = 1; // mu
  if (rpid == 10003) fp1 = 3; // K
  if (rpid == 10004) fp1 = 4; // p
  if (rpid == 10005) {
    fp2 = 3;
    fp1 = 4;
  }

  fTimeProp = new TH1F("timeprop", ";measured time [ns];entries [#]", 1000, 0, 80);
  fTimeDiff = new TH2F("timediff", ";measured time [ns];t_{meas}-t_{calc} [ns]", 500, 0, 100, 150, -5, 5);
  fChRing = new TH2F("Time4", "4", 200, -1, 1, 200, -1, 1);
  fTrackAngle0 = new TH1F("fTrackAngle0", ";#Delta [mrad];entries [#]", 500, -20, 20);
  fTrackAngle1 = new TH1F("fTrackAngle1", ";#Delta [mrad];entries [#]", 500, -20, 20);
  fTrackAngle2 = new TH1F("fTrackAngle2", ";#Delta [mrad];entries [#]", 500, -20, 20);

  fFindTime = new TH1F("ft", ";t_{measured}-t_{calculated} [ns];entries [#]", 2000, -100, 100);
  fFindTimeRes = new TH1F("ftr", "ftr", 100, -2, 2);
  fdtt = new TH2F("dtt", ";t_{meas}-t_{calc} [ns];#theta_{l} [deg]", 1000, -2, 2, 1000, 0, 90);
  fdtl = new TH2F("dtl", ";t_{meas}-t_{calc} [ns];path length [m]", 1000, -2, 2, 1000, 0, 15);
  fdtp = new TH2F("dtp", ";#theta_{l} [deg];path length [m]", 1000, 0, 90, 1000, 0, 15);
  fhChromL = new TH2F("chroml", ";\\ (t_{measured}-t_{calculated})/L_{path};#theta_{C} [rad]", 100,
                      -0.0002, 0.0002, 100, -0.03, 0.03);

  fChain->SetBranchAddress("PrtEvent", &fEvent);

  fFit = new TF1("fgaus", "[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +[3]", 0.35, 0.9);
  fChromCor = new TF1("fchrom1", "x<87? 23.15+0.41837*x : 93.82 - 0.435507*x", 15, 165);
   
  fSpect = new TSpectrum(10);
  fnX1 = TVector3(1, 0, 0);
  fnY1 = TVector3(0, 1, 0);

  int col[] = {kRed + 1, kBlue + 1, kBlack};
  for (int i = 0; i < 3; i++) {
    fTimeDiffR[i] = new TH1F(Form("TimeDiff_%d", i), ";t_{measured}-t_{calculated} [ns];entries [#]",
                            500, -10, 10);
    fTimeDiffR[i]->SetLineColor(col[i]);
  }
  for (int h = 0; h < 5; h++) {
    hthetac[h] = new TH1F(Form("thetac_%d", h), ";#theta_{C} [rad];entries [#]", 200, 0.75, 0.9);
    hthetacd[h] =
      new TH1F(Form("thetacd_%d", h), ";#Delta#theta_{C} [mrad];entries [#]", 200, -60, 60);
    hnph_gr[h] = new TH1F(Form("nph_gr_%d", h), ";detected photons [#];entries [#]", 220, 0, 220);
    hnph_ti[h] = new TH1F(Form("nph_ti_%d", h), ";detected photons [#];entries [#]", 220, 0, 220);
    fFunc[h] = new TF1(Form("gaus_%d", h), "[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])", 0.7, 0.9);

    int c = ft.color(h);
    hthetac[h]->SetLineColor(c);
    hthetacd[h]->SetLineColor(c);
    hnph_gr[h]->SetLineColor(c);
    hnph_ti[h]->SetLineColor(c);
    fFunc[h]->SetLineColor(c);
  }

  for (int i = 0; i < 20; i++) {
    fFindTimeA[i] =
      new TH1F(Form("fta_%d", i), ";t_{measured}-t_{calculated} [ns];entries [#]", 1000, -10, 10);
  }

  for (int i = 0; i < 28; i++) {
    fPmt_a[i] = new TH1F(Form("pmt_a_%d", i), Form("pmt %d;#theta_{C} [rad];entries [#]", i), 150,
                         -0.05, 0.05);
    fPmt_td[i] = new TH1F(Form("pmt_td_%d", i), Form("pmt_td%d;t_{m}-t_{c} [ns]", i), 150, -5, 5);
    fPmt_tr[i] = new TH1F(Form("pmt_tr_%d", i), Form("pmt_tr%d;t_{m}-t_{c} [ns]", i), 150, -5, 5);
  }

  double range = 100;
  if (fMomentum < 5) range = 500;
  if (fMomentum < 3) range = 500;
  if (fMomentum < 0.6) range = 500;
  if (fMomentum > 8) range = 200;
  if (fMomentum > 12) range = 100;
  if (fp1 == 0) {
    range = 200;
    if (fMomentum >= 2) range = 100;
  }

  // // pixels scan tmp
  // if(pixels >= 16) range = 500;
  // if(frun->getLens() == 10) range = 1500;
  // double range_gr = range;
  // if(pixels <= 16) range_gr = 200;
  // if(pixels > 16) range_gr = 500;
  // if(pixels > 40) range_gr = 1000;
  // if(pixels > 20 && fPhysList > 0) range_gr = 1000;
  // if(fMomentum > 8) range = 500;
  // if(fMomentum > 8) range_gr = 200;

  for (int h = 0; h < 5; h++) {
    TString la = ";ln L(K) - ln L(#pi);entries [#]";
    fLnDiffGr[h] = new TH1F(Form("LnDiffGr_%d", h), la, 400, -range, range);
    fLnDiffTi[h] = new TH1F(Form("LnDiffTi_%d", h), la, 400, -range, range);
    fLnDiffNn[h] = new TH1F(Form("LnDiffNn_%d", h), la, 400, -range, range);
    fLnDiffGr[h]->SetLineColor(ft.color(h));
    fLnDiffTi[h]->SetLineColor(ft.color(h));
    fLnDiffNn[h]->SetLineColor(ft.color(h));
    fLnDiffTi[h]->SetMarkerColor(ft.color(h) + 1);

    if (fMethod == 4) {
      for (int i = 0; i < fmaxch; i++) {
        fTime[h][i] =
          new TH1F(Form("h_%d_%d", h, i), "pdf;LE time [ns]; entries [#]", 2000, 0, 100);
      }
    }

    fSigma[h] = 0.007;
    if (fp1 == 0) { // electron
      fSigma[0] = 0.005;
    }
  }
  
  // read lut
  if (!gSystem->AccessPathName(lutfile)) {
    std::cout << "-I- reading  " << lutfile << std::endl;
    fFile = new TFile(lutfile);
    fTree = (TTree *)fFile->Get("prtlut");   
    fLut = new TClonesArray("PrtLutNode");
    fTree->SetBranchAddress("LUT", &fLut);
    fTree->GetEntry(0);
    fGeomReco = true;
  } else {
    std::cout << "-I- lut file not found  " << lutfile << std::endl;
    fGeomReco = false;
  }
  
  fTimeImaging = (fMethod == 4) ? true : false;

  // read pdf
  fPdfPath = pdffile;
  if (fPdfPath == "") {
    fPdfPath = infile;
    fPdfPath.ReplaceAll(".root", ".pdf.root");
  }

  if (fMethod == 2) {
    if (!gSystem->AccessPathName(fPdfPath)) {
      std::cout << "-I- reading  " << fPdfPath << std::endl;
      TFile pdfFile(fPdfPath);
      int binfactor = (int)(fTimeRes * 1000 / 50. + 0.1);
      std::cout << "binfactor " << binfactor << std::endl;
      
      for (int h : {fp1, fp2}) {
        for (int i = 0; i < fmaxch; i++) {
          // auto hpdf = (TH1F *)pdfFile.Get(Form("h_%d_%d",h, i));
          fTime[h][i] = (TH1F *)pdfFile.Get(Form("h_%d_%d", h, i));
          fTime[h][i]->SetDirectory(0);
          if (fTimeRes > 0) fTime[h][i]->Rebin(binfactor);
          // if (sigma > 0) hpdf->Rebin(binfactor);
          // // hpdf->Smooth();
          // fPdf[h][i] = new TGraph(hpdf);
          // fPdf[h][i]->SetBit(TGraph::kIsSortedX);
          fTimeImaging = true;
        }
      }
    } else {
      std::cout << "-I- pdf file not found  " << fPdfPath << std::endl;    
      fTimeImaging = false;
    }
  }

  // per-pmt corrections
  fCorrFile = infile + "_corr.root";
  if (fCorrType > 0 && fCorrType != 5) {
    if (!gSystem->AccessPathName(fCorrFile)) {
      std::cout << "-I- reading  " << fCorrFile << std::endl;
      int pmt, lvl;
      double c_ca, c_td, c_tr, spr[5];
      TChain ch;
      ch.SetName("corrections");
      ch.Add(fCorrFile);
      ch.SetBranchAddress("lvl", &lvl);
      ch.SetBranchAddress("pmt", &pmt);
      ch.SetBranchAddress("c_ca", &c_ca);
      ch.SetBranchAddress("c_td", &c_td);
      ch.SetBranchAddress("c_tr", &c_tr);
      ch.SetBranchAddress("spr", &spr);

      std::cout << "lvl  pmt    c_ca    c_td    c_tr  spr_pi    spr_k" << std::endl;
      for (int i = 0; i < ch.GetEntries(); i++) {
        ch.GetEvent(i);
        cor[pmt].ca = (fabs(c_ca) < 0.017) ? c_ca : 0;
	if(fCorrType == 2){
	  cor[pmt].td = (fabs(c_td) < 0.5) ? c_td : 0;
	  cor[pmt].tr = (fabs(c_tr) < 0.5) ? c_tr : 0;
	}
        for (int h = 0; h < 5; h++) {
          fSigma[h] = 0.001 * fabs(spr[h]) * 0.95;
          if (fSigma[h] > 0.02) fSigma[h] = 0.01;
        }
        cout.precision(2);
        std::cout << std::fixed << std::setw(3) << lvl << std::setw(5) << pmt << std::setw(8)
                  << 1000 * cor[pmt].ca << std::setw(8) << cor[pmt].td << std::setw(8)
                  << cor[pmt].tr << std::setw(8) << 1000 * fSigma[2] << std::setw(8)
                  << 1000 * fSigma[3] << std::endl;
      }
      fCorrLevel = lvl;
      cout.precision(-1);
    } else {
      std::cout << "-I- corr file not found  " << fCorrFile << std::endl;
    }
  }

  // read neural network model
#ifdef AI
  fNNet = true;
  if (!gSystem->AccessPathName(fNNPath)) {
    std::cout << "-I- reading  " << fNNPath << std::endl;
    fNNmodel = new cppflow::model(fNNPath.Data());
    for(auto s : (*fNNmodel).get_operations() ){
      std::cout << "s " << s << std::endl;   
    }  
  } else {
    fNNet = false;
    std::cout << "-I- neural net model not found  " << fNNPath << std::endl;
  }
#endif

}

// -----   Destructor   ----------------------------------------------------
PrtReco::~PrtReco() {}

//-------------- Loop over tracks ------------------------------------------
void PrtReco::Run(int start, int end) {
  TVector3 dird, dir, momInBar(0, 0, 1), posInBar;
  double mom = fMomentum;
  int tofPid(fp1), distPid(0), likePid(0);
  int eff_total[4] = {0}, eff_nn[4] = {0};
  bool reflected = kFALSE;
  gStyle->SetOptFit(111);

  int nsEvents(0), barid(0);

  TString outFile = PrtManager::Instance()->getOutName();
  double cangle[5] = {0}, spr[5] = {0}, trr[5] = {0}, nph_gr[5] = {0}, nph_gr_err[5] = {0},
         nph_ti[5] = {0}, nph_ti_err[5] = {0}, par5(0), par6(0), ctimeRes(0), trackRes(0), test1(0),
         test2(0), test3(0), sep_gr(0), sep_gr_err(0), sep_ti(0), sep_ti_err(0), sep_nn(0),
         sep_nn_err(0), epi_rejection1(0), epi_rejection2(0), epi_rejection3(0), track_res0(0),
         track_res1(0), track_res2(0);

  ft.set_palette(1);
  ft.create_maps();
  ft.init_digi();

  test1 = frun->getTest1();
  test2 = frun->getTest2();
  test3 = frun->getTest3();
  double dark_noise = frun->getDarkNoise();
  trackRes = frun->getBeamSize();
  double trackingResTheta = frun->getTrackingResTheta();
  double trackingResPhi = frun->getTrackingResPhi();
  fPhi = frun->getPhi();
  fTheta = frun->getTheta();
  if (fPhi >= 990) barid = fPhi - 990;
  std::cout << "-I- tracking resulution: dtheta = " << trackingResTheta
            << " dphi = " << trackingResPhi << " time res " << fTimeRes << std::endl;

  int nEvents = fChain->GetEntries();
  if (end == 0) end = nEvents;

  int pdfstart = 5000;
  if (end > pdfstart) end = pdfstart;
  if (fMethod == 4) {
    start = pdfstart;
    pdfstart = nEvents;
    end = nEvents;
  }

  std::cout << "-I- run started for [" << start << "," << end << "]" << std::endl;

  for (int ievent = start; ievent < nEvents && ievent < pdfstart; ievent++) {
    fChain->GetEntry(ievent);
    // theta = (fEvent->getMomentum().Angle(TVector3(0, 0, -1))) * TMath::RadToDeg();

    int pid = fEvent->getPid();    


    if (ievent % 1000 == 0)
      std::cout << "event # " << ievent << " has " << fEvent->getHits().size() << " hits"
                << std::endl;
    double minChangle = 0.35;
    double maxChangle = 0.9;
    double m = 0.001 * fEvent->getMomentum().Mag();
    TVector3 mom_vertex = fEvent->getMomentum().Unit();
    TVector3 mom_before = fEvent->getMomentumBefore().Unit();
    TVector3 mom_after = fEvent->getMomentumAfter().Unit();

    double theta_diff0 = 1000 * (mom_vertex.Theta() - mom_before.Theta());
    double theta_diff2 = 1000 * (mom_vertex.Theta() - mom_after.Theta());
    fTrackAngle0->Fill(theta_diff0);
    fTrackAngle2->Fill(theta_diff2);

    // // post-dirc tracking layer
    // TVector3 pa = fEvent->getPositionAfter();
    // TVector3 ma = fEvent->getMomentumAfter().Unit();
    // TVector3 pb = fEvent->getPosition();

    // TVector3 diff = pa - (pa + TVector3(15, 0, 0));
    // double prod1 = diff.Dot(TVector3(1, 0, 0));
    // double prod2 = ma.Dot(TVector3(1, 0, 0));
    // TVector3 ppa = pa - ma * (prod1 / prod2);

    // double tsr = 0.05;
    // pb += TVector3(0, gRandom->Gaus(0, tsr), gRandom->Gaus(0, tsr));
    // ppa += TVector3(0, gRandom->Gaus(0, 0.05), gRandom->Gaus(0, 0.05));

    // rotatedmom = (ppa - pb);

    // track already smeared during simulation at tracking layer. rotatedmom is direction at vertex
    // rotatedmom.SetTheta(gRandom->Gaus(rotatedmom.Theta(), trackingResTheta));
    // rotatedmom.SetPhi(gRandom->Gaus(rotatedmom.Phi(), trackingResPhi));

    for (int i = 0; i < 5; i++) {
      fAngle[i] = acos(sqrt(m * m + ft.mass(i) * ft.mass(i)) / m / 1.4738); // 1.4738 = 370 = 3.35
      fFunc[i]->SetParameter(0, 1);
      fFunc[i]->SetParameter(1, fAngle[i]);
      fFunc[i]->SetParameter(2, fSigma[i]);
    }

    // double stime = FindStartTime(fEvent);
    geom_reco(fEvent, mom_vertex, fRingFit);

    if (fRingFit) {
      
      double x0(0), y0(0), a(fAngle[2]);
      a = 0.5 * (fAngle[fp1] + fAngle[fp2]);

      FitRing(x0, y0, a);
      
      TVector3 corr(x0, y0, 1 - TMath::Sqrt(x0 * x0 + y0 * y0));
      corr = corr.Unit();

      glob_i = 0;
      glob_gr.Set(0);
      
      if (0) {
        // ft.add_canvas(Form("ring_%d",ievent), 1200, 1200);
        fChRing->SetStats(0);
        fChRing->GetXaxis()->SetTitle("#theta_{c}sin(#varphi_{c})");
        fChRing->GetYaxis()->SetTitle("#theta_{c}cos(#varphi_{c})");
        fChRing->SetTitle(Form("%1.0f#circ polar angle", fTheta));
        fChRing->Draw("colz");
        // fChRing->SetMaximum(1);

        TLegend *leg = new TLegend(0.32, 0.42, 0.67, 0.59);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        // leg->AddEntry((TObject *)0, Form("Entries %0.0f", fChRing->GetEntries()), "");
        leg->AddEntry((TObject *)0, Form("#Delta#theta_{c} %1.2f [mrad]", corr.Theta() * 1000), "");
        leg->AddEntry((TObject *)0, Form("#Delta#varphi_{c} %1.2f [rad]", corr.Phi()), "");
        leg->Draw("same");

        TArc *arc = new TArc(x0, y0, fAngle[2]);
        arc->SetLineColor(kBlack);
        arc->SetLineWidth(2);
        arc->SetFillStyle(0);
        arc->Draw("same");
	    
        gPad->Update();
        gPad->WaitPrimitive();

        // ft.save_canvas("reco", 1, 0, 1);
        fChRing->Reset();
      }
      
      TVector3 mom_fitted = mom_vertex;
      mom_fitted.RotateY(-corr.Theta());
      mom_fitted.Rotate(corr.Phi() + TMath::PiOver2(), mom_vertex);

      // corr.RotateUz(uv);
      // corr.Print();
      // TVector3 oo = mom_vertex;
      // std::cout << "oo.Theta() " << oo.Theta() <<  " oo.Theta() " << corr.Theta() << std::endl;
      // oo.SetTheta(oo.Theta() + corr.Theta() );
      // // oo.RotateX(-corr.Theta());
      // // oo.Rotate(corr.Phi(), mom_vertex);
      // oo = corr;
      double theta_diff1 = 1000 * (mom_fitted.Theta() - mom_before.Theta());
      fTrackAngle1->Fill(theta_diff1);

      geom_reco(fEvent, mom_fitted);
    }

    if (fTimeImaging) time_imaging(fEvent);
    if (fNNet) nn_reco(fEvent, mom_vertex);

    double sum_nph = 0;
    // if (m < 2.5) {  // photon yield likelihood
    //   TF1 *f_pi = new TF1("gaus", "gaus", 0, 150);
    //   f_pi->SetParameters(1, 50, 8.7); // fix me
    //   TF1 *f_p = new TF1("gaus", "gaus", 0, 150);
    //   f_p->SetParameters(1, 50, 7.1);

    //   double lh_nph_p = f_p->Eval(fNph_ti[pid]);
    //   double lh_nph_pi = f_pi->Eval(fNph_ti[pid]);
    //   sum_nph = lh_nph_p - lh_nph_pi;
    // }

    if (fVerbose == 1) {
      ft.add_canvas("ff", 800, 400);
      if (hthetac[fp1]->GetMaximum() > 0) hthetac[fp1]->Scale(1 / hthetac[fp1]->GetMaximum());
      hthetac[fp1]->Draw("hist");
      fFunc[fp1]->Draw("same");
      fFunc[fp2]->Draw("same");
      ft.wait_primitive("ff");
      // ft.canvasDel("ff");

      FindPeak(cangle, spr);
      hthetac[fp1]->Reset();
    }

    if (++nsEvents >= end) break;
  }

  { // calclulate efficiencies
    // double eff = ft.calculate_efficiency(fLnDiffGr[fp1],fLnDiffGr[fp2]);
    // std::cout << "GR eff = " << eff << std::endl;
    // if (eff_total[3] > 0) std::cout << "NN eff = " << eff_nn[3] / (float) eff_total[3] << std::endl;
    // std::cout << "eff_total " << eff_total[3] <<  " eff_nn " << eff_nn[3] << std::endl;

    double eff_gr = ft.calculate_efficiency(fLnDiffGr[fp1], fLnDiffGr[fp2]);
    std::cout << "Eff GR = " << eff_gr << std::endl;
    double eff_ti = ft.calculate_efficiency(fLnDiffTi[fp1], fLnDiffTi[fp2]);
    std::cout << "Eff TI = " << eff_ti << std::endl;
    double eff_nnt = ft.calculate_efficiency(fLnDiffNn[fp1], fLnDiffNn[fp2]);
    std::cout << "Eff NN = " << eff_nnt << std::endl;

    if (eff_total[2] > 0) std::cout << "Eff NN (pi) = " << eff_nn[2] / (float)eff_total[2] << std::endl;
    if (eff_total[3] > 0) std::cout << "Eff NN (K) = " << eff_nn[3] / (float)eff_total[3] << std::endl;
  }
  
  if (fMethod == 4) { // create pdf
    std::cout << "saving pdfs into " << fPdfPath << std::endl;

    TFile efile(fPdfPath, "RECREATE");
    for (int h : {fp1, fp2}) {
      for (int i = 0; i < fmaxch; i++) {
        fTime[h][i]->Scale(1 / (double)fTotal[h]);
        fTime[h][i]->Write();
      }
    }
    efile.Write();
    efile.Close();
    std::cout << "totals: " << fTotal[fp1] << " " << fTotal[fp2] << std::endl;
  }

  if (fMethod == 2) {
    FindPeak(cangle, spr);
 
    for (int h = 0; h < 5; h++) {
      if (hnph_gr[h]->Integral() < 20) continue;
      hnph_gr[h]->Fit("gaus", "SQ", "", 5, 250);
      auto f = hnph_gr[h]->GetFunction("gaus");
      if (f) {
        nph_gr[h] = f->GetParameter(1);
        nph_gr_err[h] = f->GetParError(1);
      }

      if (hnph_ti[h]->Integral() < 20) continue;
      hnph_ti[h]->Fit("gaus", "SQ", "", 5, 250);
      f = hnph_ti[h]->GetFunction("gaus");
      if (f) {
        nph_ti[h] = f->GetParameter(1);
        nph_ti_err[h] = f->GetParError(1);
      }
    }

    TF1 *ff;
    double m1 = 0, m2 = 0, s1 = 100, s2 = 100, dm1 = 0, dm2 = 0, ds1 = 0, ds2 = 0;
    if (fLnDiffGr[fp2]->Integral() > 10) {
      fLnDiffGr[fp2]->Fit("gaus", "SQ");
      ff = fLnDiffGr[fp2]->GetFunction("gaus");
      if (ff) {
        m1 = ff->GetParameter(1);
        s1 = ff->GetParameter(2);
        dm1 = ff->GetParError(1);
        ds1 = ff->GetParError(2);
      }
      if (fp1 == 0 && mom < 1.5) { // handle tails
        fLnDiffGr[fp2]->Fit("gaus", "SQ", "", m1 - 2.0 * s1, 500);
        ff = fLnDiffGr[fp2]->GetFunction("gaus");
        if (ff) {
          m1 = ff->GetParameter(1);
          s1 = ff->GetParameter(2);
          dm1 = ff->GetParError(1);
          ds1 = ff->GetParError(2);
        }
      }
    }
    if (fLnDiffGr[fp1]->Integral() > 10) {
      fLnDiffGr[fp1]->Fit("gaus", "SQ");
      ff = fLnDiffGr[fp1]->GetFunction("gaus");
      if (ff) {
        m2 = ff->GetParameter(1);
        s2 = ff->GetParameter(2);
        dm2 = ff->GetParError(1);
        ds2 = ff->GetParError(2);
      }
      if (fp1 == 0 && mom < 1.5) { /// handle tails
        fLnDiffGr[fp1]->Fit("gaus", "SQ", "", -500, m2 + 2.0 * s2);
        ff = fLnDiffGr[fp1]->GetFunction("gaus");
        if (ff) {
          m2 = ff->GetParameter(1);
          s2 = ff->GetParameter(2);
          dm2 = ff->GetParError(1);
          ds2 = ff->GetParError(2);
        }
      }
    }
    sep_gr = (fabs(m1 - m2)) / (0.5 * (s1 + s2));

    double e1, e2, e3, e4;
    e1 = 2 / (s1 + s2) * dm1;
    e2 = -2 / (s1 + s2) * dm2;
    e3 = -2 * (m1 - m2) / ((s1 + s2) * (s1 + s2)) * ds1;
    e4 = -2 * (m1 - m2) / ((s1 + s2) * (s1 + s2)) * ds2;
    sep_gr_err = sqrt(e1 * e1 + e2 * e2 + e3 * e3 + e4 * e4);

    if (fTimeImaging) {
      m1 = 0, m2 = 0;

      epi_rejection1 = CalcRejection(fLnDiffTi[fp1], fLnDiffTi[fp2], 0.90);
      epi_rejection2 = CalcRejection(fLnDiffTi[fp1], fLnDiffTi[fp2], 0.95);
      epi_rejection3 = CalcRejection(fLnDiffTi[fp1], fLnDiffTi[fp2], 0.98);

      if (fLnDiffTi[fp2]->Integral() > 10) {
        fLnDiffTi[fp2]->Fit("gaus", "Q");
        ff = fLnDiffTi[fp2]->GetFunction("gaus");
        if (ff) {
          m1 = ff->GetParameter(1);
          s1 = ff->GetParameter(2);
          dm1 = ff->GetParError(1);
          ds1 = ff->GetParError(2);
          ff->SetLineColor(kBlack);
        }

        if (fp1 == 0 && mom < 1.5) { /// handle tails
          fLnDiffTi[fp2]->Fit("gaus", "S", "", -500, m1 + 1.5 * s1);
          ff = fLnDiffTi[fp2]->GetFunction("gaus");
          m2 = ff->GetParameter(1);
          s2 = ff->GetParameter(2);
          dm2 = ff->GetParError(1);
          ds2 = ff->GetParError(2);
        }
      }

      if (fLnDiffTi[fp1]->Integral() > 10) {
        fLnDiffTi[fp1]->Fit("gaus", "Q");
        ff = fLnDiffTi[fp1]->GetFunction("gaus");
        if (ff) {
          m2 = ff->GetParameter(1);
          s2 = ff->GetParameter(2);
          dm2 = ff->GetParError(1);
          ds2 = ff->GetParError(2);
          ff->SetLineColor(kBlack);
        }
        if (fp1 == 0 && mom < 1.5) { /// handle tails
          fLnDiffTi[fp1]->Fit("gaus", "S", "", m2 - 1.5 * s2, 500);
          ff = fLnDiffTi[fp1]->GetFunction("gaus");
          m2 = ff->GetParameter(1);
          s2 = ff->GetParameter(2);
          dm2 = ff->GetParError(1);
          ds2 = ff->GetParError(2);
        }
      }

      // double epid = 0, pimisid = 0;
      // {
      //   auto f1 = fLnDiffTi[fp1]->GetFunction("gaus");
      //   auto f2 = fLnDiffTi[fp2]->GetFunction("gaus");
      //   double mid = 0; // 0.5 * (f1->GetParameter(1) + f2->GetParameter(1));
      //   epid = f1->Integral(mid, 300);
      //   pimisid = f2->Integral(mid, 300);
      //   epi_rejection = epid / pimisid;
      //   std::cout << "F " << epid << " " << pimisid << std::endl;
      // }

      std::cout << "rejection " << epi_rejection1 << std::endl;
      std::cout << "rejection " << epi_rejection2 << std::endl;
      std::cout << "rejection " << epi_rejection3 << std::endl;

      sep_ti = (fabs(m1 - m2)) / (0.5 * (s1 + s2));

      e1 = 2 / (s1 + s2) * dm1;
      e2 = -2 / (s1 + s2) * dm2;
      e3 = -2 * (m1 - m2) / ((s1 + s2) * (s1 + s2)) * ds1;
      e4 = -2 * (m1 - m2) / ((s1 + s2) * (s1 + s2)) * ds2;
      sep_ti_err = sqrt(e1 * e1 + e2 * e2 + e3 * e3 + e4 * e4);
    }

    if (1) {
      m1 = 0;
      m2 = 0;
      dm1 = 0;
      dm2 = 0;
      if (fLnDiffNn[fp2]->Integral() > 10) {
        fLnDiffNn[fp2]->Fit("gaus", "Q");
        ff = fLnDiffNn[fp2]->GetFunction("gaus");
        if (ff) {
          m1 = ff->GetParameter(1);
          s1 = ff->GetParameter(2);
          dm1 = ff->GetParError(1);
          ds1 = ff->GetParError(2);
          ff->SetLineColor(kBlack);
        }
      }

      if (fLnDiffNn[fp1]->Integral() > 10) {
        fLnDiffNn[fp1]->Fit("gaus", "Q");
        ff = fLnDiffNn[fp1]->GetFunction("gaus");
        if (ff) {
          m2 = ff->GetParameter(1);
          s2 = ff->GetParameter(2);
          dm2 = ff->GetParError(1);
          ds2 = ff->GetParError(2);
          ff->SetLineColor(kBlack);
        }
      }

      sep_nn = (fabs(m1 - m2)) / (0.5 * (s1 + s2));

      e1 = 2 / (s1 + s2) * dm1;
      e2 = -2 / (s1 + s2) * dm2;
      e3 = -2 * (m1 - m2) / ((s1 + s2) * (s1 + s2)) * ds1;
      e4 = -2 * (m1 - m2) / ((s1 + s2) * (s1 + s2)) * ds2;
      sep_nn_err = sqrt(e1 * e1 + e2 * e2 + e3 * e3 + e4 * e4);
    }

    std::cout << Form("%3d : SPR = %2.2f N_gr = %2.2f +/- %2.2f  N_ti = %2.2f +/- %2.2f",
                      ft.pdg(fp1), spr[fp1], nph_gr[fp1], nph_gr_err[fp1], nph_ti[fp1],
                      nph_ti_err[fp1])
              << std::endl;
    std::cout << Form("%3d : SPR = %2.2f N_gr = %2.2f +/- %2.2f  N_ti = %2.2f +/- %2.2f",
                      ft.pdg(fp2), spr[fp2], nph_gr[fp2], nph_gr_err[fp2], nph_ti[fp2],
                      nph_ti_err[fp2])
              << std::endl;
    std::cout << Form("SEP GR = %2.2f +/- %2.2f ", sep_gr, sep_gr_err) << std::endl;
    std::cout << Form("SEP TI = %2.2f +/- %2.2f ", sep_ti, sep_ti_err) << std::endl;
  }

  if (!fVerbose) gROOT->SetBatch(1);

  if (0) { // draw start time
    ft.add_canvas(Form("ctimeres_%d", int(fTheta + 0.01)), 800, 400);
    fFindTimeRes->Draw();

    double rr[20];
    for (int i = 0; i < 20; i++) {
      TGaxis::SetMaxDigits(3);
      ft.add_canvas(Form("cta_%d", i), 800, 400);
      rr[i] = ft.fit(fFindTimeA[i], 6, 20, 4, 1, 0).Y();
      fFindTimeA[i]->Draw();
    }
    for (int i = 0; i < 20; i++) {
      std::cout << (i ? "," : "") << rr[i] << std::endl;
    }
    ctimeRes = ft.fit(fFindTimeRes).Y();

    gStyle->SetOptStat(0);
    ft.add_canvas("fdtt", 800, 500);
    fdtt->Draw("colz");
    fdtt->SetMaximum(0.8 * fdtt->GetMaximum());

    ft.add_canvas("fdtl", 800, 500);
    fdtl->Draw("colz");
    fdtl->SetMaximum(0.8 * fdtl->GetMaximum());

    // prt_fitslices(fdtl,-2,2,2,2,0)->Draw("pl same");
    // prt_fitslices(fdtl,-2,2,2,2,2)->Draw("pl same");
    // prt_fitslices(fdtl,-2,2,2,2,3)->Draw("pl same");

    ft.add_canvas("fdtp", 800, 500);
    fdtp->Draw("colz");
    fdtp->SetMaximum(0.8 * fdtp->GetMaximum());

    if (fVerbose > 1) {
      gPad->Modified();
      gPad->Update();
      gPad->WaitPrimitive();
    }
  }

  { // track resolution with cherenkov ring fit
    auto r = fTrackAngle0->Fit("gaus","SQ");
    if (r > -1) track_res0 = r->Parameter(2);
    TF1 *mgaus = new TF1("mgaus","gaus");
    mgaus->SetParLimits(1,-0.7,0.7);
    r = fTrackAngle1->Fit("mgaus","SQ");
    if (r > -1) track_res1 = r->Parameter(2);
    r = fTrackAngle2->Fit("gaus","SQ");
    if (r > -1) track_res2 = r->Parameter(2);
    std::cout << "track_res0 " << track_res0 <<  " track_res1 " << track_res1 <<" track_res2 " << track_res2 << std::endl;
  }

  { // chromatic corrections 
    fhChromL->SetStats(0);
    auto g = ft.fit_slices_x(fhChromL, -0.00008, 0.00008, 0.01, 2, 0);
    g->SetMarkerStyle(20);
    g->SetLineColor(kBlack);
    g->SetMarkerSize(1.5);
    auto r = g->Fit("pol1", "S");
    if (r > -1) test3 = r->Parameter(1);
    std::cout << "test3 " << test3 << std::endl;
    
    fhChromL->Draw("colz");
    g->Draw("PL");
  }

  { // tree
    // outFile.ReplaceAll("reco_", Form("reco_%d_", frun->getId()));
    TFile file(outFile, "recreate");
    TTree tree("reco", "reco");
    tree.Branch("mom", &mom, "mom/D");
    tree.Branch("theta", &fTheta, "fTheta/D");
    tree.Branch("phi", &fPhi, "fPhi/D");
    tree.Branch("barid", &barid, "barid/I");
    tree.Branch("tofPid", &tofPid, "tofPid/I");
    tree.Branch("distPid", &distPid, "distPid/I");
    tree.Branch("likePid", &likePid, "likePid/I");
    tree.Branch("spr", &spr, "spr[5]/D");
    tree.Branch("trr", &trr, "trr[5]/D");
    tree.Branch("nph_gr", &nph_gr, "nph_gr[5]/D");
    tree.Branch("nph_gr_err", &nph_gr_err, "nph_gr_err[5]/D");
    tree.Branch("nph_ti", &nph_ti, "nph_ti[5]/D");
    tree.Branch("nph_ti_err", &nph_ti_err, "nph_ti_err[5]/D");
    tree.Branch("cangle", &cangle, "cangle[5]/D");
    tree.Branch("sep_gr", &sep_gr, "sep_gr/D");
    tree.Branch("sep_gr_err", &sep_gr_err, "sep_gr_err/D");
    tree.Branch("sep_ti", &sep_ti, "sep_ti/D");
    tree.Branch("sep_ti_err", &sep_ti_err, "sep_ti_err/D");
    tree.Branch("sep_nn", &sep_nn, "sep_nn/D");
    tree.Branch("sep_nn_err", &sep_nn_err, "sep_nn_err/D");
    tree.Branch("track_res0", &track_res0, "track_res0/D");
    tree.Branch("track_res1", &track_res1, "track_res1/D");
    tree.Branch("track_res2", &track_res2, "track_res2/D");

    tree.Branch("trackres", &trackRes, "trackRes/D");
    tree.Branch("trackingResTheta", &trackingResTheta, "trackingResTheta/D");
    tree.Branch("timeres", &fTimeRes, "fTimeRes/D");
    tree.Branch("timecut", &fTimeCut, "fTimeCut/D");
    tree.Branch("ctimeres", &ctimeRes, "ctimeRes/D");
    tree.Branch("test1", &test1, "test1/D");
    tree.Branch("test2", &test2, "test2/D");
    tree.Branch("test3", &test3, "test3/D");
    tree.Branch("dark_noise", &dark_noise, "dark_noise/D");
    tree.Branch("par5", &par5, "par5/D");
    tree.Branch("par6", &par6, "par6/D");
    tree.Branch("epi_rejection1", &epi_rejection1, "epi_rejection1/D");
    tree.Branch("epi_rejection2", &epi_rejection2, "epi_rejection2/D");
    tree.Branch("epi_rejection3", &epi_rejection3, "epi_rejection3/D");

    tree.Fill();
    tree.Write();
    std::cout << "-I- writing " << outFile << std::endl;
  }

  if (fVerbose > 1) {
    TString nid = Form("_%d_%1.2f_%1.4f_%1.2f", fp1, frun->getTheta(), test1, mom);
    TGaxis::SetMaxDigits(3);

    { // cherenkov angle
      ft.add_canvas("tangle" + nid, 800, 400);
      ft.normalize(hthetac, 5);

      hthetac[fp1]->SetTitle(Form("theta %1.2f", fTheta));
      hthetac[fp1]->Draw("");
      hthetac[fp2]->Draw("same");
      drawTheoryLines(mom);

      // ft.add_canvas("tangled" + nid, 800, 400);
      // ft.normalize(hthetacd, 5);
      // hthetacd[fp2]->SetTitle(Form("theta %1.2f", fTheta));
      // hthetacd[fp2]->Draw("");
      // hthetacd[fp1]->Draw("same");
    }

    { // nph
      ft.add_canvas("nph" + nid, 800, 400);
      ft.normalize(hnph_gr, 5);
      hnph_gr[fp1]->SetStats(0);
      hnph_gr[fp1]->Draw();
      hnph_gr[fp2]->Draw("same");

      hnph_ti[fp1]->Draw("same");
      hnph_ti[fp2]->Draw("same");
    }

    { // sep
      fLnDiffGr[fp1]->SetStats(0);
      fLnDiffGr[fp2]->SetStats(0);
      fLnDiffTi[fp1]->SetStats(0);
      fLnDiffTi[fp2]->SetStats(0);
      fLnDiffNn[fp1]->SetStats(0);
      fLnDiffNn[fp2]->SetStats(0);

      ft.add_canvas("lh_gr" + nid, 800, 400);
      ft.normalize(fLnDiffGr, 5);
      fLnDiffGr[fp2]->SetName(Form("s_%2.2f", sep_gr));
      fLnDiffGr[fp2]->SetTitle(Form("GR separation = %2.2f s.d.", sep_gr));
      TString lhtitle = "ln L(" + ft.lname(fp1) + ") - ln L(" + ft.lname(fp2) + ")";
      fLnDiffGr[fp2]->GetXaxis()->SetTitle(lhtitle);

      fLnDiffGr[fp2]->Draw();
      fLnDiffGr[fp1]->Draw("same");

      if (fTimeImaging) {
        fLnDiffTi[fp2]->GetXaxis()->SetTitle(lhtitle);
        ft.add_canvas("lh_ti" + nid, 800, 400);
        ft.normalize(fLnDiffTi[fp1], fLnDiffTi[fp2]);
        fLnDiffTi[fp2]->SetTitle(Form("TI separation = %2.2f s.d.", sep_ti));
        // fLnDiffTi[fp2]->SetLineColor(kBlue + 1);
        fLnDiffTi[fp2]->SetMarkerStyle(20);
        fLnDiffTi[fp2]->SetMarkerSize(0.85);
        // fLnDiffTi[fp2]->SetMarkerColor(kBlue + 1);
        fLnDiffTi[fp2]->Draw("E");

        // fLnDiffTi[fp1]->SetLineColor(kRed + 1);
        fLnDiffTi[fp1]->SetMarkerStyle(20);
        fLnDiffTi[fp1]->SetMarkerSize(0.85);
        // fLnDiffTi[fp1]->SetMarkerColor(kRed + 1);
        fLnDiffTi[fp1]->SetName(Form("s_%2.2f", sep_ti));
        fLnDiffTi[fp1]->Draw("E same");
      }
#ifdef AI
      ft.add_canvas("lh_nn" + nid, 800, 400);
      fLnDiffNn[fp2]->SetTitle(Form("NN separation = %2.2f s.d.", sep_nn));
      fLnDiffNn[fp2]->GetXaxis()->SetTitle(lhtitle);
      fLnDiffNn[fp2]->Draw();
      fLnDiffNn[fp1]->Draw("same");
#endif

    }

    { // chromatic corrections
      ft.add_canvas("chroml" + nid, 800, 400);
      fhChromL->SetStats(0);
      auto g = ft.fit_slices_x(fhChromL,-0.00008,0.00008,0.03,2,0);      
      g->SetMarkerStyle(20);
      g->SetLineColor(kBlack);
      g->SetMarkerSize(1.5);
      g->Fit("pol1");
      fhChromL->Draw("colz");
      g->Draw("PL");
    }

    { // hp
      auto cdigi = ft.draw_digi(0, 0);
      cdigi->SetName("hp" + nid);
      ft.add_canvas(cdigi);
    }

    if (fRingFit) { // cherenkov ring

      ft.add_canvas("ring" + nid, 800, 800);

      fChRing->SetStats(0);
      fChRing->GetXaxis()->SetTitle("#theta_{c}sin(#varphi_{c})");
      fChRing->GetYaxis()->SetTitle("#theta_{c}cos(#varphi_{c})");
      fChRing->SetTitle(Form("%1.0f#circ polar angle", fTheta));
      fChRing->Draw("colz");
      double x0(0), y0(0), a = 0.5 * (fAngle[2] + fAngle[3]);
      FitRing(x0, y0, a);
      TVector3 corr(x0, y0, 1 - TMath::Sqrt(x0 * x0 + y0 * y0));

      TLegend *leg = new TLegend(0.32, 0.42, 0.67, 0.59);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry((TObject *)0, Form("Entries %0.0f", fChRing->GetEntries()), "");
      leg->AddEntry((TObject *)0, Form("#Delta#theta_{c} %1.2f [mrad]", corr.Theta() * 1000), "");
      leg->AddEntry((TObject *)0, Form("#Delta#varphi_{c} %1.2f [rad]", corr.Phi()), "");
      leg->Draw();

      TArc *arc = new TArc(x0, y0, fAngle[2]);
      arc->SetLineColor(kBlack);
      arc->SetLineWidth(2);
      arc->SetFillStyle(0);
      arc->Draw();
      glob_i = 0;
      glob_gr.Set(0);
    }

    if (fCorrType > 0 && fCorrLevel < 1) { // corrections
      std::cout << "-I- writing " << fCorrFile << std::endl;

      TFile fc(fCorrFile, "recreate");
      TTree *tc = new TTree("corrections", "corrections");
      int pmt, lvl = 0;
      double c_ca = 0, c_td = 0, c_tr = 0, fitrange = 0.01;
      tc->Branch("lvl", &lvl, "lvl/I");
      tc->Branch("pmt", &pmt, "pmt/I");
      tc->Branch("c_ca", &c_ca, "c_ca/D");
      tc->Branch("c_td", &c_td, "c_td/D");
      tc->Branch("c_tr", &c_tr, "c_tr/D");
      tc->Branch("spr", &spr, "spr[5]/D");

      TString nid = Form("cor_a%d", lvl);
      ft.add_canvas(nid, 1600, 1000);
      auto c = ft.get_canvas(nid);
      c->Divide(6, 4);
      nid = Form("cor_t%d", lvl);
      ft.add_canvas(nid, 1600, 1000);
      auto c2 = ft.get_canvas(nid);
      c2->Divide(6, 4);
      c->cd(1);
      for (pmt = 0; pmt < fnpmt; pmt++) {
        if (fPmt_a[pmt]->GetEntries() < 20) continue;

	  c->cd(pmt + 1);
	  c_td = -ft.fit(fPmt_td[pmt], 5, 20, 1).X();
	  c_tr = -ft.fit(fPmt_tr[pmt], 5, 20, 1).X();

	  fPmt_td[pmt]->Draw();
	  fPmt_tr[pmt]->SetLineColor(2);
          fPmt_tr[pmt]->Draw("same");

	  c2->cd(pmt + 1);
          c_ca = -ft.fit(fPmt_a[pmt], fitrange, 20, 0.01).X();
	  lvl = 2;
          tc->Fill();
          fPmt_a[pmt]->Draw();
	  std::cout <<"c_ca " << c_ca << "c_td " << c_td << " c_tr " << c_tr << std::endl;
      }

      tc->Write();
      fc.Write();
      fc.Close();
    }

    { // time
      ft.add_canvas("tdiff" + nid, 800, 400);
      ft.normalize(fTimeDiffR, 3);
      for (int i = 0; i < 3; i++) {
        if (fTimeDiffR[i]->GetEntries() > 100) {
          fTimeDiffR[i]->SetStats(0);
          fTimeDiffR[i]->SetTitle(Form("theta %1.2f", fTheta));
          fTimeDiffR[i]->Draw((i == 0) ? "h" : "hsame");
        }
      }

      ft.add_canvas("diff" + nid, 800, 400);
      fTimeDiff->SetStats(0);
      fTimeDiff->Draw("colz");
      
      ft.add_canvas("time" + nid, 800, 400);
      fTimeProp->Draw();
    }

    { // track smearing
      ft.add_canvas("track_smear" + nid, 800, 400);
      ft.normalize(fTrackAngle0, fTrackAngle1);     
      fTrackAngle1->SetLineColor(kRed);
      fTrackAngle0->GetXaxis()->SetRangeUser(-10, 10);
      fTrackAngle0->Draw();
      fTrackAngle1->Draw("same");
    }

    TString filedir = fCorrFile;
    if (filedir.Contains("/")) {
      filedir.Remove(filedir.Last('/'));
      filedir += "/";
    } else filedir = "";

    ft.save_canvas(filedir + "reco", 0, 0, 0);
    
    if (fVerbose > 2) ft.wait_primitive("lh_gr" + nid, "none");
  }

  // delete fTime; // abort now to save time (for small pixels)
}

void PrtReco::FindPeak(double (&cangle)[5], double (&spr)[5]) {
  for (int h = 0; h < 5; h++) {
    spr[h] = 0;
    cangle[h] = 0;

    if (hthetac[h]->Integral() > 20) {
      gROOT->SetBatch(1);
      int nfound = fSpect->Search(hthetac[h], 1, "", 0.9); // 0.6
      if (nfound > 0) cangle[h] = fSpect->GetPositionX()[0];
      else cangle[h] = hthetac[h]->GetXaxis()->GetBinCenter(hthetac[h]->GetMaximumBin());

      fFit->SetParameters(100, cangle[h], 0.005, 10);
      fFit->FixParameter(2, 0.005);                   // width
      // fFit->SetParLimits(2, 0.001, 0.02);      
      hthetac[h]->Fit("fgaus", "Q", "", cangle[h] - 3.5 * fSigma[h], cangle[h] + 3.5 * fSigma[h]);
      fFit->ReleaseParameter(2); // width
      hthetac[h]->Fit("fgaus", "MQ", "", cangle[h] - 3.5 * fSigma[h], cangle[h] + 3.5 * fSigma[h]);
      cangle[h] = fFit->GetParameter(1);
      spr[h] = fFit->GetParameter(2) * 1000;
      if (fVerbose > 2) gROOT->SetBatch(0);
    }
  }
}

void glob_circleFcn(int &, double *, double &f, double *par, int) {
  int np = glob_gr.GetN();
  f = 0;
  double *x = glob_gr.GetX();
  double *y = glob_gr.GetY();
  for (int i = 0; i < np; i++) {
    double u = x[i] + par[0];
    double v = y[i] + par[1];
    double dr = par[2] - TMath::Sqrt(u * u + v * v);
    f += dr * dr;
  }
}

void glob_circleFcn2(int &, double *, double &f, double *par, int) {
  int np = glob_gr.GetN();
  f = 0;
  double *x = glob_gr.GetX();
  double *y = glob_gr.GetY();
  for (int i = 0; i < np; i++) {
    double u = x[i] + par[0];
    double v = y[i] + par[1];
    double dr = par[2] - TMath::Sqrt(u * u + v * v);
     if (dr < 0.01) f += dr * dr;
     else f += fabs(dr);
  }
}

void PrtReco::FitRing(double &x0, double &y0, double &theta) {

  // TGraph ff_gr;
  // int ff_i(0);
  // int np = glob_gr.GetN();
  // double *x = glob_gr.GetX();
  // double *y = glob_gr.GetY();
  // for (int i = 0; i < np; i++) {
  //   if (fabs(theta - TMath::Sqrt(x[i] * x[i] + y[i] * y[i])) < 0.01) {
  //     ff_gr.SetPoint(ff_i++, x[i], y[i]);
  //   }
  // }
  // glob_gr = ff_gr;
  
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

  // Fit a circle to the graph points
  TVirtualFitter::SetDefaultFitter("Minuit"); // default is Minuit
  TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
  fitter->SetPrecision(0.00000001);
  fitter->SetMaxIterations(1000);
  double arglist[] = {-1};
  fitter->ExecuteCommand("SET PRINT", arglist, 1);
 
  fitter->SetFCN(glob_circleFcn);
  fitter->SetParameter(0, "x0", 0, 0.0005, -0.015, 0.015);
  fitter->SetParameter(1, "y0", 0, 0.0005, -0.015, 0.015);
  fitter->SetParameter(2, "R", theta, 0.001, theta - 0.015, theta + 0.015);

  // fitter->FixParameter(0);
  // fitter->FixParameter(1);
  // fitter->FixParameter(2);

  fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  // fitter->SetFCN(glob_circleFcn2);
  // fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  x0 = fitter->GetParameter(0);
  y0 = fitter->GetParameter(1);
  theta = fitter->GetParameter(2);
}

int PrtReco::FindPdg(double mom, double cangle) {
  double tdiff, diff = 100;
  int minid = 0;
  for (int i = 0; i < 5; i++) {
    tdiff = fabs(cangle - acos(sqrt(mom * mom + ft.mass(i) * ft.mass(i)) / mom /
                               1.46907)); // 1.46907 - fused silica
    if (tdiff < diff) {
      diff = tdiff;
      minid = i;
    }
  }
  return ft.pdg(minid);
}

double PrtReco::FindStartTime(PrtEvent *evt) {

  TVector3 dir, dird, cdir = evt->getMomentum().Unit();
  double tangle, bartime, luttheta, htime, ctime, evtime, dirz;
  bool reflected;
  double shift = 0;
  double speed = 0.1985;
  double lenz = 2100 - evt->getPosition().Z();

  for (PrtHit hit : evt->getHits()) {

    int mcp = hit.getPmt();
    int pix = hit.getPixel();
    int ch = ft.map_pmtpix[mcp][pix];
    htime = hit.getLeadTime() + shift;
    dirz = hit.getMomentum().Z();

    if (dirz < 0) reflected = kTRUE;
    else reflected = kFALSE;    
    if (reflected) lenz = 2 * fRadiatorL - lenz;

    PrtLutNode *node = (PrtLutNode *)fLut->At(ch);
    int size = node->Entries();

    for (int i = 0; i < size; i++) {
      // if(fabs(path-node->GetPathId(i))>0.001) continue;
      dird = node->GetEntry(i);
      evtime = node->GetTime(i);

      for (int u = 0; u < 4; u++) {
        if (u == 0) dir = dird;
        if (u == 1) dir.SetXYZ(dird.X(), -dird.Y(), dird.Z());
        if (u == 2) dir.SetXYZ(-dird.X(), dird.Y(), dird.Z());
        if (u == 3) dir.SetXYZ(-dird.X(), -dird.Y(), dird.Z());
        if (reflected) dir.SetXYZ(dir.X(), dir.Y(), -dir.Z());
        if (dir.Angle(fnX1) < fCriticalAngle || dir.Angle(fnY1) < fCriticalAngle) continue;

        luttheta = dir.Theta();
        if (luttheta > TMath::PiOver2()) luttheta = TMath::Pi() - luttheta;
        bartime = lenz / cos(luttheta) / (1000 * speed);

        ctime = fabs((bartime + evtime));

        tangle = cdir.Angle(dir);

        if (tangle > fAngle[3] - 0.02 && tangle < fAngle[2] + 0.02) {
          fFindTime->Fill(htime - ctime);
          fdtt->Fill(htime - ctime, luttheta * TMath::RadToDeg());
          fdtl->Fill(htime - ctime, speed * htime);
          fdtp->Fill(luttheta * TMath::RadToDeg(), speed * htime);
          int bin = 0.2 * htime;
          if (bin < 20) fFindTimeA[bin]->Fill(htime - ctime);
        }

        if (tangle > 0.35 && tangle < 0.9) {
          hthetac[2]->Fill(tangle);
        }
      }
    }
  }

  if (!fVerbose) gROOT->SetBatch(1);

  // if (fVerbose == 1) ft.add_canvas(Form("hstime_%d", gggg), 800, 400);
  double mean = ft.fit(fFindTime, 3, 20, 2, 1, 0, "QN").X();
  // if(fVerbose==1){
  //   gStyle->SetOptStat(1001111);
  //   fFindTime->Draw();
  //   // tf.waitPrimitive(Form("hstime_%d",gggg));
  //   tf.save_canvas("data/reco",0,1);
  //   gggg++;
  // }

  fFindTimeRes->Fill(mean + shift);
  fFindTime->Reset();
  return mean;
}

void PrtReco::drawTheoryLines(double mom) {
  gPad->Update();

  for (int i = 0; i < 5; i++) {
    fAngle[i] =
      acos(sqrt(mom * mom + ft.mass(i) * ft.mass(i)) / mom / 1.4738); // 1.4738 = 370 = 3.35
  }

  TLine *line = new TLine(0, 0, 0, 1000);
  line->SetX1(fAngle[fp2]);
  line->SetX2(fAngle[fp2]);
  line->SetY1(gPad->GetUymin());
  line->SetY2(gPad->GetUymax());
  line->SetLineColor(ft.color(fp2));
  line->Draw();

  TLine *line1 = new TLine(0, 0, 0, 1000);
  line1->SetX1(fAngle[fp1]);
  line1->SetX2(fAngle[fp1]);
  line1->SetY1(gPad->GetUymin());
  line1->SetY2(gPad->GetUymax());
  line1->SetLineColor(ft.color(fp1));
  line1->Draw();
}

double PrtReco::CalcRejection(TH1F *h1, TH1F *h2, double eff) {
  auto ax1 = h1->GetXaxis();
  auto ax2 = h2->GetXaxis();
  double range = 450;
  double id_total = h1->Integral(ax1->FindBin(-range), ax1->FindBin(range));
  double b;
  for (b = -range; b < range; b += 0.5) {
    double id = h1->Integral(ax1->FindBin(-range), ax1->FindBin(b));
    if (id / id_total > 1 - eff) break;
  }

  double id = h1->Integral(ax1->FindBin(b), ax1->FindBin(range));
  double misid = h2->Integral(ax2->FindBin(b), ax2->FindBin(range));
  if (misid == 0) misid = 0.001;
  return id / misid;
}

void PrtReco::geom_reco(PrtEvent *event, TVector3 mom, bool ringfit) {

  double sum1(0), sum2(0), noise(0.2);
  double evtime, bartime, len, lenz, tdiff, tangle, luttheta, luttime;
  int nph(0), pid = event->getPid();

  for (auto hit : event->getHits()) {

    double hittime = hit.getLeadTime() + gRandom->Gaus(0, fTimeRes);
    double hittime_nc = hittime;
    double dirz = hit.getMomentum().Z();

    int mcp = hit.getPmt();
    int pix = hit.getPixel();
    int ch = hit.getChannel();

    TVector3 dir, dird, dir0 = hit.getMomentum().Unit();    
    double lenz = 0.5 * fRadiatorL - event->getPosition().Z();
    bool reflected = false;
    
    if (dirz < 0) reflected = true;
    else reflected = false;

    if (fabs(dirz) < 1E-6) { // dark noise hit
      if (fTheta > 99) reflected = false;
      else if (fTheta < 81) reflected = true;
      else {
        if (hittime < 42) reflected = false;
        else reflected = true;
      }
    }

    if (reflected) {
      lenz = 2 * fRadiatorL - lenz;
      hittime += cor[mcp].tr;
    } else {
      hittime += cor[mcp].td;
    }

    PrtLutNode *node;
    int size = 0;
    if (fGeomReco) {
      node = (PrtLutNode *)fLut->At(ch);
      size = node->Entries();
    }

    bool isGoodHit_gr(false);

    // double fAngle =  event->GetAngle()-90;
    // TVector3 mom = momInBar;
    // mom.RotateY(-fAngle/180.*TMath::Pi());
    // std::cout<<"fAngle   "<<fAngle <<std::endl;
    // mom.Print();

    Long_t hpath = hit.getPathInPrizm();
    TString spath = Form("%ld", hpath);
    // if(spath.Length()>8) continue;
    // if(!spath.EqualTo("87")) continue;
    // if(spath.Contains("1")) continue;

    for (int i = 0; i < size; i++) {
      dird = node->GetEntry(i);
      evtime = node->GetTime(i);

      Long_t lpath = node->GetPathId(i);
      // TString slpath = Form("%ld", lpath);
      bool ipath = (hpath == lpath) ? 1 : 0;
      // if(!slpath.Contains("4")) continue;
      // if(!ipath) continue;
      // if(lpath!=387) continue;
      // if(node->GetNRefl(i)>8) continue;

      for (int u = 0; u < 4; u++) {
        if (u == 0) dir = dird;
        if (u == 1) dir.SetXYZ(dird.X(), -dird.Y(), dird.Z());
        if (u == 2) dir.SetXYZ(-dird.X(), dird.Y(), dird.Z());
        if (u == 3) dir.SetXYZ(-dird.X(), -dird.Y(), dird.Z());
        if (reflected) dir.SetXYZ(dir.X(), dir.Y(), -dir.Z());
        if (dir.Angle(fnX1) < fCriticalAngle || dir.Angle(fnY1) < fCriticalAngle) continue;

        luttheta = dir.Theta();
        if (luttheta > TMath::PiOver2()) luttheta = TMath::Pi() - luttheta;

        len = lenz / cos(luttheta);
        bartime = len / 199.5; // 198.5
	luttime = bartime + evtime;
        tdiff = hittime - luttime;
	tangle = mom.Angle(dir);
	
	if(!ringfit){
	  fTimeProp->Fill(hittime);
	  fTimeDiffR[reflected]->Fill(tdiff);
	  if (ipath) fTimeDiffR[2]->Fill(tdiff);
	}
	
        tangle += cor[mcp].ca; // per-PMT angle correction;

        if (fabs(tdiff) < 2 && fPhysList < 10)
          tangle -= fChromCor->Eval(fTheta) * (hittime_nc - luttime) / len; // chromatic correction

        if (fabs(tdiff) > fTimeCut + luttime * 0.035) continue;

        if (ringfit) {
	  double chringcut = 0.005;
	  if(frun->getTrackingResTheta() > 0.002) chringcut = 0.01;
          if (fabs(tangle - fAngle[fp2]) > chringcut && fabs(tangle - fAngle[fp1]) > chringcut) continue;

          TVector3 rdir = TVector3(-dir.X(), dir.Y(), dir.Z());
          // rdir.RotateX(0.004);
          // rdir.Rotate(1.57 + TMath::PiOver2(), mom);

          rdir.RotateUz(mom);
          double cphi = rdir.Phi();

          // if(tangle*TMath::Cos(cphi)<0) continue;
          fChRing->Fill(tangle * TMath::Sin(cphi), tangle * TMath::Cos(cphi));
          // fChRing->Fill(tangle * TMath::Sin(cphi), -tangle * TMath::Cos(cphi));
          glob_gr.SetPoint(glob_i, tangle * TMath::Sin(cphi), tangle * TMath::Cos(cphi));
          glob_i++;

        } else {
          fTimeDiff->Fill(hittime, tdiff);
          hthetac[pid]->Fill(tangle);
          hthetacd[pid]->Fill((tangle - fAngle[pid]) * 1000);
          fhChromL->Fill(tdiff / len, (tangle - fAngle[pid]));
          fPmt_a[mcp]->Fill(tangle - fAngle[pid]);
          if (reflected) fPmt_tr[mcp]->Fill(tdiff);
          else fPmt_td[mcp]->Fill(tdiff);

          if (fabs(tangle - fAngle[fp2]) > 0.05 && fabs(tangle - fAngle[fp1]) > 0.05) continue;

          isGoodHit_gr = true;

          sum1 += -TMath::Log(fFunc[fp1]->Eval(tangle) + noise);
          sum2 += -TMath::Log(fFunc[fp2]->Eval(tangle) + noise);
        }
      }
    }

    if (isGoodHit_gr || (!fGeomReco)) {
      nph++;
      if (frun->getPid() == 10005) {
        if (pid == 3) ft.fill_digi(mcp, pix);
      } else if (pid == 2) ft.fill_digi(mcp, pix);
    }
  }

  if (!ringfit) {
    double sum_gr = sum1 - sum2;
    if (sum_gr != 0) fLnDiffGr[pid]->Fill(sum_gr);
    if (nph > 1) hnph_gr[pid]->Fill(nph);
  }
}

void PrtReco::time_imaging(PrtEvent *event) {

  int pid = event->getPid();
  double sum1(0), sum2(0), noise(0.5e-5);
  int nph(0);
 
  for (auto hit : event->getHits()) {
    
    double t = hit.getLeadTime() + gRandom->Gaus(0, fTimeRes);
    int mcp = hit.getPmt();
    int pix = hit.getPixel();
    int ch = hit.getChannel();

    if (fMethod == 2) {
      nph++;

      double lh1 = fTime[fp1][ch]->GetBinContent(fTime[fp1][ch]->FindBin(t));
      double lh2 = fTime[fp2][ch]->GetBinContent(fTime[fp2][ch]->FindBin(t));
      // double lh1 = fPdf[fp1][ch]->Eval(t);
      // double lh2 = fPdf[fp2][ch]->Eval(t);

      if (lh1 < 0) lh1 = 0;
      if (lh2 < 0) lh2 = 0;

      sum1 += TMath::Log((lh1 + noise));
      sum2 += TMath::Log((lh2 + noise));

      if (0) {
        TString x = (sum1 > sum2) ? " <====== PROTON" : "";
        std::cout << Form("f %1.6f s %1.6f mcp %d pix %d   pid %d", sum1, sum2, mcp, pix, pid)
                  << "  " << x << std::endl;

        ft.add_canvas("ctemp", 800, 400);
        // ft.normalize(fPdf4[ch],fPdf2[ch]);
        fPdf[fp2][ch]->SetLineColor(2);
        fPdf[fp2][ch]->SetLineColor(4);
        fPdf[fp2][ch]->Draw("APL");
        fPdf[fp1][ch]->SetTitle(Form("mcp=%d  pix=%d", mcp, pix));
        fPdf[fp1][ch]->GetXaxis()->SetTitle("LE time [ns]");
        fPdf[fp1][ch]->GetYaxis()->SetTitle("PDF value");
        fPdf[fp1][ch]->GetXaxis()->SetRangeUser(0, 40);
        fPdf[fp1][ch]->Draw("PL same");
        gPad->Update();
        TLine *gLine = new TLine(0, 0, 0, 1000);
        gLine->SetLineWidth(2);
        gLine->SetX1(t);
        gLine->SetX2(t);
        gLine->SetY1(gPad->GetUymin());
        gLine->SetY2(gPad->GetUymax());
        gLine->Draw();
        gPad->Update();
        gPad->WaitPrimitive();
      }
    }

    if (fMethod == 4) {
      double w = 1;
      fTotal[pid]++;
      fTime[pid][ch]->Fill(t, w);
    }
  }

  double sum_nph = 0;
  
  if (nph > 1) hnph_ti[pid]->Fill(nph);
  double sum_ti = 1.5 * (sum1 - sum2) + 30 * sum_nph;
  if (fabs(sum_ti) > 0.1) fLnDiffTi[pid]->Fill(1.0 * sum_ti);
}

void PrtReco::nn_reco(PrtEvent *event, TVector3 mom) {
#ifdef AI
  std::vector<double> vinput(6144, 0.0);
  std::vector<int> vinput2d(6144 * 350, 0);
  std::vector<int> vinputInd(200 * 3, 0);

  int pid = event->getPid();
 
  for (auto hit : event->getHits()) {

    double hittime = hit.getLeadTime() + gRandom->Gaus(0, fTimeRes);
    int mcp = hit.getPmt();
    int pix = hit.getPixel();
    int ch = hit.getChannel();

    int tpix = pix - 1;
    int tx = int(16 * (mcp % 6) + tpix % 16);
    int ty = int(16 * (mcp / 6) + tpix / 16);
    int tc = 16 * 6 * ty + tx;
    if (tc > -1) vinput[tc] = hittime;

    if (hittime < 70) {
      int ic = ch * 350 + int(5 * hittime);
      vinput2d[ic] = 1;
      if (fNph_gr[pid] < 200) {
        vinputInd[fNph_gr[pid] * 3 + 0] = 0;
        vinputInd[fNph_gr[pid] * 3 + 1] = ch;
        vinputInd[fNph_gr[pid] * 3 + 2] = int(5 * hittime);
      }
    }
  }

  if (fNNet) { // neural network

    // std::vector<int> input(6144,0);
    // input = cppflow::cast(input, TF_UINT8, TF_FLOAT);

    // Creates a tensor from the vector with shape [X_dim, Y_dim]

    // auto input = cppflow::tensor(vinput2d, {1, 6144, 350});
    auto input = cppflow::tensor(vinputInd, {1, 200, 3});
    input = cppflow::cast(input, TF_INT32, TF_INT32);

    // auto input = cppflow::tensor(vinput, {1, 64, 96, 1});
    // input = cppflow::cast(input, TF_FLOAT, TF_INT64);
    auto output = (*fNNmodel)(input);
    auto t = cppflow::arg_max(output, 1).get_tensor();
    float *ll = static_cast<float *>(TF_TensorData(output.get_tensor().get()));
    int nn_pid = static_cast<int *>(TF_TensorData(t.get()))[0];
    // if (nn_pid == 0) nn_pid = 2;
    // if (nn_pid == 1) nn_pid = 3;
    double dll = 5 * (ll[fp2] - ll[fp1]);

    if (fabs(0.0 - dll) > 0.3) fLnDiffNn[pid]->Fill(dll);
    // Show the predicted class
    // std::cout << output << std::endl;
    // std::cout << "PID " << pid << " nn " << nn_pid << std::endl;
    if (pid == nn_pid) eff_nn[nn_pid]++;
    eff_total[pid]++;
  }
#endif
}
