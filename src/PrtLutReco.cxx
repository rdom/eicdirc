// -----------------------------------------
// PrtLutReco.cpp
//
// Created on: 13.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------

#include "PrtLutReco.h"

using std::cout;
using std::endl;

TH1F *fHist1 = new TH1F("Time1", "1", 1000, 0, 20);
TH1F *fHist2 = new TH1F("Time2", "2", 1000, -10, 10);
TH1F *fHistDiff[3];
TH2F *fDiff =
  new TH2F("diff", ";measured time [ns];t_{measured}-t_{calc} [ns]", 500, 0, 100, 150, -5, 5);
TH2F *fHist4 = new TH2F("Time4", "4", 200, -1, 1, 200, -1, 1);
TH2F *fHist5 = new TH2F("Time5", "5", 200, -1, 1, 200, -1, 1);
TH1F *fFindTime = new TH1F("ft", ";t_{measured}-t_{calculated} [ns];entries [#]", 2000, -100, 100);
TH1F *fFindTimeA[20];
TH1F *fFindTimeRes = new TH1F("ftr", "ftr", 100, -2, 2);
TH2F *fdtt =
  new TH2F("dtt", ";t_{measured}-t_{calculated} [ns];#theta_{l} [deg]", 1000, -2, 2, 1000, 0, 90);
TH2F *fdtl =
  new TH2F("dtl", ";t_{measured}-t_{calculated} [ns];path length [m]", 1000, -2, 2, 1000, 0, 15);
TH2F *fdtp = new TH2F("dtp", ";#theta_{l} [deg];path length [m]", 1000, 0, 90, 1000, 0, 15);
TH2F *fhChrom =
  new TH2F("chrom", ";t_{measured}-t_{calculated} [ns];#theta_{C} [rad]", 100, -2, 2, 100, -30, 30);
TH2F *fhChromL = new TH2F("chroml", ";(t_{measured}-t_{calculated})/L_{path};#theta_{C} [rad]", 100,
                          -0.0002, 0.0002, 100, -30, 30);
TH1F *fHistMcp[28];
double fCorr[28];

int gg_i(0);
TGraph gg_gr;

// -----   Default constructor   -------------------------------------------
PrtLutReco::PrtLutReco(TString infile, TString lutfile, TString pdffile, int verbose) {
  fVerbose = verbose;
  fCriticalAngle = asin(1.00028 / 1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
  fp1 = 3;
  fp2 = 2;

  fChain = new TChain("data");
  fChain->Add(infile);
  fEvent = new PrtEvent();
  frun = PrtManager::Instance()->getRun();
  ft = PrtTools(frun);
  fmaxch = ft.maxdircch();
  fnpmt = frun->getNpmt();
  fnpix = frun->getNpix();
  fMethod = frun->getRunType();
  fStudyId = frun->getStudy();
  fRadiator = frun->getRadiator();
  int fileId = frun->getId();
  fPk = 4;
  if (frun->getPid() == 10002) fPk = 3; //  kaon

  fChain->SetBranchAddress("PrtEvent", &fEvent);

  fFile = new TFile(lutfile);
  fTree = (TTree *)fFile->Get("prtlut");
  fLut = new TClonesArray("PrtLutNode");
  fTree->SetBranchAddress("LUT", &fLut);
  fTree->GetEntry(0);

  fFit = new TF1("fgaus", "[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +[3]", 0.35, 0.9);
  fSpect = new TSpectrum(10);

  int col[] = {kRed + 1, kBlue + 1, kBlack};
  for (int i = 0; i < 3; i++) {
    fHistDiff[i] = new TH1F(Form("TimeDiff_%d", i), ";t_{measured}-t_{calculated} [ns];entries [#]",
                            500, -10, 10);
    fHistDiff[i]->SetLineColor(col[i]);
  }
  for (int h = 0; h < 5; h++) {
    hthetac[h] = new TH1F(Form("thetac_%d", h), ";#theta_{C} [rad];entries [#]", 200, 0.75, 0.9);
    hthetacd[h] =
      new TH1F(Form("thetacd_%d", h), ";#Delta#theta_{C} [mrad];entries [#]", 200, -60, 60);
    hnph[h] = new TH1F(Form("nph_%d", h), ";detected photons [#];entries [#]", 150, 0, 150);
    fFunc[h] = new TF1(Form("gaus_%d", h), "[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])", 0.7, 0.9);

    int c = ft.color(h);
    hthetac[h]->SetLineColor(c);
    hthetacd[h]->SetLineColor(c);
    hnph[h]->SetLineColor(c);
    fFunc[h]->SetLineColor(c);
  }

  for (int i = 0; i < 20; i++) {
    fFindTimeA[i] =
      new TH1F(Form("fta_%d", i), ";t_{measured}-t_{calculated} [ns];entries [#]", 1000, -10, 10);
  }

  for (int i = 0; i < 28; i++) {
    fHistMcp[i] = new TH1F(Form("fHistMcp_%d", i),
                           Form("fHistMcp_%d;#theta_{C} [rad];entries [#]", i), 50, -0.05, 0.05);
  }

  // read corrections
  fCorrFile = infile + "_corr.root";
  for (int i = 0; i < fnpmt; i++) fCorr[i] = 0;
  if (!gSystem->AccessPathName(fCorrFile)) {
    std::cout << "------- reading  " << fCorrFile << std::endl;
    int pmt;
    double corr, cspr[5];
    TChain ch;
    ch.SetName("corr");
    ch.Add(fCorrFile);
    ch.SetBranchAddress("pmt", &pmt);
    ch.SetBranchAddress("corr", &corr);
    ch.SetBranchAddress("cspr", &cspr);
    for (int i = 0; i < ch.GetEntries(); i++) {
      ch.GetEvent(i);
      fCorr[pmt] = (fabs(corr) < 0.007) ? corr : 0.00001;
      fSigma = 0.001 * cspr[2] * 0.9;
      std::cout << "pmt " << pmt << "  " << fCorr[pmt] << " spr = " << fSigma << std::endl;
    }
  } else {
    std::cout << "------- corr file not found  " << fCorrFile << std::endl;
  }

  cout << "-I- PrtLutReco: Intialization successfull" << endl;
}

// -----   Destructor   ----------------------------------------------------
PrtLutReco::~PrtLutReco() {}

//-------------- Loop over tracks ------------------------------------------
void PrtLutReco::Run(int start, int end) {
  TVector3 dird, dir, momInBar(0, 0, 1), posInBar;
  double mom, tangle, tdiff, evtime, bartime, lenz, dirz, luttheta, hitTime;
  int tofPid(0), distPid(0), likePid(0);
  bool reflected = kFALSE;
  gStyle->SetOptFit(111);

  TVector3 fnX1 = TVector3(1, 0, 0);
  TVector3 fnY1 = TVector3(0, 1, 0);
  int nsHits(0), nsEvents(0);

  TString outFile = PrtManager::Instance()->getOutName();
  double theta(0), phi(0), cangle[5] = {0}, spr[5] = {0}, trr[5] = {0}, nph[5] = {0},
                           nph_err[5] = {0}, par5(0), par6(0), timeRes(0), ctimeRes(0), test1(0),
                           test2(0), test3(0), sep(0), sep_err(0);

   
  ft.set_palette(1);
  ft.create_maps();  
  ft.init_digi();  

  outFile.ReplaceAll("reco_", Form("reco_%d_", frun->getId()));
  TFile file(outFile, "recreate");
  TTree tree("reco", "reco");
  tree.Branch("mom", &mom, "mom/D");
  tree.Branch("tofPid", &tofPid, "tofPid/I");
  tree.Branch("distPid", &distPid, "distPid/I");
  tree.Branch("likePid", &likePid, "likePid/I");
  tree.Branch("spr", &spr, "spr[5]/D");
  tree.Branch("trr", &trr, "trr[5]/D");
  tree.Branch("nph", &nph, "nph[5]/D");
  tree.Branch("nph_err", &nph_err, "nph_err[5]/D");
  tree.Branch("cangle", &cangle, "cangle[5]/D");
  tree.Branch("sep", &sep, "sep/D");
  tree.Branch("sep_err", &sep_err, "sep_err/D");
  tree.Branch("par5", &par5, "par5/D");
  tree.Branch("par6", &par6, "par6/D");
  tree.Branch("timeres", &timeRes, "timeRes/D");
  tree.Branch("ctimeres", &ctimeRes, "ctimeRes/D");
  tree.Branch("test1", &test1, "test1/D");
  tree.Branch("test2", &test2, "test2/D");
  tree.Branch("test3", &test3, "test3/D");
  tree.Branch("theta", &theta, "theta/D");
  tree.Branch("phi", &phi, "phi/D");
  
  test1 = frun->getTest1();
  test2 = frun->getTest2();
  test3 = frun->getTest3();
  fMethod = frun->getRunType();
  timeRes = frun->getTimeSigma();
  double timeCut = 0.5;
  int pid = frun->getPid();
  fp1 = 2; // pi
  if(pid == 10001)  fp1 = 0; // e
  if(pid == 10002)  fp1 = 3; // K
  if(pid == 10003)  fp1 = 1; // mu
  
  if( fCor_level == 0) timeRes = 1.5;

  int nEvents = fChain->GetEntries();
  if (end == 0) end = nEvents;

  std::cout << "Run started for [" << start << "," << end << "]" << std::endl;

  int fEvId = 2030; // 2032

  for (int ievent = 0; ievent < nEvents; ievent++) {
    fChain->GetEntry(ievent);
    theta = (fEvent->getMomentum().Angle(TVector3(0, 0, -1))) * TMath::RadToDeg();
    mom = fEvent->getMomentum().Mag() / 1000.;
    pid = fEvent->getPid();
    int tnph[5] = {0};

    if (ievent % 1000 == 0)
      std::cout << "Event # " << ievent << " has " << fEvent->getHits().size() << " hits" << std::endl;
    double minChangle = 0.35;
    double maxChangle = 0.9;
    TVector3 rotatedmom = fEvent->getMomentum().Unit();
    double sum1(0), sum2(0), noise(0.2);

    // // track smearing
    // TVector3 init = rotatedmom;
    // rotatedmom.RotateY(gRandom->Gaus(0, test1));
    // rotatedmom.Rotate(TMath::Pi(), init);

    if (fSigma < 0.003) fSigma = 0.007;

    if (ievent == 0) {
      double range = 160;
      if (mom > 1.5) range = 100;
      if (mom > 2) range = 60;
      if (mom < 1) range = 300;
      if (mom < 0.7) range = 400;
      if (mom < 0.6) range = 500;

      if (fp2 == 3 && theta < 120) range = 50;
      if (fp2 == 3 && mom < 4) range = 500;

      for (int h = 0; h < 5; h++) {
        fLnDiff[h] =
          new TH1F(Form("LnDiff_%d", h), ";ln L(K) - ln L(#pi);entries [#]", 100, -range, range);
        fLnDiff[h]->SetLineColor(ft.color(h));
      }
    }

    // if(fp1==0 || fp2==0){ //electron
    //   fSigma=0.0057;
    //   fSigma=0.0045;
    //   if(theta<40 ) fSigma=0.007;
    //   if(mom<0.7) fSigma=0.008;
    //   if(mom>2) fSigma=0.006;
    // }

    // if(fp2==3){ // kaon
    //   if(theta<120) fSigma=0.004;
    // }

    for (int i = 0; i < 5; i++) {
      fAngle[i] =
        acos(sqrt(mom * mom + ft.mass(i) * ft.mass(i)) / mom / 1.4738); // 1.4738 = 370 = 3.35
      fFunc[i]->SetParameter(0, 1);
      fFunc[i]->SetParameter(1, fAngle[i]);
      fFunc[i]->SetParameter(2, fSigma);
      // if(i==0) fFunc[i]->SetParameter(2,0.0074);
      // if(i==2) fFunc[i]->SetParameter(2,0.0064);
    }

    // double stime = FindStartTime(fEvent);

   
    
    for (auto hit : fEvent->getHits()) {
   
      hitTime = hit.getLeadTime() + gRandom->Gaus(0, 0.1);
      dirz = hit.getMomentum().Z();
      int mcp = hit.getPmt();
      int pix = hit.getPixel();
      int ch = hit.getChannel();

      TVector3 dir0 = hit.getMomentum().Unit();
      TVector3 cz = TVector3(-rotatedmom.X(), rotatedmom.Y(), rotatedmom.Z());
      TVector3 cd = TVector3(-dir0.X(), dir0.Y(), dir0.Z());
      TVector3 unitdir1 = rotatedmom.Unit();
      TVector3 unitdir2 = rotatedmom.Unit();
      cz.RotateUz(unitdir1);
      cd.RotateUz(unitdir2);
      
      lenz = 2100 - fEvent->getPosition().Z();      
      double phi0 = cd.Phi();
      if (dirz < 0) {	
        reflected = true;
        lenz = 2 * 4200 - lenz;
      } else {
        reflected = false;
      }
      
      double theta0 = rotatedmom.Angle(dir0);
      fHist5->Fill(theta0 * TMath::Sin(phi0), theta0 * TMath::Cos(phi0));
      
      PrtLutNode *node = (PrtLutNode *)fLut->At(ch);
      int size = node->Entries();
      bool isGoodHit(false);
      
      // double fAngle =  fEvent->GetAngle()-90;
      // TVector3 rotatedmom = momInBar;
      // rotatedmom.RotateY(-fAngle/180.*TMath::Pi());
      // std::cout<<"fAngle   "<<fAngle <<std::endl;
      // rotatedmom.Print();

      Long_t hpath = hit.getPathInPrizm();
      TString spath = Form("%ld", hpath);
      // if(spath.Length()>8) continue;

      // std::cout<<ch<<" "<<mcp<<" ========================= spath "<<spath<<std::endl;

      // if(!spath.EqualTo("87")) continue;
      // if(spath.Contains("1")) continue;

      for (int i = 0; i < size; i++) {
        dird = node->GetEntry(i);

        Long_t lpath = node->GetPathId(i);
        TString slpath = Form("%ld", lpath);
        bool ipath = 0;
        if (hpath == lpath) ipath = 1;
        // if(!slpath.Contains("4")) continue;
        // std::cout<<"slpath "<<slpath<<std::endl;
        // if(!ipath) continue;
        // if(lpath!=387) continue;
        // if(node->GetNRefl(i)>8) continue;

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

          bartime = lenz / cos(luttheta) / 198.5; // 198

          fHist1->Fill(hitTime);
          double luttime = bartime + evtime;
	  
          tdiff = hitTime - luttime;
          fHistDiff[reflected]->Fill(tdiff);
          if (ipath) fHistDiff[2]->Fill(tdiff);

          tangle = rotatedmom.Angle(dir) + fCorr[mcp]; // 45;
          // if(tangle>TMath::PiOver2()) tangle = TMath::Pi()-tangle;

          if (fabs(tdiff) < 2) tangle -= 0.012 * tdiff; // chromatic correction
	  if (fabs(tdiff) > timeCut + luttime * 0.035) continue;
          fDiff->Fill(hitTime, tdiff);

          // if(theta<50){
          //   if(fabs(tdiff)<1.5) tangle -= 0.005*tdiff; // chromatic correction
          // }
          // if(fabs(theta-90)<10){
          //   if(reflected && fabs(tdiff)<1) tangle -= 0.0128*tdiff;
          // }
          // if(fabs(theta-145)<10) if(fabs(tdiff)<1.5) tangle += 0.02*tdiff;

          // if(tangle>TMath::Pi()/2.) tangle = TMath::Pi()-tangle;
          // if(fabs(0.8218-tangle)>0.002) continue;
          // if(fabs(0.83-tangle)>0.003) continue;

          hthetac[pid]->Fill(tangle);
          hthetacd[pid]->Fill((tangle - fAngle[pid]) * 1000);
          fHistMcp[mcp]->Fill(tangle - fAngle[pid]);
          fhChrom->Fill(tdiff, (tangle - fAngle[pid]) * 1000);
          fhChromL->Fill(tdiff / (lenz / cos(luttheta)), (tangle - fAngle[pid]) * 1000);

          if (fabs(tangle - fAngle[fp2]) > 0.03 && fabs(tangle - fAngle[fp1]) > 0.03) continue;

          if (tangle > minChangle && tangle < maxChangle) {
            TVector3 rdir = TVector3(-dir.X(), dir.Y(), dir.Z());
            TVector3 unitdir3 = rotatedmom.Unit();
            rdir.RotateUz(unitdir3);
            double cphi = rdir.Phi();
            // if(tangle*TMath::Cos(cphi)<0) continue;
            fHist4->Fill(tangle * TMath::Sin(cphi), tangle * TMath::Cos(cphi));
            gg_gr.SetPoint(gg_i, tangle * TMath::Sin(cphi), tangle * TMath::Cos(cphi));
            gg_i++;
          }

          isGoodHit = true;

          sum1 += -TMath::Log(fFunc[fp1]->Eval(tangle) + noise);
          sum2 += -TMath::Log(fFunc[fp2]->Eval(tangle) + noise);
        }
      }

      if (isGoodHit) {
        nsHits++;
        tnph[pid]++;
	if (pid == 2) ft.fill_digi(mcp,pix);
      }
    }

    double sum = sum1 - sum2;

    if (sum != 0) fLnDiff[pid]->Fill(sum);
    if (tnph[pid] > 1) hnph[pid]->Fill(tnph[pid]);

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

  if (fMethod == 2) {
    FindPeak(cangle, spr);

    for (int h = 0; h < 5; h++) {
      if (hnph[h]->GetEntries() < 20) continue;
      hnph[h]->Fit("gaus", "", "S", 5, 150);
      auto f = hnph[h]->GetFunction("gaus");
      if (f) {
        nph[h] = f->GetParameter(1);
        nph_err[h] = f->GetParError(1);
      }
    }
    
    TF1 *ff;
    double m1 = 0, m2 = 0, s1 = 0, s2 = 0, dm1 = 0, dm2 = 0, ds1 = 0, ds2 = 0;
    if (fLnDiff[fp2]->GetEntries() > 10) {
      fLnDiff[fp2]->Fit("gaus", "S");
      ff = fLnDiff[fp2]->GetFunction("gaus");
      m1 = ff->GetParameter(1);
      s1 = ff->GetParameter(2);
      dm1 = ff->GetParError(1);
      ds1 = ff->GetParError(2);

      if (fp1 == 0 && mom < 1.5) { // handle tails
        fLnDiff[fp2]->Fit("gaus", "S", "", m1 - 1.8 * s1, 500);
        ff = fLnDiff[fp2]->GetFunction("gaus");
        m1 = ff->GetParameter(1);
        s1 = ff->GetParameter(2);
        dm1 = ff->GetParError(1);
        ds1 = ff->GetParError(2);
      }
    }
    if (fLnDiff[fp1]->GetEntries() > 10) {
      fLnDiff[fp1]->Fit("gaus", "S");
      ff = fLnDiff[fp1]->GetFunction("gaus");
      m2 = ff->GetParameter(1);
      s2 = ff->GetParameter(2);
      dm2 = ff->GetParError(1);
      ds2 = ff->GetParError(2);

      if (fp1 == 0 && mom < 1.5) { /// handle tails
        fLnDiff[fp1]->Fit("gaus", "S", "", -500, m2 + 1.8 * s2);
        ff = fLnDiff[fp1]->GetFunction("gaus");
        m2 = ff->GetParameter(1);
        s2 = ff->GetParameter(2);
        dm2 = ff->GetParError(1);
        ds2 = ff->GetParError(2);
      }
    }
    sep = (fabs(m1 - m2)) / (0.5 * (s1 + s2));

    double e1, e2, e3, e4;
    e1 = 2 / (s1 + s2) * dm1;
    e2 = -2 / (s1 + s2) * dm2;
    e3 = -2 * (m1 - m2) / ((s1 + s2) * (s1 + s2)) * ds1;
    e4 = -2 * (m1 - m2) / ((s1 + s2) * (s1 + s2)) * ds2;
    sep_err = sqrt(e1 * e1 + e2 * e2 + e3 * e3 + e4 * e4);

    std::cout << Form("%3d : SPR = %2.2f N = %2.2f +/- %2.2f", ft.pdg(fp1), spr[fp1], nph[fp1],
                      nph_err[fp1])
              << std::endl;
    std::cout << Form("%3d : SPR = %2.2f N = %2.2f +/- %2.2f", ft.pdg(fp2), spr[fp2], nph[fp2],
                      nph_err[fp2])
              << std::endl;
    std::cout << Form("SEP = %2.2f +/- %2.2f ", sep, sep_err) << std::endl;
  }

  if (!fVerbose) gROOT->SetBatch(1);

  if (0) { // draw start time
    ft.add_canvas(Form("ctimeres_%d", int(theta + 0.01)), 800, 400);
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

  tree.Fill();
  tree.Write();

  if (fVerbose > 1) {
    TString nid = Form("_%1.2f_%1.4f_%1.2f", theta, test1, mom);

    { // cherenkov angle
      ft.add_canvas("tangle" + nid, 800, 400);
      ft.normalize(hthetac, 5);

      hthetac[fp1]->SetTitle(Form("theta %1.2f", theta));
      hthetac[fp1]->Draw("");
      hthetac[fp2]->Draw("same");
      drawTheoryLines(6);

      ft.add_canvas("tangled" + nid, 800, 400);
      ft.normalize(hthetacd, 5);
      hthetacd[fp2]->SetTitle(Form("theta %1.2f", theta));
      hthetacd[fp2]->Draw("");
      hthetacd[fp1]->Draw("same");
    }

    { // nph
      ft.add_canvas("nph" + nid, 800, 400);
      ft.normalize(hnph, 5);
      hnph[fp1]->SetStats(0);
      hnph[fp1]->Draw();
      hnph[fp2]->Draw("same");
    }

    { // sep
      fLnDiff[fp1]->SetStats(0);
      fLnDiff[fp2]->SetStats(0);

      ft.add_canvas("lh" + nid, 800, 400);
      ft.normalize(fLnDiff, 5);
      fLnDiff[fp2]->SetName(Form("s_%2.2f", sep));
      fLnDiff[fp2]->SetTitle(Form("separation = %2.2f s.d.", sep));
      fLnDiff[fp2]->GetXaxis()->SetTitle("ln L(" + ft.lname(fp2) + ") - ln L(" + ft.lname(fp1) +
                                         ")");
      fLnDiff[fp2]->Draw();
      fLnDiff[fp1]->Draw("same");
    }

    { // chromatic corrections
      ft.add_canvas("chrom" + nid, 800, 400);
      fhChrom->SetStats(0);
      fhChrom->Draw("colz");
      ft.add_canvas("chroml" + nid, 800, 400);
      fhChromL->SetStats(0);
      fhChromL->Draw("colz");
    }

    { // hp
      auto cdigi = ft.draw_digi(0, 0);
      cdigi->SetName("hp" + nid);
      ft.add_canvas(cdigi);
    }

    { // cherenkov ring
      ft.add_canvas("ring" + nid, 500, 500);

      fHist4->SetStats(0);
      fHist4->GetXaxis()->SetTitle("#theta_{c}sin(#varphi_{c})");
      fHist4->GetYaxis()->SetTitle("#theta_{c}cos(#varphi_{c})");
      fHist4->SetTitle(Form("Calculated from LUT, #theta = %1.2f#circ", theta));
      fHist4->Draw("colz");
      double x0(0), y0(0);
      FitRing(x0, y0, fAngle[2]);
      TVector3 corr(x0, y0, 1 - TMath::Sqrt(x0 * x0 + y0 * y0));
      // std::cout<<"Tcorr "<< corr.Theta()*1000<< "  Pcorr "<< corr.Phi() <<std::endl;

      TLegend *leg = new TLegend(0.32, 0.42, 0.67, 0.59);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry((TObject *)0, Form("Entries %0.0f", fHist4->GetEntries()), "");
      leg->AddEntry((TObject *)0, Form("#Delta#theta_{c} %f [mrad]", corr.Theta() * 1000), "");
      leg->AddEntry((TObject *)0, Form("#Delta#varphi_{c} %f [rad]", corr.Phi()), "");
      leg->Draw();

      TArc *arc = new TArc(x0, y0, fAngle[2]);
      arc->SetLineColor(kRed);
      arc->SetLineWidth(1);
      arc->SetFillStyle(0);
      arc->Draw();
      gg_i = 0;
      gg_gr.Set(0);
    }

    { // corrections
      if (fabs(fCorr[0]) < 0.00000001 && fabs(fCorr[10]) < 0.00000001 &&
          fabs(fCorr[20]) < 0.00000001) {
        std::cout << "Writing " << fCorrFile << std::endl;

        TFile fc(fCorrFile, "recreate");
        TTree *tc = new TTree("corr", "corr");
        int pmt;
        double corr;
        tc->Branch("pmt", &pmt, "pmt/I");
        tc->Branch("corr", &corr, "corr/D");
        tc->Branch("cspr", &spr, "cspr[5]/D");

        for (int i = 0; i < fnpmt; i++) {
          if (fHistMcp[i]->GetEntries() < 20) continue;
          ft.add_canvas(Form("r_tangle_%d", i), 800, 400);

          corr = -ft.fit(fHistMcp[i], 0.008, 20, 0.01).X();
          // fHistMcp[i]->Fit("gaus","MQ","",-0.01,0.01).X();
          // auto f = fHistMcp[i]->GetFunction("gaus");
          // if(f)
          {
            pmt = i;
            //  corr= -f->GetParameter(1);
            tc->Fill();
            std::cout << "if(mcpid==" << i << ") tangle += " << corr << ";" << std::endl;
            fHistMcp[i]->Draw();
          }
        }

        tc->Write();
        fc.Write();
        fc.Close();
      }
    }

    { // time
      ft.add_canvas("tdiff" + nid, 800, 400);
      ft.normalize(fHistDiff, 3);
      for (int i = 0; i < 3; i++) {
        if (fHistDiff[i]->GetEntries() > 100) {
          fHistDiff[i]->SetStats(0);
          fHistDiff[i]->SetTitle(Form("theta %1.2f", theta));
          fHistDiff[i]->Draw((i == 0) ? "h" : "hsame");
        }
      }

      ft.add_canvas("diff" + nid, 800, 400);
      fDiff->SetStats(0);
      fDiff->Draw("colz");
    }

    TString filedir = fCorrFile;
    if (filedir.Contains("/")) {
      filedir.Remove(filedir.Last('/'));
      filedir += "/";
    } else filedir = "";

    ft.save_canvas(filedir + "reco", 1, 0, 0);

    if (fVerbose > 2) ft.wait_primitive("lh" + nid, "none");
  }
}

void PrtLutReco::FindPeak(double (&cangle)[5], double (&spr)[5]) {
  for (int h = 0; h < 5; h++) {
    spr[h] = 0;
    cangle[h] = 0;

    if (hthetac[h]->GetEntries() > 20) {
      gROOT->SetBatch(1);
      int nfound = fSpect->Search(hthetac[h], 1, "", 0.9); // 0.6
      if (nfound > 0) cangle[h] = fSpect->GetPositionX()[0];
      else cangle[h] = hthetac[h]->GetXaxis()->GetBinCenter(hthetac[h]->GetMaximumBin());

      fFit->SetParameters(100, cangle[h], 0.005, 10); // peak
      fFit->SetParameter(2, 0.005);                   // width
      fFit->FixParameter(2, 0.008);                   // width
      hthetac[h]->Fit("fgaus", "Q", "", cangle[h] - 3.5 * fSigma, cangle[h] + 3.5 * fSigma);
      fFit->ReleaseParameter(2); // width
      hthetac[h]->Fit("fgaus", "MQ", "", cangle[h] - 3.5 * fSigma, cangle[h] + 3.5 * fSigma);
      cangle[h] = fFit->GetParameter(1);
      spr[h] = fFit->GetParameter(2) * 1000;
      if (fVerbose > 2) gROOT->SetBatch(0);
    }
  }
}

void circleFcn(int &, double *, double &f, double *par, int) {
  int np = gg_gr.GetN();
  f = 0;
  double *x = gg_gr.GetX();
  double *y = gg_gr.GetY();
  for (int i = 0; i < np; i++) {
    double u = x[i] + par[0];
    double v = y[i] + par[1];
    double dr = par[2] - TMath::Sqrt(u * u + v * v);
    f += dr * dr;
  }
}

void circleFcn2(int &, double *, double &f, double *par, int) {
  int np = gg_gr.GetN();
  f = 0;
  double *x = gg_gr.GetX();
  double *y = gg_gr.GetY();
  for (int i = 0; i < np; i++) {
    double u = x[i] + par[0];
    double v = y[i] + par[1];
    double dr = par[2] - TMath::Sqrt(u * u + v * v);
    if (dr > 0.07) f += dr * dr;
    else f += fabs(dr);
  }
}

void PrtLutReco::FitRing(double &x0, double &y0, double &theta) {

  TGraph ff_gr;
  int ff_i(0);
  int np = gg_gr.GetN();
  double *x = gg_gr.GetX();
  double *y = gg_gr.GetY();
  for (int i = 0; i < np; i++) {
    if (fabs(theta - TMath::Sqrt(x[i] * x[i] + y[i] * y[i])) < 0.05) {
      ff_gr.SetPoint(ff_i, x[i], y[i]);
      ff_i++;
    }
  }
  gg_gr = ff_gr;

  // Fit a circle to the graph points
  TVirtualFitter::SetDefaultFitter("Minuit"); // default is Minuit
  TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
  fitter->SetPrecision(0.00000001);
  fitter->SetMaxIterations(1000);
  double arglist[] = {-1};
  fitter->ExecuteCommand("SET PRINT", arglist, 1);

  fitter->SetFCN(circleFcn);
  fitter->SetParameter(0, "x0", 0.03, 0.01, -0.05, 0.05);
  fitter->SetParameter(1, "y0", 0, 0.01, -0.05, 0.05);
  fitter->SetParameter(2, "R", theta, 0.01, theta - 0.05, theta + 0.05);

  // fitter->FixParameter(0);
  // fitter->FixParameter(1);

  fitter->FixParameter(2);
  fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  // fitter->SetFCN(circleFcn2);
  // fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  x0 = fitter->GetParameter(0);
  y0 = fitter->GetParameter(1);
  theta = fitter->GetParameter(2);
}

int PrtLutReco::FindPdg(double mom, double cangle) {
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

int gggg = 0;
double PrtLutReco::FindStartTime(PrtEvent *evt) {
  TVector3 fnX1 = TVector3(1, 0, 0);
  TVector3 fnY1 = TVector3(0, 1, 0);
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
    if (reflected) lenz = 2 * 4200 - lenz;

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

  if (fVerbose == 1) ft.add_canvas(Form("hstime_%d", gggg), 800, 400);
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

void PrtLutReco::drawTheoryLines(double mom) {
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
