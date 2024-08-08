// -----------------------------------------
// PrtLutReco.h
//
// Created on: 13.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------
// Class for reconstruction in DIRC using look-up table method

#ifndef PrtLutReco_h
#define PrtLutReco_h 1

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TArc.h"
#include "TVirtualFitter.h"

#include "PrtRun.h"
#include "PrtEvent.h"
#include "PrtHit.h"
#include "PrtManager.h"
#include "PrtLutNode.h"
#include "PrtTools.h"

#ifdef AI
#include "cppflow/cppflow.h"
#endif

class PrtLutReco {

 public:
  // Standard constructors
  PrtLutReco(TString infile, TString lutfile, TString pdffile, TString nnfile, int verbose = 0);

  // Destructor
  ~PrtLutReco();
  void Run(int start = 0, int end = 0);
  void drawTheoryLines(double mom = 6);

 private:
  void FindPeak(double (&cherenkovreco)[5], double (&spr)[5]);
  int FindPdg(double mom, double cangle);
  void FitRing(double &x0, double &y0, double &theta);
  double FindStartTime(PrtEvent *e);
  double CalcRejection(TH1F *h1, TH1F *h2, double eff);

  PrtRun *frun;
  PrtTools ft;
  int fmaxch, fnpmt, fnpix;
  int fDetectorID;
  double fMomentum, fBboxNum, fPipehAngle, fDphi, fBarPhi;
  int fMethod;
  int fPhysList;
  int fRadiator;
  double fRadiatorL;
  int fStudyId;
  bool fTimeImaging;
  bool fGeomReco;
  bool fRingFit;
  bool fNNet;
  int fCorrType;
  int fCorrLevel;

  TClonesArray *fLut;
  TClonesArray *fTrackInfoArray;

  TFile *fFile;
  TTree *fTree;
  TChain *fChain;

  PrtEvent *fEvent;
  PrtHit fHit;

  // Set the parameters to the default values.
  void SetDefaultParameters();

  // Verbosity level
  int fVerbose;
  int fBatch;
  int nevents;
  TString fInputFile;
  int fp1;
  int fp2;
  int fCor_level;
  double fCor_angle[48] = {0}, fCor_time[48] = {0}, fCor_time_refl[2] = {0};
  double fCorrSpr;
  TString fCorrPath;
  TString fPdfPath;
  TString fNNPath;
  
  TF1 *fFit, *fChromCor;  
  TSpectrum *fSpect;
  double fSigma[5];
  TH1F *hthetac[5];
  TH1F *hthetacd[5];
  TH1F *hnph_gr[5], *hnph_ti[5];
  double fAngle[5];
  TF1 *fFunc[5];
  TH1F *fLnDiffGr[5];
  TH1F *fLnDiffTi[5];
  TH1F *fLnDiffNn[5];
  double fCriticalAngle;
  TString fCorrFile;
  TH1F *fTime[5][99000]; //7000
  // std::vector<std::vector<TH1F*>> fTime;//(5, std::vector<TH1F*>(99000));
  TGraph *fPdf[5][99000];
  TH1F *fPmt_a[28], *fPmt_td[28], *fPmt_tr[28], *fFindTimeA[20], *fHistDiff[3];
  TH1F *fHist1, *fHist2, *fTrackAngle0, *fTrackAngle1, *fTrackAngle2, *fFindTime, *fFindTimeRes;
  TH2F *fDiff, *fHist4, *fHist5, *fdtt, *fdtl, *fdtp, *fhChromL;
  
#ifdef AI
  cppflow::model *fNNmodel;
#endif
  
};

#endif
