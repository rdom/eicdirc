// -----------------------------------------
// PrtLutReco.h
//
// Created on: 13.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------
// Class for reconstruction in DIRC using look-up table method
 
#ifndef PrtLutReco_h
#define PrtLutReco_h 1

#include "PrtEvent.h"
#include "PrtHit.h"

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TSpectrum.h"

class PrtLutReco{

public:

  // Standard constructors
  PrtLutReco(TString infile, TString lutfile, int verbose=0);

  // Destructor
  ~PrtLutReco();
  void Run(int start=0, int end=0);
  void drawTheoryLines(double mom=6);

private:
  bool FindPeak(double& cherenkovreco, double& spr);
  int FindPdg(double mom, double cangle);
  void FitRing(double& x0, double& y0, double& theta);
  double FindStartTime(PrtEvent *e);
  int fDetectorID;  
  double fBboxNum,fPipehAngle,fDphi,fBarPhi;
  double fSigma;
  
  TClonesArray *fLut;
  TClonesArray *fTrackInfoArray;

  TFile *fFile; 
  TTree *fTree;
  TChain *fChain;

  PrtEvent* fEvent;
  PrtHit fHit;
  
  // Set the parameters to the default values.
  void SetDefaultParameters();
  
  // Verbosity level
  int fVerbose;
  int nevents;
  int fMethod;
  TString fInputFile;
  int fRpid;

  TF1 *fFit;
  TSpectrum *fSpect;
  TH1F *hthetac[5];
  TH1F *hthetacd[5];
  TH1F *hnph[5];
  double fAngle[5];
  TF1 *fFunc[5];
  TH1F *fLnDiff[5];
  double fCriticalAngle;
  TString fCorrFile;
};

#endif
