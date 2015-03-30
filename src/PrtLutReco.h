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
  PrtLutReco(TString infile, TString lutfile);

  // Destructor
  ~PrtLutReco();
  void Run(Int_t start=0, Int_t end=0);

private:
  Bool_t FindPeak(Double_t& cherenkovreco, Double_t& spr,Int_t a);
  Int_t FindPdg(Double_t mom, Double_t cangle);
  void FitRing(Double_t& x0, Double_t& y0, Double_t& theta);
  Int_t fDetectorID;  
  Double_t fBboxNum,fPipehAngle,fDphi,fBarPhi;

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
  Int_t fVerbose;
  Int_t nevents;
  TString fInputFile;
  TH1F *fHist;
  TF1 *fFit;
  TSpectrum *fSpect;

};

#endif
