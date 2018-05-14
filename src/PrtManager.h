// -----------------------------------------
// PrtManager.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtManager_h
#define PrtManager_h

#include "globals.hh"

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "TClonesArray.h"

#include "PrtEvent.h"
#include "PrtHit.h"
#include "PrtTrackInfo.h"

class PrtManager
{
  static PrtManager* fInstance;
  TFile *fRootFile;
  TTree *fTree;
  PrtEvent *fEvent;
  PrtTrackInfo *fTrackInfo;
  PrtHit *fHit;
  TH1F *fHist;

public:
  PrtManager(G4String outfile, G4int runtype);
  ~PrtManager(){};
  static PrtManager* Instance(G4String outfile="hits.root", G4int runtype=0);
  void Save()             { fRootFile->Write(); }
  void Fill();
  void FillLut();
  void AddEvent(PrtEvent event);
  void AddHit(PrtHit hit);
  void AddTrackInfo(PrtTrackInfo trackinfo);
  PrtEvent* Event(){ return fEvent; }
  
  // Mutators
  void SetRunType(int val){ fRunType = val; }
  void SetPhysList(int val){ fPhysList = val; }
  void SetGeometry(int val){ fGeometry = val; }
  void SetEvType(int val){ fEvType = val; }
  void SetBeamDimension(double val){ fBeamDimension = val; }
  void SetRadiator(int val){ fRadiator = val; }
  void SetLens(int val){ fLens = val; }
  void SetMcpLayout(int val){ fMcpLayout = val; }
  void SetAngle(double val){ fAngle = val; }
  void SetRadiatorL(double val){ fRadiatorL = val; }
  void SetRadiatorW(double val){ fRadiatorW = val; }
  void SetRadiatorH(double val){ fRadiatorH = val; }
  void SetParticle(int val){ fParticle = val; }
  void SetMomentum(TVector3 val){ fMomentum = val; if(fRunType==0 || fRunType==10) fEvent->SetMomentum(fMomentum);}
  void SetCurrentCherenkov(double val){ fCurrentCherenkov = val; }
  void SetShift(double val){ fShift = val; }
  void SetTest1(double val){ fTest1 = val; }
  void SetTest2(double val){ fTest2 = val; }
  void SetDisplayOpt(int val){ fDispalyOpt = val; }
  void SetTimeRes(double val){ fTimeRes = val; }
  void SetMix(int val){fMix = val;}
  
  
  // Accessors
  int GetRunType(){ return fRunType; }
  int GetPhysList(){ return fPhysList; }
  int GetGeometry(){ return fGeometry; }
  int GetEvType(){ return fEvType; }
  double GetBeamDimension(){ return fBeamDimension; }
  int GetRadiator(){ return fRadiator; }
  int GetLens(){ return fLens; }
  int GetMcpLayout(){ return fMcpLayout; }
  double GetAngle(){ return fAngle; }
  double GetRadiatorL(){ return fRadiatorL; }
  double GetRadiatorW(){ return fRadiatorW; }
  double GetRadiatorH(){ return fRadiatorH; }
  int GetParticle(){ return fParticle; }
  TVector3 GetMomentum(){ return fMomentum; }
  double GetCurrentCherenkov(){ return fCurrentCherenkov; }
  double GetShift(){ return fShift; }
  double GetTest1(){ return fTest1; }
  double GetTest2(){ return fTest2; }
  int GetDisplayOpt(){ return fDispalyOpt; }
  double GetTimeRes(){ return fTimeRes; }
  TTree *GetTree(){ return fTree; }
  TString GetOutName(){return fOutName;}
  int GetMix(){ return fMix; }

  
private: 
  int fRunType;
  int fPhysList;
  int fGeometry;
  int fEvType;
  int fRadiator;
  int fLens;
  int fMcpLayout;
  double fAngle;
  double fRadiatorL;
  double fRadiatorW;
  double fRadiatorH;
  int fParticle;
  double fBeamDimension;
  TVector3 fMomentum;
  TClonesArray *fLut;
  TClonesArray *fTrackInfoArray;
  double fCurrentCherenkov;
  double fShift;
  double fTest1;
  double fTest2;
  int fDispalyOpt;
  double fTimeRes;
  TString fOutName;
  int fMix;
  
};

#endif
