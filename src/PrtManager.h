#ifndef PrtManager_h
#define PrtManager_h

#include "globals.hh"

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>

#include "PrtEvent.h"
#include "PrtHit.h"

class PrtManager
{
  static PrtManager* fInstance;
  TFile *fRootFile;
  TTree *fTree;
  PrtEvent *fEvent;
  PrtHit *fHit;
  TH1F *fHist;

public:
  PrtManager(G4String outfile);
  ~PrtManager(){};
  static PrtManager* Instance(G4String outfile="hits.root");
  void Save()             { fRootFile->Write(); }
  void Fill();
  void AddEvent(PrtEvent event);
  void AddHit(PrtHit hit);
  PrtEvent* Event(){ return fEvent; }
  
  // Mutators
  void SetPhysList(int val){ fPhysList = val; }
  void SetGeometry(int val){ fGeometry = val; }
  void SetLens(int val){ fLens = val; }
  void SetAngle(double val){ fAngle = val; }
  void SetParticle(int val){ fParticle = val; }
  void SetMomentum(double val){ fMomentum = val; }
  
  // Accessors
  int GetPhysList(){ return fPhysList; }
  int GetGeometry(){ return fGeometry; }
  int GetLens(){ return fLens; }
  double GetAngle(){ return fAngle; }
  int GetParticle(){ return fParticle; }
  double GetMomentum(){ return fMomentum; }
  
private: 
  int fPhysList;
  int fGeometry;
  int fLens;
  double fAngle;
  int fParticle;
  double fMomentum;
};

#endif
