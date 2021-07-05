// -----------------------------------------
// PrtManager.h
//
// author  : r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtManager_h
#define PrtManager_h

#include "globals.hh"

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "TClonesArray.h"

#include "PrtRun.h"
#include "PrtEvent.h"
#include "PrtHit.h"
#include "PrtLutNode.h"

class PrtManager {
  static PrtManager *fInstance;
  TFile *fRootFile;
  TTree *fRunTree, *fTree;
  PrtRun *fRun;
  PrtEvent *fEvent;
  PrtHit *fHit;
  TH1F *fHist;

 public:
  PrtManager(TString outfile, PrtRun *run);
  ~PrtManager(){};
  static PrtManager *Instance(TString outfile = "hits.root", PrtRun *run = nullptr);
  void save();
  void fill();
  void fillLut();
  void addEvent(PrtEvent event);
  void addHit(PrtHit hit, TVector3 localpos, TVector3 vertex = TVector3(0, 0, 0));
  void addHit(PrtHit hit);
  PrtEvent *getEvent() { return fEvent; }

  // mutators
  void setRun(PrtRun *v) { fRun = v; }
  void setMomentum(TVector3 v) { fMomentum = v; }
  void setDisplayOpt(int v) { fDisplayOpt = v; }

  // accessors
  PrtRun *getRun() { return fRun; }
  TString getOutName() { return fOutName; }
  int getDisplayOpt(){ return fDisplayOpt; }
  TTree *getTree(){ return fTree; }
  
 private:
  TString fOutName;
  int fRunType;
  int fDisplayOpt;
  TClonesArray *fLut;
  
  TVector3 fMomentum;
  TVector3 fnX1;
  TVector3 fnY1;
  double fCriticalAngle;
};

#endif
