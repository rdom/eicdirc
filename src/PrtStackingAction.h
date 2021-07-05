#ifndef PrtStackingAction_H
#define PrtStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

#include "TGraph.h"
#include "TRandom.h"

class PrtStackingAction : public G4UserStackingAction {
 public:
  PrtStackingAction();
  virtual ~PrtStackingAction();

 public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track *aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();

 private:
  G4int fScintillationCounter;
  G4int fCerenkovCounter;
  int fRunType;
  TGraph *fDetEff[2];
};

#endif
