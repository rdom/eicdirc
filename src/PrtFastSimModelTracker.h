#ifndef PrtFastSimModelTracker_h
#define PrtFastSimModelTracker_h 1

#include "G4VFastSimulationModel.hh"
#include "G4Step.hh"
#include "G4Navigator.hh"

#include "PrtRun.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TH2F.h"

class PrtFastSimModelTracker : public G4VFastSimulationModel {
 public:
  PrtFastSimModelTracker(G4String aModelName, G4Region *aEnvelope);

  ~PrtFastSimModelTracker();

  virtual G4bool IsApplicable(const G4ParticleDefinition &aParticle);

  virtual G4bool ModelTrigger(const G4FastTrack &aFastTrack);

  virtual void DoIt(const G4FastTrack &aFastTrack, G4FastStep &aFastStep);

  double get_res(TGraph *gg[3], double theta, double mom);
  
 private:
  PrtRun *fRun;
  double fTrackingRes;
  TGraph *grtheta[3], *grphi[3];
  TH2F *fTrMapTheta, *fTrMapPhi;
};

#endif
