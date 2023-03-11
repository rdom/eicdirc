#ifndef PrtFastSimModelTracker_h
#define PrtFastSimModelTracker_h 1

#include "G4VFastSimulationModel.hh"
#include "G4Step.hh"
#include "G4Navigator.hh"

#include "PrtRun.h"
#include "TGraph.h"

class PrtFastSimModelTracker : public G4VFastSimulationModel {
 public:
  PrtFastSimModelTracker(G4String aModelName, G4Region *aEnvelope);

  ~PrtFastSimModelTracker();

  /// Checks if this model should be applied to this particle type.
  /// @param aParticle A particle definition (type).
  virtual G4bool IsApplicable(const G4ParticleDefinition &aParticle);

  /// Checks if the model should be applied taking into account the kinematics
  /// of a track.
  /// @param aFastTrack A track.
  virtual G4bool ModelTrigger(const G4FastTrack &aFastTrack);

  /// Calculates the final position (at the outer boundary of the tracking detector)
  /// of a particle with the momentum at the entrance of the tracking detector.
  /// Smears the particle momentum and saves it, together with the tracking detector
  /// resolution and efficiency to the PrtPrimaryParticleInformation.
  /// @param aFastTrack A track.
  /// @param aFastStep A step.
  virtual void DoIt(const G4FastTrack &aFastTrack, G4FastStep &aFastStep);

  double get_res(TGraph *gg[3], double theta, double mom);
  
 private:
  PrtRun *fRun;
  double fTrackingRes;
  TGraph *grtheta[3], *grphi[3];
};

#endif
