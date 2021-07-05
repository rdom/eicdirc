#ifndef PrtOpBoundaryProcess_h
#define PrtOpBoundaryProcess_h

#include "globals.hh"
#include "G4OpBoundaryProcess.hh"

class PrtOpBoundaryProcess : public G4OpBoundaryProcess {
 public:
  PrtOpBoundaryProcess();
  ~PrtOpBoundaryProcess(){};

 public:
  G4VParticleChange *PostStepDoIt(const G4Track &aTrack, const G4Step &aStep);

 private:
  int fLensId;
  int fRunType;
  int fEvType;
};

#endif /*PrtOpBoundaryProcess_h*/
