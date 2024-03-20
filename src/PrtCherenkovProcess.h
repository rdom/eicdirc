#ifndef PrtCherenkovProcess_h
#define PrtCherenkovProcess_h

#include "globals.hh"
#include "G4Cerenkov.hh"

class PrtCherenkovProcess : public G4Cerenkov
{
public:
  PrtCherenkovProcess(const G4String& processName, G4ProcessType type = fElectromagnetic);
  ~PrtCherenkovProcess(){};

public:
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

private:
  int fLensId;

  G4double GetAverageNumberOfPhotons(const G4double charge,
				     const G4double beta,
				     const G4Material *aMaterial,
				     G4MaterialPropertyVector* Rindex) const;
};

#endif /*PrtCherenkovProcess_h*/
