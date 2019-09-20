// -----------------------------------------
// PrtPrimaryGeneratorAction class
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtPrimaryGeneratorAction_h
#define PrtPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class PrtPrimaryGeneratorMessenger;

class PrtPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrtPrimaryGeneratorAction();
  virtual ~PrtPrimaryGeneratorAction();

public:
  virtual void GeneratePrimaries(G4Event*);

  void SetOptPhotonPolar();
  void SetOptPhotonPolar(G4double);

private:
  G4ParticleGun* fParticleGun;
  G4ParticleDefinition* fParticleK;
  G4ParticleDefinition* fParticlePi;
  G4ParticleDefinition* fParticleE;
  G4ParticleDefinition* fParticleMu;
  PrtPrimaryGeneratorMessenger* fGunMessenger;
};

#endif
