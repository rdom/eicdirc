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

#include "PrtRun.h"
#include "TGraph.h"

class G4ParticleGun;
class G4Event;
class PrtPrimaryGeneratorMessenger;

class PrtPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  PrtPrimaryGeneratorAction();
  virtual ~PrtPrimaryGeneratorAction();

 public:
  virtual void GeneratePrimaries(G4Event *);

  void SetOptPhotonPolar();
  void SetOptPhotonPolar(G4double);

 private:
  PrtRun *fRun;
  int fRunType, fGeomType, fPid, fPdg, iter;
  double fRadiatorL, fRadiatorW, fRadiatorH;
  double fTracking;
  G4ParticleGun *fParticleGun;
  G4ParticleDefinition *fParticleOP, *fParticle[9];
  PrtPrimaryGeneratorMessenger *fGunMessenger;
  int fCurrentEvent;
};

#endif
