#include "PrtPrimaryGeneratorAction.h"
#include "PrtPrimaryGeneratorMessenger.h"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "PrtManager.h"

PrtPrimaryGeneratorAction::PrtPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {

  int n_particle = 1;
  fRun = PrtManager::Instance()->getRun();
  double mom = fRun->getMomentum();
  fRunType = fRun->getRunType();
  fRadiatorL = fRun->getRadiatorL();
  fRadiatorW = fRun->getRadiatorW();
  fRadiatorH = fRun->getRadiatorH();
  
  fParticleGun = new G4ParticleGun(n_particle);

  // create a messenger for this class
  fGunMessenger = new PrtPrimaryGeneratorMessenger(this);

  // default kinematic
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  fParticle[4] = particleTable->FindParticle("proton");
  fParticle[3] = particleTable->FindParticle("kaon+");
  fParticle[2] = particleTable->FindParticle("pi+");
  fParticle[1] = particleTable->FindParticle("mu+");
  fParticle[0] = particleTable->FindParticle("e-");
  fParticleOP = particleTable->FindParticle("opticalphoton");

  fParticleGun->SetParticleDefinition(fParticle[2]);
  fParticleGun->SetParticleTime(0.0 * ns);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0 * cm, 0.0 * cm, -0.5 * cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  fParticleGun->SetParticleMomentum(G4ThreeVector(0, 0, mom * GeV));

  iter = 0;
  fPid = 4;
}

PrtPrimaryGeneratorAction::~PrtPrimaryGeneratorAction() {
  delete fParticleGun;
  delete fGunMessenger;
}

void PrtPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {

  PrtManager::Instance()->addEvent(PrtEvent());
  int pdg = fRun->getPid();
  double theta = 180-fRun->getTheta();
  double phi = fRun->getTheta();
  double zpos = fRun->getBeamZ();
  double ypos = fRun->getBeamX();  
  
  if (pdg > 0) {
    if (pdg == 2212) fPid = 4;
    else if (pdg == 321) fPid = 3;
    else if (pdg == 211) fPid = 2;
    else if (pdg == 10001 && fPid > 2) fPid = 2;
    else if (pdg == 10001 && fPid == 2) fPid = 4;
    else if (pdg == 10002 && fPid > 2) fPid = 2;
    else if (pdg == 10002 && fPid == 2) fPid = 3;

    PrtManager::Instance()->getEvent()->setPid(fPid);
    fParticleGun->SetParticleDefinition(fParticle[fPid]);
  } else {
    fParticleGun->SetParticleDefinition(fParticleOP);
  }

  if (fRunType == 0 || fRunType == 10 || fRunType == 5) { // simulation
    G4ThreeVector vec(0, 0, 1);
    // G4int id = anEvent->GetEventID()%5;
    //  if(id==0)  vec.setTheta(M_PI-110*deg);
    //  if(id==0)  vec.setPhi(0*deg);
    //  if(id==1)  vec.setTheta(M_PI-30*deg);
    //  if(id==1)  vec.setPhi(70*deg);
    //  if(id==2)  vec.setTheta(M_PI-140*deg);
    //  if(id==2)  vec.setPhi(180*deg);
    //  if(id==3)  vec.setTheta(M_PI-70*deg);
    //  if(id==3)  vec.setPhi(250*deg);

    double trackresolution = fRun->getBeamSize();
    if (theta > 0) {
      if (trackresolution < 0.00001) {
        vec.setTheta(theta*TMath::DegToRad());
      } else {
        // smear track resolution
        G4ThreeVector vec0 = vec;
        vec0.setTheta(theta);
        vec.setTheta(G4RandGauss::shoot(theta, trackresolution));
        vec.rotate(2 * M_PI * G4UniformRand(), vec0);
      }
    } else {
      double theta = M_PI * G4UniformRand();
      theta = acos((cos(30 * deg) - cos(150 * deg)) * G4UniformRand() + cos(150 * deg));

      G4ThreeVector vec0 = vec;
      vec0.setTheta(M_PI - theta);

      theta = G4RandGauss::shoot(theta, trackresolution);
      vec.setTheta(M_PI - theta);
      vec.rotate(2 * M_PI * G4UniformRand(), vec0);

      if (fRunType != 5) PrtManager::Instance()->getEvent()->setMomentum(TVector3(vec.x(),vec.y(),vec.z()));
    }

    if (fRun->getEv() == 1) {
      if (fRunType == 5) {
        vec.setTheta(theta + G4UniformRand() * 0.09 - 0.045);
        vec.rotateZ(G4UniformRand() * 0.032);
      }
      if (fRunType == 0)
        vec.rotateZ(fRun->getTest1()); // 0.016 hits the middle of the BaBar bar
    }

    fParticleGun->SetParticlePosition(G4ThreeVector(0, ypos, zpos));
    if (fRunType != 5) PrtManager::Instance()->getEvent()->setPosition(TVector3(0, ypos, zpos));

    fParticleGun->SetParticleMomentumDirection(vec);
  }

  if (fRunType == 1) { // LUT generation

    double barShift = 0; // 390/12./2;
    if (fRun->getEv() == 1) barShift = 0.5 * 35;

    fParticleGun->SetParticlePosition(G4ThreeVector(0, barShift, fRadiatorL / 2. - 0.5));
    G4ThreeVector vec(0, 0, -1);
    vec.setTheta(acos(G4UniformRand()));
    vec.setPhi(2 * M_PI * G4UniformRand());

    fParticleGun->SetParticleMomentumDirection(vec);
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);

  G4ThreeVector dir = fParticleGun->GetParticleMomentumDirection();
  dir *= fParticleGun->GetParticleMomentum();
  PrtManager::Instance()->getEvent()->setMomentum(TVector3(dir.x(),dir.y(),dir.z()));
}

void PrtPrimaryGeneratorAction::SetOptPhotonPolar() {
  double angle = G4UniformRand() * 360.0 * deg;
  SetOptPhotonPolar(angle);
}

void PrtPrimaryGeneratorAction::SetOptPhotonPolar(double angle) {

  if (fParticleGun->GetParticleDefinition()->GetParticleName() != "opticalphoton") {
    G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
              "the particleGun is not an opticalphoton "
           << fParticleGun->GetParticleDefinition()->GetParticleName() << G4endl;
    return;
  }

  G4ThreeVector normal(1., 0., 0.);
  G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
  G4ThreeVector product = normal.cross(kphoton);
  double modul2 = product * product;

  G4ThreeVector e_perpend(0., 0., 1.);
  if (modul2 > 0.) e_perpend = (1. / std::sqrt(modul2)) * product;
  G4ThreeVector e_paralle = e_perpend.cross(kphoton);

  G4ThreeVector polar = std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
  fParticleGun->SetParticlePolarization(polar);
}
