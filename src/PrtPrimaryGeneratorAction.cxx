#include "PrtPrimaryGeneratorAction.h"
#include "PrtPrimaryGeneratorMessenger.h"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4RunManager.hh"

#include "PrtManager.h"


PrtPrimaryGeneratorAction::PrtPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {

  int n_particle = 1;
  fRun = PrtManager::Instance()->getRun();
  double mom = fRun->getMomentum();
  fRunType = fRun->getRunType();
  fTracking = fRun->getTrackingResTheta();

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
  fCurrentEvent = -1;  
}

PrtPrimaryGeneratorAction::~PrtPrimaryGeneratorAction() {
  delete fParticleGun;
  delete fGunMessenger;
}

void PrtPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {

  PrtManager::Instance()->addEvent(PrtEvent());
  int pdg = fabs(fRun->getPid());
  double theta = (180 - fRun->getTheta()) * TMath::DegToRad();
  // double theta2 = (180 - fRun->getTest1()) * TMath::DegToRad();
  double phi = fRun->getPhi() * TMath::DegToRad();
  double zpos = fRun->getBeamZ();
  double ypos = fRun->getBeamX();
  fRadiatorL = fRun->getRadiatorL();
  fRadiatorW = fRun->getRadiatorW();
  fRadiatorH = fRun->getRadiatorH();
  fGeomType = fRun->getGeometry();
  fCurrentEvent++;

  if (pdg != 0) {
    if (pdg == 2212) fPid = 4;
    else if (pdg == 11) fPid = 0;
    else if (pdg == 13) fPid = 1;
    else if (pdg == 211) fPid = 2;
    else if (pdg == 321) fPid = 3;
    else if (pdg == 10000 && fPid != 2) fPid = 2;
    else if (pdg == 10000) fPid = 0;
    else if (pdg == 10001 && fPid != 2) fPid = 2;
    else if (pdg == 10001) fPid = 1;
    else if (pdg == 10003 && fPid != 2) fPid = 2;
    else if (pdg == 10003) fPid = 3;
    else if (pdg == 10004 && fPid != 2) fPid = 2;
    else if (pdg == 10004) fPid = 4;
    else if (pdg == 10005 && fPid != 3) fPid = 3;
    else if (pdg == 10005) fPid = 4;
    fParticleGun->SetParticleDefinition(fParticle[fPid]);
  } else {
    fParticleGun->SetParticleDefinition(fParticleOP);
  }
  PrtManager::Instance()->getEvent()->setPid(fPid);
  G4ThreeVector vec(0, 0, 1);
  if (fRunType == 0 || fRunType == 10 || fRunType == 5) { // simulation
    
    if (fGeomType < 2) ypos = 0.5 * fRadiatorW;
    
    fParticleGun->SetParticlePosition(G4ThreeVector(0, ypos, zpos));

    // // second track
    // vec.setTheta(theta2);
    // vec.setPhi(phi);
    // fParticleGun->SetParticleMomentumDirection(vec);
    // fParticleGun->GeneratePrimaryVertex(anEvent);

    if(fCurrentEvent > 10000) fTracking = 0.0005; // small smearing for PDFs

    if (theta > 0 && theta < M_PI) {

      vec.setTheta(theta);
      vec.setPhi(phi);
      
      // if (fTracking < 0.00001) {
      //   vec.setTheta(theta);
      //   vec.setPhi(phi);
      // } else if (fTracking < 10) {

      //   // // smear track resolution
      //   // vec.setTheta(theta);
      //   // vec.setPhi(phi);
      //   // G4ThreeVector vec0 = vec;
      //   // vec.setTheta(G4RandGauss::shoot(theta, fTracking));
      //   // vec.rotate(M_PI * G4UniformRand(), vec0.unit());

      //   vec.setTheta(G4RandGauss::shoot(theta, fTracking));
      //   vec.setPhi(G4RandGauss::shoot(phi, fTracking));
	
      //   std::cout << "track resolution dtheta = " << fTracking << " mrad, dphi = " << fTracking
      //             << " mrad" << std::endl;
      // // } else {	
      // //   double dtheta = 0.001 * get_res(grtheta, theta, 0.001 * fParticleGun->GetParticleMomentum());
      // //   double dphi = 0.001 * get_res(grphi, theta, 0.001 * fParticleGun->GetParticleMomentum());
      // // 	fRun->setTrackingResTheta(dtheta);
      // // 	fRun->setTrackingResPhi(dphi);       

      // //   vec.setTheta(G4RandGauss::shoot(theta, dtheta));
      // //   vec.setPhi(G4RandGauss::shoot(phi, dphi));
      // // 	std::cout << "track resolution dtheta = " << dtheta <<  " mrad, dphi = " << dphi << " mrad" << std::endl;
      // }
    } else {
      theta = M_PI * G4UniformRand();
      theta = acos((cos(30 * deg) - cos(150 * deg)) * G4UniformRand() + cos(150 * deg));

      G4ThreeVector vec0 = vec;
      vec0.setTheta(M_PI - theta);

      theta = G4RandGauss::shoot(theta, fTracking);
      vec.setTheta(M_PI - theta);
      vec.rotate(2 * M_PI * G4UniformRand(), vec0);
      // vec.setPhi(2 * M_PI * G4UniformRand());
    }

    if (fRun->getEv() == 1) {
      if (fRunType == 5) {
        vec.setTheta(theta + G4UniformRand() * 0.09 - 0.045);
        vec.rotateZ(G4UniformRand() * 0.032);
      }
      if (fRunType == 0) vec.rotateZ(0.016); // 0.016 hits the middle of the BaBar bar
    }  
    fParticleGun->SetParticleMomentumDirection(vec);
  }

  if (fRunType == 1) { // LUT generation

    double barShift = 0.5 * fRadiatorW; // 390/12./2;
    if (fRun->getEv() == 1) barShift = 0.5 * 35;
    if (fGeomType == 2) barShift = 0;

    fParticleGun->SetParticlePosition(G4ThreeVector(0, barShift, 0.5 * fRadiatorL - 0.5));
    G4ThreeVector v(0, 0, -1);
    v.setTheta(acos(G4UniformRand()));
    v.setPhi(2 * M_PI * G4UniformRand());

    fParticleGun->SetParticleMomentumDirection(v);
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);

  G4ThreeVector dir = fParticleGun->GetParticleMomentumDirection();
  dir *= fParticleGun->GetParticleMomentum();

  PrtManager::Instance()->getEvent()->setMomentum(TVector3(dir.x(), dir.y(), dir.z()));
  PrtManager::Instance()->setMomentum(TVector3(dir.x(), dir.y(), dir.z())); 

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
