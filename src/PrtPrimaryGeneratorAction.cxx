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
#include "G4RandomDirection.hh"

#include "PrtManager.h"


PrtPrimaryGeneratorAction::PrtPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {

  int n_particle = 1;
  fRun = PrtManager::Instance()->getRun();
  double mom = fRun->getMomentum();
  fRunType = fRun->getRunType();
  fTracking = fRun->getTrackingResTheta();

  fParticleGun = new G4ParticleGun(n_particle);
  fGunMessenger = new PrtPrimaryGeneratorMessenger(this);

  // default kinematic
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  
  fParticle[0] = particleTable->FindParticle("e-");
  fParticle[1] = particleTable->FindParticle("mu-");
  fParticle[2] = particleTable->FindParticle("pi+");
  fParticle[3] = particleTable->FindParticle("kaon+");
  fParticle[4] = particleTable->FindParticle("proton");
  
  fParticle[5] = particleTable->FindParticle("e+");
  fParticle[6] = particleTable->FindParticle("mu+");
  fParticle[7] = particleTable->FindParticle("pi-");
  fParticle[8] = particleTable->FindParticle("kaon-");  

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
  int pdg = fRun->getPid();
  double theta = (180 - fRun->getTheta()) * TMath::DegToRad();
  // double theta2 = (180 - fRun->getTest1()) * TMath::DegToRad();
  double phi = fRun->getPhi() * TMath::DegToRad();
  double zpos = fRun->getBeamZ();
  double ypos = fRun->getBeamX();
  fRadiatorL = fRun->getRadiatorL();
  fRadiatorW = fRun->getRadiatorW();
  fRadiatorH = fRun->getRadiatorH();
  fGeomType = fRun->getGeometry();
  fField = fRun->getField();
  fCurrentEvent++;

  if (pdg != 0) {
    if (pdg == 2212) fPid = 4;
    else if (pdg == 11) fPid = 0;
    else if (pdg == 13) fPid = 1;
    else if (pdg == 211) fPid = 2;
    else if (pdg == 321) fPid = 3;
    else if (pdg == -11) fPid = 5;
    else if (pdg == -13) fPid = 6;
    else if (pdg == -211) fPid = 7;
    else if (pdg == -321) fPid = 8;    
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
    
    // // second track
    // vec.setTheta(theta2);
    // vec.setPhi(phi);
    // fParticleGun->SetParticleMomentumDirection(vec);
    // fParticleGun->GeneratePrimaryVertex(anEvent);

    // WJL target bar code -------------------------
    int kTargetBar = -1;
    double fTargetBarY = 0;
    double fBarsGap = 0.15; //!!!! ALSO DEFINED IN PrtDetectorConstruction
    double radius = 770.5;  // ePIC
    double acharge = fParticleGun->GetParticleCharge();
    double momentum = fParticleGun->GetParticleMomentum();
    if (fRun->getPhi() >= 990 && fRun->getPhi() <= 999) {
      //
      kTargetBar = fRun->getPhi() - 990; // [0,9]
      fTargetBarY = ((double)(kTargetBar - 5)) * fRadiatorW +
                    fRadiatorW / 2.; // mm  !!!! ASSUMES TEN BARS, also IGNORES GAPS!!!!
      double partheta = (180. - fRun->getTheta()); // deg
      double distheta = fabs(partheta - 90.);      // deg
      // field=2 MARCO 1.7T, theta scalings to handle radial component of field
      double Bval = (-1.75 + (1.75 - 1.6) * (distheta / 60.)) * CLHEP::tesla;
      double pt = momentum * sin(theta);
      double R = 1000. * (pt / CLHEP::GeV) / 0.3 / fabs(Bval / CLHEP::tesla); // meters->mm
      if (fField == 0) R = 1E12;                                              // no field
      double x1 = 0;                                                          // primary vertex
      double y1 = 0;                                                          // primary vertex
      double x2 = radius + fRadiatorH / 2.;                                   // target point
      double y2 = fTargetBarY; // target point !!!! HERE ASSUMES TEN BARS PER BOX
      double xa = (x2 - x1) / 2.;
      double ya = (y2 - y1) / 2.;
      double x0 = x1 + xa;
      double y0 = y1 + ya;
      double a = sqrt(xa * xa + ya * ya);
      double b = sqrt(R * R - a * a);
      double x3 = x0 + b * ya / a;
      double y3 = y0 - b * xa / a;
      double x4 = x0 - b * ya / a;
      double y4 = y0 + b * xa / a;
      // reverse order to get phi angle of perpendicular to ray from PV to (x3,y3), NEG
      double ang3 = atan2(x3, y3);
      // reverse order to get phi angle of perpendicular to ray from PV to (x4,y4), POS
      double ang4 = atan2(x4, y4);
      if (acharge > 0.) phi = -ang4;      // pos
      else if (acharge < 0.) phi = M_PI - ang3; // neg
      //
    } else {
      ypos = 0.5 * fRadiatorW;
      phi = fRun->getPhi() * TMath::DegToRad();
    }
    // WJL end target bar code -------------------------

    fParticleGun->SetParticlePosition(G4ThreeVector(0, ypos, zpos));
     
    if (theta > 0 && theta < M_PI) {
      vec.setTheta(theta);
      vec.setPhi(phi);      
    } else {
      // theta = M_PI * G4UniformRand();
      theta = acos((cos(25 * deg) - cos(160 * deg)) * G4UniformRand() + cos(160 * deg));
      theta = M_PI - theta;
	    
      vec.setTheta(theta);
      vec.setPhi(2 * M_PI * G4UniformRand());

      // 3 tracks for visualization
      vec.setTheta(25 * deg);
      vec.setPhi(0.5 * deg);
      fParticleGun->SetParticleMomentumDirection(vec);
      fParticleGun->GeneratePrimaryVertex(anEvent);
      
      vec.setTheta(120 * deg);
      vec.setPhi(M_PI / 3. + 7 *deg);
      fParticleGun->SetParticleMomentumDirection(vec);
      fParticleGun->GeneratePrimaryVertex(anEvent);

      vec.setTheta(150 * deg);
      vec.setPhi(-84 * deg);
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
    int targetbar = 5;
    int bars = fRun->getRadiator();
    if (fRun->getPhi() >= 990 && fRun->getPhi() <= 999)  targetbar = fRun->getPhi() - 990; // [0,9]
    double barShift = (-0.5 * bars + targetbar) * fRadiatorW + 0.5 * fRadiatorW;

    if (fRun->getEv() == 1) barShift = 0.5 * 35;
    if (fGeomType == 2) barShift = 0;

    fParticleGun->SetParticlePosition(G4ThreeVector(-200, barShift, 0.5 * fRadiatorL + 630 - 0.2));
    G4ThreeVector v(0, 0, -1);
    v.setTheta(acos(G4UniformRand()));
    v.setPhi(2 * M_PI * G4UniformRand());

    fParticleGun->SetParticlePolarization(G4RandomDirection());
    fParticleGun->SetParticleMomentumDirection(v);
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);
  G4ThreeVector dir = fParticleGun->GetParticleMomentumDirection();
  dir *= fParticleGun->GetParticleMomentum();

  // PrtManager::Instance()->getEvent()->setMomentum(TVector3(dir.x(), dir.y(), dir.z()));
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
