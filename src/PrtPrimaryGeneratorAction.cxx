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


  // set tracking resolution

  double mombins[] = {0.75, 1.25, 1.75, 2.5,  3.5,  4.5,  5.5,  6.5,  7.5,  8.5, 9.5,
                      10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5};

  double rtheta_0[] = {5.06794,  3.86015,  2.15975,  1.70208,  1.13201,  0.949914, 0.722017,
                       0.627616, 0.601248, 0.482372, 0.456024, 0.424919, 0.38758,  0.330883,
                       0.314847, 0.298489, 0.276163, 0.282891, 0.287801, 0.251317, 0.232265};
  double rtheta_1[] = {4.30083,  2.28747,  1.74149,  1.2333,   0.884043, 0.679373, 0.583677,
                       0.484924, 0.424464, 0.376458, 0.326366, 0.30053,  0.280166, 0.267392,
                       0.24391,  0.222152, 0.20277,  0.199871, 0.193495, 0.176591, 0.17664};
  double rtheta_2[] = {5.45008,  4.13439,  2.78373,  1.94307,  1.65759,  1.24321,  1.03413,
                       0.844061, 0.764819, 0.697782, 0.591453, 0.537856, 0.521389, 0.434031,
                       0.446836, 0.380808, 0.370561, 0.358088, 0.343058, 0.297783, 0.315344};

  double rphi_0[] = {18.6745,  12.426,   3.19723,  2.67958,  1.97109,  1.82281,  1.44746,
                     1.11069,  1.08232,  1.01937,  0.895427, 0.787919, 0.70771,  0.707664,
                     0.614802, 0.594137, 0.525847, 0.560032, 0.532071, 0.466807, 0.47712};
  double rphi_1[] = {6.54607,  3.47465,  2.92379,  2.04704,  1.48889,  1.17779,  0.912445,
                     0.812745, 0.697406, 0.629205, 0.544928, 0.495619, 0.460813, 0.444274,
                     0.402091, 0.362215, 0.361622, 0.334873, 0.328851, 0.301573, 0.29414};
  double rphi_2[] = {16,       11.1076,  5.4352,   3.55541, 2.73494,  2.07148,  1.98559,
                     1.49008,  1.38194,  1.26842,  1.14616, 1.02384,  1.01611,  0.943076,
                     0.871505, 0.746419, 0.735032, 0.68612, 0.670443, 0.631896, 0.575671};

  grtheta[0] = new TGraph(21, mombins, rtheta_0);
  grtheta[1] = new TGraph(21, mombins, rtheta_1);
  grtheta[2] = new TGraph(21, mombins, rtheta_2);

  grphi[0] = new TGraph(21, mombins, rphi_0);
  grphi[1] = new TGraph(21, mombins, rphi_1);
  grphi[2] = new TGraph(21, mombins, rphi_2);
  
}

double PrtPrimaryGeneratorAction::get_res(TGraph *g[3], double theta, double mom){

  double x[3], y[3];
  x[0] = 2.0 * atan(exp(-1 * 1.5)) * TMath::RadToDeg();  // 25
  x[1] = 2.0 * atan(exp(0)) * TMath::RadToDeg();         // 90
  x[2] = 2.0 * atan(exp(-1 * -1.5)) * TMath::RadToDeg(); // 155

  y[0] = g[0]->Eval(mom);
  y[1] = g[1]->Eval(mom);
  y[2] = g[2]->Eval(mom);

  TGraph * gr = new TGraph(3,x,y);    
  return gr->Eval(theta);;
};

PrtPrimaryGeneratorAction::~PrtPrimaryGeneratorAction() {
  delete fParticleGun;
  delete fGunMessenger;
}

void PrtPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {

  PrtManager::Instance()->addEvent(PrtEvent());
  int pdg = fabs(fRun->getPid());
  double theta = (180 - fRun->getTheta()) * TMath::DegToRad();
  double theta2 = (180 - fRun->getTest1()) * TMath::DegToRad();
  double phi = fRun->getPhi() * TMath::DegToRad();
  double zpos = fRun->getBeamZ();
  double ypos = fRun->getBeamX();
  fRadiatorL = fRun->getRadiatorL();
  fRadiatorW = fRun->getRadiatorW();
  fRadiatorH = fRun->getRadiatorH();
  fGeomType = fRun->getGeometry();

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

    double trackresolution = fRun->getBeamSize();

    if (theta > 0 && theta < M_PI) {
      if (trackresolution < 0.00001) {
        vec.setTheta(theta);
        vec.setPhi(phi);
      } else if (trackresolution < 10) {
        // smear track resolution
        G4ThreeVector vec0 = vec;
        vec0.setTheta(theta);
        vec0.setPhi(phi);
        vec.setTheta(G4RandGauss::shoot(theta, trackresolution));
        vec.setPhi(phi);
        vec.rotate(2 * M_PI * G4UniformRand(), vec0);
      } else {
	std::cout << "fParticleGun->GetParticleMomentum() " << fParticleGun->GetParticleMomentum() * MeV << std::endl;
	
        double dtheta = get_res(grtheta, theta, 0.001 * fParticleGun->GetParticleMomentum());
        double dphi = get_res(grphi, theta, 0.001 * fParticleGun->GetParticleMomentum());
	std::cout << "track resolution dtheta = " << dtheta <<  " mrad, dphi = " << dphi << "mrad" << std::endl;

        vec.setTheta(G4RandGauss::shoot(theta, dtheta));
        vec.setPhi(G4RandGauss::shoot(phi, dphi));
      }
    } else {
      theta = M_PI * G4UniformRand();
      theta = acos((cos(30 * deg) - cos(150 * deg)) * G4UniformRand() + cos(150 * deg));

      G4ThreeVector vec0 = vec;
      vec0.setTheta(M_PI - theta);

      theta = G4RandGauss::shoot(theta, trackresolution);
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
    G4ThreeVector vec(0, 0, -1);
    vec.setTheta(acos(G4UniformRand()));
    vec.setPhi(2 * M_PI * G4UniformRand());

    fParticleGun->SetParticleMomentumDirection(vec);
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
