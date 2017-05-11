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
 : G4VUserPrimaryGeneratorAction(), 
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  //create a messenger for this class
  fGunMessenger = new PrtPrimaryGeneratorMessenger(this);

  //default kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  fParticleK = particleTable->FindParticle("kaon+");
  fParticlePi = particleTable->FindParticle("pi+");

  fParticleGun->SetParticleDefinition(fParticlePi);
  fParticleGun->SetParticleTime(0.0*ns);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0*cm,0.0*cm,0.0*cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  fParticleGun->SetParticleEnergy(500.0*keV);
}

PrtPrimaryGeneratorAction::~PrtPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
}

void PrtPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  G4double x,y,z;
  G4double angle = PrtManager::Instance()->GetAngle();
  G4double radiatorL = PrtManager::Instance()->GetRadiatorL();
  G4double radiatorW = PrtManager::Instance()->GetRadiatorW();
  G4double radiatorH = PrtManager::Instance()->GetRadiatorH();

  if(PrtManager::Instance()->GetMixPiK()){
    if(PrtManager::Instance()->GetParticle()==211 || PrtManager::Instance()->GetParticle()==0){
      fParticleGun->SetParticleDefinition(fParticleK);
      PrtManager::Instance()->SetParticle(321);
    }else{
      fParticleGun->SetParticleDefinition(fParticlePi);
      PrtManager::Instance()->SetParticle(211);
    }
  }
  
  PrtManager::Instance()->AddEvent(PrtEvent());
  if(PrtManager::Instance()->GetRunType() == 0 || PrtManager::Instance()->GetRunType() == 10){ // simulation
    G4ThreeVector vec(0,0,1);
    //G4int id = anEvent->GetEventID()%5;
    // if(id==0)  vec.setTheta(M_PI-110*deg);
    // if(id==0)  vec.setPhi(0*deg);
    // if(id==1)  vec.setTheta(M_PI-30*deg);
    // if(id==1)  vec.setPhi(70*deg);
    // if(id==2)  vec.setTheta(M_PI-140*deg);
    // if(id==2)  vec.setPhi(180*deg);
    // if(id==3)  vec.setTheta(M_PI-70*deg);
    // if(id==3)  vec.setPhi(250*deg);

    // // else{
    if(angle>0)  vec.setTheta(angle);
    else{
      G4double theta = M_PI*G4UniformRand();
      if(theta>140*deg) theta-=40*deg;
      if(theta<20*deg) theta+=20*deg;
      vec.setTheta(M_PI-theta);
      vec.setPhi(2*M_PI*G4UniformRand());
    }
    // // }
   
    fParticleGun->SetParticleMomentumDirection(vec);
  }

  if(PrtManager::Instance()->GetBeamDinsion() == 1){ // smearing and divergence
    G4double sigma = 1*cm;
    z = fParticleGun->GetParticlePosition().z();
    x = G4RandGauss::shoot(0,sigma);
    y = G4RandGauss::shoot(0,sigma);
    
    fParticleGun->SetParticlePosition(G4ThreeVector(x,y,z));
    PrtManager::Instance()->Event()->SetPosition(TVector3(x,y,z));
  }
  if(PrtManager::Instance()->GetRunType() == 1){ // LUT generation
    
    G4double barShift=0; // 390/12./2;
    fParticleGun->SetParticlePosition(G4ThreeVector(0,barShift,radiatorL/2.-0.1));
    G4double angle = -G4UniformRand()*M_PI;
    G4ThreeVector vec(0,0,-1);
    vec.setTheta(acos(G4UniformRand()));
    vec.setPhi(2*M_PI*G4UniformRand());
    
    //  vec.rotateY(-M_PI/2.);
    fParticleGun->SetParticleMomentumDirection(vec);
  }
  if(PrtManager::Instance()->GetRunType() == 5){ // calibration light
    G4double shift = PrtManager::Instance()->GetShift();
    
    fParticleGun->SetParticlePosition(G4ThreeVector(-1250/2.+0.1-shift,0,5+tan(45*M_PI/180.)*shift+25));
    G4double angle = -G4UniformRand()*M_PI;
    G4ThreeVector vec(0,0,1);
    vec.setTheta(acos(G4UniformRand()));
    vec.setPhi(2*M_PI*G4UniformRand());
    
    vec.rotateY(-M_PI/2.);
    fParticleGun->SetParticleMomentumDirection(vec);
  }
  fParticleGun->GeneratePrimaryVertex(anEvent);

  G4ThreeVector dir = fParticleGun->GetParticleMomentumDirection();
  dir *= fParticleGun->GetParticleMomentum();
  PrtManager::Instance()->SetMomentum(TVector3(dir.x(),dir.y(),dir.z()));
}

void PrtPrimaryGeneratorAction::SetOptPhotonPolar()
{
 G4double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}

void PrtPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
 if (fParticleGun->GetParticleDefinition()->GetParticleName()!="opticalphoton")
   {
     G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
       "the particleGun is not an opticalphoton " << 
       fParticleGun->GetParticleDefinition()->GetParticleName()<< G4endl;
     return;
   }

 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton);
 G4double modul2       = product*product;
 
 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
 
 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 fParticleGun->SetParticlePolarization(polar);
}
