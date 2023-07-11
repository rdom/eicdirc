#include "PrtCherenkovProcess.h"

#include "PrtManager.h"

#include "G4Poisson.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpticalPhoton.hh"

PrtCherenkovProcess::PrtCherenkovProcess(const G4String &processName, G4ProcessType type)
  : G4Cerenkov(processName, type) {}

G4VParticleChange *PrtCherenkovProcess::PostStepDoIt(const G4Track &aTrack, const G4Step &aStep)

// This routine is called for each tracking Step of a charged particle
// in a radiator. A Poisson-distributed number of photons is generated
// according to the Cerenkov formula, distributed evenly along the track
// segment and uniformly azimuth w.r.t. the particle direction. The
// parameters are then transformed into the Master Reference System, and
// they are added to the particle change.

{
  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
  const G4Material *aMaterial = aTrack.GetMaterial();

  G4StepPoint *pPreStepPoint = aStep.GetPreStepPoint();
  G4StepPoint *pPostStepPoint = aStep.GetPostStepPoint();

  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  G4double t0 = pPreStepPoint->GetGlobalTime();

  G4MaterialPropertiesTable *MPT = aMaterial->GetMaterialPropertiesTable();
  if (!MPT) return pParticleChange;

  G4MaterialPropertyVector *Rindex = MPT->GetProperty(kRINDEX);
  if (!Rindex) return pParticleChange;

  G4double charge = aParticle->GetDefinition()->GetPDGCharge();
  G4double beta = (pPreStepPoint->GetBeta() + pPostStepPoint->GetBeta()) * 0.5;

  G4double MeanNumberOfPhotons = GetAverageNumberOfPhotons(charge, beta, aMaterial, Rindex);

  if (MeanNumberOfPhotons <= 0.0) {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }

  G4double step_length = aStep.GetStepLength();
  MeanNumberOfPhotons = MeanNumberOfPhotons * step_length;
  G4int fNumPhotons = (G4int)G4Poisson(MeanNumberOfPhotons);
  
  if (fNumPhotons <= 0) { //==!fStackingFlag)
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }

  G4double Pmin = Rindex->Energy(0);
  G4double Pmax = Rindex->GetMaxEnergy();

  // d/ monochromatic photons
  // fNumPhotons = 13; // scaled for 30 deg
  Pmin = 3.18 * 1E-6;
  Pmax = Pmin;
  
  ////////////////////////////////////////////////////////////////
  aParticleChange.SetNumberOfSecondaries(fNumPhotons);

  if (1) { //fTrackSecondariesFirst
    if (aTrack.GetTrackStatus() == fAlive) aParticleChange.ProposeTrackStatus(fSuspend);
  }

  ////////////////////////////////////////////////////////////////

  G4double dp = Pmax - Pmin;

  G4double nMax = Rindex->GetMaxValue();
  G4double BetaInverse = 1. / beta;

  G4double maxCos = BetaInverse / nMax;
  G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);

  G4double beta1 = pPreStepPoint->GetBeta();
  G4double beta2 = pPostStepPoint->GetBeta();

  G4double MeanNumberOfPhotons1 = GetAverageNumberOfPhotons(charge, beta1, aMaterial, Rindex);
  G4double MeanNumberOfPhotons2 = GetAverageNumberOfPhotons(charge, beta2, aMaterial, Rindex);

  for (G4int i = 0; i < fNumPhotons; ++i) {
    // Determine photon energy
    G4double rand;
    G4double sampledEnergy, sampledRI;
    G4double cosTheta, sin2Theta;

    // sample an energy
    do {
      rand = G4UniformRand();
      sampledEnergy = Pmin + rand * dp;
      sampledRI = Rindex->Value(sampledEnergy);
      cosTheta = BetaInverse / sampledRI;      

      if (cosTheta > 1) {
        std::cout << "Warning - PrtCherenkovProcess:  cosTheta " << cosTheta << std::endl;
        return pParticleChange;
      }
      
      sin2Theta = (1.0 - cosTheta) * (1.0 + cosTheta);
      rand = G4UniformRand();

      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while (rand * maxSin2 > sin2Theta);

    // Create photon momentum direction vector. The momentum direction is still
    // with respect to the coordinate system where the primary particle
    // direction is aligned with the z axis
    rand = G4UniformRand();
    G4double phi = twopi * rand;
    G4double sinPhi = std::sin(phi);
    G4double cosPhi = std::cos(phi);
    G4double sinTheta = std::sqrt(sin2Theta);
    G4ParticleMomentum photonMomentum(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);

    // Rotate momentum direction back to global reference system
    photonMomentum.rotateUz(p0);

    // Determine polarization of new photon
    G4ThreeVector photonPolarization(cosTheta * cosPhi, cosTheta * sinPhi, -sinTheta);

    // Rotate back to original coord system
    photonPolarization.rotateUz(p0);

    // Generate a new photon:
    auto aCerenkovPhoton = new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), photonMomentum);

    aCerenkovPhoton->SetPolarization(photonPolarization);
    aCerenkovPhoton->SetKineticEnergy(sampledEnergy);

    G4double NumberOfPhotons, N;

    do {
      rand = G4UniformRand();
      NumberOfPhotons = MeanNumberOfPhotons1 - rand * (MeanNumberOfPhotons1 - MeanNumberOfPhotons2);
      N = G4UniformRand() * std::max(MeanNumberOfPhotons1, MeanNumberOfPhotons2);
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while (N > NumberOfPhotons);

    G4double delta = rand * aStep.GetStepLength();
    G4double deltaTime =
      delta / (pPreStepPoint->GetVelocity() +
               rand * (pPostStepPoint->GetVelocity() - pPreStepPoint->GetVelocity()) * 0.5);

    G4double aSecondaryTime = t0 + deltaTime;
    G4ThreeVector aSecondaryPosition = x0 + rand * aStep.GetDeltaPosition();

    // Generate new G4Track object:
    G4Track *aSecondaryTrack = new G4Track(aCerenkovPhoton, aSecondaryTime, aSecondaryPosition);

    aSecondaryTrack->SetTouchableHandle(aStep.GetPreStepPoint()->GetTouchableHandle());
    aSecondaryTrack->SetParentID(aTrack.GetTrackID());
    // aSecondaryTrack->SetCreatorModelID(secID);
    aParticleChange.AddSecondary(aSecondaryTrack);
  }

  if (verboseLevel > 1) {
    G4cout << "\n Exiting from G4Cerenkov::DoIt -- NumberOfSecondaries = "
           << aParticleChange.GetNumberOfSecondaries() << G4endl;
  }

  return pParticleChange;
}

G4double PrtCherenkovProcess::GetAverageNumberOfPhotons(const G4double charge, const G4double beta,
                                                        const G4Material *aMaterial,
                                                        G4MaterialPropertyVector *Rindex) const {

  constexpr G4double Rfact = 369.81 / (eV * cm);
  if (beta <= 0.0) return 0.0;
  G4double BetaInverse = 1. / beta;

  // Vectors used in computation of Cerenkov Angle Integral:
  // 	- Refraction Indices for the current material
  //	- new G4PhysicsFreeVector allocated to hold CAI's
  std::size_t materialIndex = aMaterial->GetIndex();

  // Retrieve the Cerenkov Angle Integrals for this material
  G4PhysicsVector *CerenkovAngleIntegrals = ((*thePhysicsTable)(materialIndex));

  std::size_t length = CerenkovAngleIntegrals->GetVectorLength();
  if (0 == length) return 0.0;

  // Min and Max photon energies
  G4double Pmin = Rindex->Energy(0);
  G4double Pmax = Rindex->GetMaxEnergy();

  // Min and Max Refraction Indices
  G4double nMin = Rindex->GetMinValue();
  G4double nMax = Rindex->GetMaxValue();

  // Max Cerenkov Angle Integral
  G4double CAImax = (*CerenkovAngleIntegrals)[length - 1];

  G4double dp, ge;
  // If n(Pmax) < 1/Beta -- no photons generated
  if (nMax < BetaInverse) {
    dp = 0.0;
    ge = 0.0;
  }
  // otherwise if n(Pmin) >= 1/Beta -- photons generated
  else if (nMin > BetaInverse) {
    dp = Pmax - Pmin;
    ge = CAImax;
  }
  // If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then we need to find a P such
  // that the value of n(P) == 1/Beta. Interpolation is performed by the
  // GetEnergy() and Value() methods of the G4MaterialPropertiesTable and
  // the Value() method of G4PhysicsVector.
  else {
    Pmin = Rindex->GetEnergy(BetaInverse);
    dp = Pmax - Pmin;

    G4double CAImin = CerenkovAngleIntegrals->Value(Pmin);
    ge = CAImax - CAImin;

    if (verboseLevel > 1) {
      G4cout << "CAImin = " << CAImin << G4endl << "ge = " << ge << G4endl;
    }
  }

  // Calculate number of photons
  G4double NumPhotons =
    Rfact * charge / eplus * charge / eplus * (dp - ge * BetaInverse * BetaInverse);

  return NumPhotons;
}
