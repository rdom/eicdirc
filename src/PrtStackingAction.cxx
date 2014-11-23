//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrtStackingAction.h"

#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrtStackingAction::PrtStackingAction()
  : G4UserStackingAction(),
    fScintillationCounter(0), fCerenkovCounter(0)
{
 // quantum efficiency data from Alex Britting, Jan 25, 2011
  // unit is percent
  // first value is at 200 nm, last at 700 nm
  // credible range start around 250nm, >= 280nm to be safe
  double fEfficiencyR[1000], fLambda[1000];
  double fEfficiency[501] = {231.84,615.36,657.4,258.78,9839.92,44.67,67.87,51.01,41.49,5.36,49.4,2.13,35.49,8.66,5.03,7.51,13.27,18.71,3.92,3.66,8.2,0.56,7.68,2.87,10.06,3.47,3.39,6.99,6.01,4.92,6.25,5.97,6.92,8.29,10.45,8.68,8.6,9.79,11.76,9.53,10.98,9.9,10.97,11.31,10.88,10.78,12.16,12.38,12.37,13.04,12.36,13.18,13.7,13.85,13.66,13.98,14.55,14.93,14.82,14.97,14.98,15.14,15.35,15.37,15.43,15.49,15.59,15.84,15.84,15.92,16.01,16.22,16.41,16.42,16.52,16.86,17.1,17.17,17.22,17.46,17.79,17.99,18.13,18.33,18.34,18.53,18.72,18.95,19.02,19.15,19.28,19.45,19.66,19.69,19.77,19.73,19.95,19.98,20.17,20.29,20.33,20.37,20.47,20.48,20.57,20.75,20.8,20.84,20.86,20.88,21.0,21.06,21.0,21.06,21.06,21.04,21.1,21.14,21.08,21.17,21.3,21.38,21.49,21.58,21.69,21.77,21.87,22.02,22.13,22.29,22.35,22.45,22.53,22.55,22.64,22.67,22.73,22.74,22.71,22.79,22.76,22.77,22.76,22.75,22.78,22.7,22.68,22.72,22.66,22.64,22.7,22.67,22.71,22.67,22.75,22.77,22.83,22.84,22.93,22.97,23.0,23.08,23.16,23.27,23.25,23.37,23.44,23.49,23.55,23.52,23.58,23.64,23.63,23.58,23.64,23.63,23.62,23.64,23.63,23.66,23.59,23.59,23.56,23.58,23.63,23.57,23.66,23.62,23.67,23.64,23.54,23.57,23.51,23.53,23.45,23.3,23.41,23.25,23.21,23.08,23.01,22.92,22.9,22.76,22.76,22.61,22.53,22.48,22.39,22.29,22.24,22.2,22.12,22.07,21.96,21.89,21.87,21.76,21.74,21.58,21.49,21.48,21.37,21.29,21.2,21.17,21.03,20.98,20.92,20.85,20.76,20.69,20.58,20.56,20.47,20.37,20.32,20.24,20.13,20.08,19.9,19.84,19.77,19.69,19.63,19.51,19.41,19.27,19.06,19.01,18.87,18.7,18.49,18.41,18.17,17.98,17.84,17.69,17.5,17.25,17.15,16.98,16.79,16.66,16.48,16.32,16.19,16.02,15.88,15.77,15.67,15.5,15.39,15.23,15.09,15.04,14.92,14.75,14.7,14.5,14.45,14.34,14.25,14.16,14.13,14.0,13.92,13.84,13.76,13.73,13.61,13.54,13.52,13.45,13.41,13.39,13.31,13.22,13.17,13.13,13.06,13.2,13.09,12.97,12.92,12.73,12.65,12.4,12.22,12.02,11.79,11.59,11.33,11.03,10.68,10.46,10.14,9.88,9.62,9.36,9.14,8.87,8.63,8.51,8.24,8.07,7.88,7.77,7.65,7.52,7.35,7.27,7.21,7.1,6.92,6.89,6.79,6.74,6.56,6.54,6.5,6.39,6.33,6.25,6.27,6.14,6.06,6.04,6.01,5.91,5.89,5.79,5.75,5.75,5.67,5.61,5.51,5.52,5.43,5.43,5.34,5.31,5.35,5.23,5.2,5.14,5.11,5.11,5.01,4.98,4.93,4.99,4.89,4.82,4.87,4.8,4.7,4.65,4.65,4.61,4.49,4.56,4.44,4.42,4.44,4.35,4.35,4.27,4.29,4.19,4.13,4.08,4.02,4.07,3.92,3.95,3.88,3.82,3.86,3.74,3.71,3.66,3.72,3.62,3.55,3.56,3.57,3.45,3.38,3.36,3.36,3.28,3.25,3.19,3.26,3.13,3.17,3.15,3.04,2.98,2.93,2.98,2.9,2.89,2.9,2.81,2.74,2.81,2.68,2.73,2.7,2.57,2.58,2.55,2.55,2.37,2.39,2.39,2.44,2.37,2.26,2.27,2.27,2.23,2.26,2.14,2.08,2.15,2.06,2.09,2.04,2.0,1.95,2.02,1.87,1.9,1.8,1.87,1.85,1.87,1.81,1.86,1.74,1.74,1.63,1.59,1.5,1.5,1.44,1.47,1.32,1.24,1.28,1.19,1.21,1.21,1.1,1.1,1.05,1.06,0.94,0.92,0.87,0.92,0.81,0.86,0.78,0.77,0.8,0.67,0.7,0.81,0.61,0.64,0.71,0.66,0.67,0.68,0.69,0.68,0.73};
  
  double fCollectionEff=0.65; // Collection Efficiency 
  double credibleLimit=280.;
  
  for(Int_t i=0; i<1000; i++){
    fLambda[i] = i;
  }
  
  for (Int_t i=0;i<1000;i++) {
    if (i<(Int_t)(credibleLimit) || i > 700){
	fEfficiencyR[i]=0;
    }else{
      // total detector efficiency
      fEfficiencyR[i]=fEfficiency[i-200]/100.*fCollectionEff;
    }
  }   
  fRand = new TRandom();
  fDetEff = new TGraph(1000, fLambda,fEfficiencyR); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrtStackingAction::~PrtStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
PrtStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  { // particle is optical photon
    if(aTrack->GetParentID()>0)
    { // particle is secondary
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
        fScintillationCounter++;
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov")
        fCerenkovCounter++;
    }
  }
  
  if(aTrack->GetDefinition()->GetParticleName()=="opticalphoton" && aTrack->GetParentID()!=1)  return fKill;


  G4String ParticleName = aTrack->GetDynamicParticle()->
    GetParticleDefinition()->GetParticleName();

  if(ParticleName == "opticalphoton") {
    // apply detector efficiency at the production stage:    
    if(true){
      Double_t lambda = 197.0*2.0*pi/(aTrack->GetMomentum().mag()*1.0E6);          
      Double_t ra = fRand->Uniform(0., 1.);
      if(ra > fDetEff->Eval(lambda)){ 
	return fKill;
      }
    }
  }

  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrtStackingAction::NewStage()
{
  // G4cout << "Number of Scintillation photons produced in this event : "
  //        << fScintillationCounter << G4endl;
  G4cout << "Number of Cerenkov photons produced in this event : "
         << fCerenkovCounter << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrtStackingAction::PrepareNewEvent()
{
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
