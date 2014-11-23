#include "PrtSteppingAction.h"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"

PrtSteppingAction::PrtSteppingAction()
: G4UserSteppingAction()
{ 
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fEventNumber = -1;
}

PrtSteppingAction::~PrtSteppingAction(){ }


void PrtSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4int eventNumber = G4RunManager::GetRunManager()->
                                              GetCurrentEvent()->GetEventID();

  if (eventNumber != fEventNumber) {
     // G4cout << " Number of Scintillation Photons in previous event: "
     //        << fScintillationCounter << G4endl;
     // G4cout << " Number of Cerenkov Photons in previous event: "
     //        << fCerenkovCounter << G4endl;
     fEventNumber = eventNumber;
     fScintillationCounter = 0;
     fCerenkovCounter = 0;
  }

  G4Track* track = step->GetTrack();
  
  int parentId = track->GetParentID();

  // std::cout<<"parentId   "<<parentId <<std::endl;
 

  //  G4cout<<step->GetPreStepPoint()->GetPhysicalVolume()->GetName()  <<" - "
  // 	<<step->GetPostStepPoint()->GetPhysicalVolume()->GetName() <<"  "<< G4endl;
  //if(step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Bar" && step->GetPostStepPoint()->GetPhysicalVolume()->GetName()=="ExpHall" ) track->SetTrackStatus(fStopAndKill);
  //if(step->GetPreStepPoint()->GetPosition().x()>10 ) track->SetTrackStatus(fStopAndKill);

  G4String ParticleName = track->GetDynamicParticle()->
                                 GetParticleDefinition()->GetParticleName();

 
  
  //std::cout<<"ParticleName "<<ParticleName <<std::endl;
  
  if (ParticleName == "opticalphoton") return;

  const std::vector<const G4Track*>* secondaries =
                                            step->GetSecondaryInCurrentStep();

  if (secondaries->size()>0) {
     for(unsigned int i=0; i<secondaries->size(); ++i) {
        if (secondaries->at(i)->GetParentID()>0) {
           if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
               == G4OpticalPhoton::OpticalPhotonDefinition()){
              if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
               == "Scintillation")fScintillationCounter++;
              if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
               == "Cerenkov")fCerenkovCounter++;
           }
        }
     }
  }
}

