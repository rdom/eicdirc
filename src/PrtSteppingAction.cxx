#include "PrtSteppingAction.h"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"

#include <G4TouchableHistory.hh>
#include "PrtManager.h"

PrtSteppingAction::PrtSteppingAction() : G4UserSteppingAction() {
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
  fEventNumber = -1;
}

PrtSteppingAction::~PrtSteppingAction() {}

void PrtSteppingAction::UserSteppingAction(const G4Step *step) {

  G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  G4Track *track = step->GetTrack();
  
  if (track->GetCurrentStepNumber() > 50000 || track->GetTrackLength() > 30000) {
    // std::cout<<"WRN: too many steps or track length > 30 m  N="
    // <<track->GetCurrentStepNumber()<<" Len = "<<track->GetTrackLength()/1000. <<std::endl;
    track->SetTrackStatus(fStopAndKill);
  }

  int parentId = track->GetParentID();
  G4String prevname = "", postvname = "";
  if (step->GetPreStepPoint()->GetPhysicalVolume() &&
      step->GetPostStepPoint()->GetPhysicalVolume()) {
    prevname = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
    postvname = step->GetPostStepPoint()->GetPhysicalVolume()->GetName();
  }
  
  if (prevname == "wMcp" && postvname == "wDirc") track->SetTrackStatus(fStopAndKill);
  if (prevname == "wFd" && postvname == "wDirc") track->SetTrackStatus(fStopAndKill);
  if (prevname == "wBBWindow" && postvname == "wDirc") track->SetTrackStatus(fStopAndKill);
  if (prevname == "wMcp" && postvname == "wPixel") track->SetTrackStatus(fStopButAlive);

}
