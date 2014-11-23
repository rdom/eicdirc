#include "G4RichTrajectory.hh"
#include "G4TrackingManager.hh"
#include "G4IdentityTrajectoryFilter.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"

#include "PrtTrackingAction.h"

PrtTrackingAction::PrtTrackingAction()
: G4UserTrackingAction()
{ 
 
}


void PrtTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Require rich trajectory...
  fpTrackingManager->SetTrajectory(new G4RichTrajectory(aTrack));
  // Activate storing of auxiliary points for smoother trajectory...
  static G4IdentityTrajectoryFilter curvedFilter;
  G4TransportationManager::GetTransportationManager()->
    GetPropagatorInField()->SetTrajectoryFilter(&curvedFilter);

}
