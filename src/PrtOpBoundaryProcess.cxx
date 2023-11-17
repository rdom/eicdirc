#include "PrtOpBoundaryProcess.h"
#include <G4TouchableHistory.hh>
#include "PrtManager.h"
#include "G4TransportationManager.hh"

PrtOpBoundaryProcess::PrtOpBoundaryProcess() : G4OpBoundaryProcess() {
  fLensId = PrtManager::Instance()->getRun()->getLens();
  fRunType = PrtManager::Instance()->getRun()->getRunType();
  fEvType = PrtManager::Instance()->getRun()->getEv();
  fRadiatorL = PrtManager::Instance()->getRun()->getRadiatorL();
}

G4VParticleChange *PrtOpBoundaryProcess::PostStepDoIt(const G4Track &aTrack, const G4Step &aStep) {
  G4StepPoint *pPreStepPoint = aStep.GetPreStepPoint();
  G4StepPoint *pPostStepPoint = aStep.GetPostStepPoint();
  G4VParticleChange *particleChange = G4OpBoundaryProcess::PostStepDoIt(aTrack, aStep);

  // int parentId = aTrack.GetParentID();
  // std::cout<<"parentId   "<<parentId <<std::endl;
  // if(parentId==1) particleChange->ProposeTrackStatus(fStopAndKill);

  double endofbar = 0.5 * fRadiatorL;

  // ideal focusing
  if (fLensId == 10 && fEvType != 3) {
    G4ThreeVector theGlobalPoint1 = pPostStepPoint->GetPosition();
    G4TouchableHistory *touchable = (G4TouchableHistory *)(pPostStepPoint->GetTouchable());
    G4ThreeVector lpoint = touchable->GetHistory()->GetTransform(1).TransformPoint(theGlobalPoint1);

    if (lpoint.getZ() < endofbar + 0.0001 && lpoint.getZ() > endofbar - 0.0001) {
      G4ThreeVector ww = pPreStepPoint->GetTouchableHandle()
                           ->GetHistory()
                           ->GetTopTransform()
                           .Inverse()
                           .TransformPoint(G4ThreeVector(0, 0, endofbar));

      // in global CS
      double newz = endofbar + 630 + 0.1; // lpoint.getZ()

      if (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() != "wGlue")
        particleChange->ProposeTrackStatus(fStopAndKill);
      else aParticleChange.ProposePosition(ww.getX(), ww.getY(), newz);
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->ComputeSafety(
        G4ThreeVector(ww.getX(), ww.getY(), newz));
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
  }

  // ideal focusing
  if (fLensId == 10 && fEvType == 3) {
    endofbar = endofbar + 630;
    G4ThreeVector gpoint = pPostStepPoint->GetPosition();
    G4TouchableHistory *touchable = (G4TouchableHistory *)(pPostStepPoint->GetTouchable());
    G4ThreeVector lpoint = touchable->GetHistory()->GetTransform(1).TransformPoint(gpoint);
    
    if (gpoint.getZ() < endofbar + 0.1 && gpoint.getZ() > endofbar - 0.1) {

      if(fEvType == 3) endofbar = endofbar - 500;
      G4ThreeVector ww = pPreStepPoint->GetTouchableHandle()
                           ->GetHistory()
                           ->GetTopTransform()
                           .Inverse()
                           .TransformPoint(G4ThreeVector(0, 0, endofbar));

      // in global CS
      double newz = endofbar + 500 + 0.1; // lpoint.getZ()
      // if (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() != "wGlue")
        // particleChange->ProposeTrackStatus(fStopAndKill);
      // else
      aParticleChange.ProposePosition(ww.getX(), ww.getY(),newz); 
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->ComputeSafety(
        G4ThreeVector(ww.getX(), ww.getY(), newz));

      // //cyl lens
      // aParticleChange.ProposePosition(ww.getX(), lpoint.getY(), newz);
      // G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->ComputeSafety(
      //   G4ThreeVector(ww.getX(), lpoint.getY(), newz));
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
  }

  if (fRunType == 1 && pPostStepPoint->GetPosition().z() < pPreStepPoint->GetPosition().z()) {
    if (fEvType != 1) particleChange->ProposeTrackStatus(fStopAndKill);
    if (pPreStepPoint->GetPosition().z() < endofbar)
      particleChange->ProposeTrackStatus(fStopAndKill);
  }  

  if (aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName() == "wExpVol" &&
      pPostStepPoint->GetPosition().z() < pPreStepPoint->GetPosition().z()) {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() == "wLens3" &&
      pPostStepPoint->GetPosition().z() < pPreStepPoint->GetPosition().z()) {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // kill photons outside bar and prizm
  if (GetStatus() == FresnelRefraction &&
      aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName() == "wDirc") {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if ((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() == "wLens1" ||
       aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() == "wLens2") &&
      aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName() == "wDirc") {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // // black edge of the lens3
  // if((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3"
  //     &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc")
  //    || (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3"
  // 	 &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wLens3")){
  //   particleChange->ProposeTrackStatus(fStopAndKill);
  // }

  if (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() == "wLens1" &&
      aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName() == "wLens1") {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }
  if (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() == "wLens2" &&
      aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName() == "wLens2") {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  return particleChange;
}
