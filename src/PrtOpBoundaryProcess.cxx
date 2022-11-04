#include "PrtOpBoundaryProcess.h"
#include <G4TouchableHistory.hh>
#include "PrtManager.h"


PrtOpBoundaryProcess::PrtOpBoundaryProcess() : G4OpBoundaryProcess() {
  fLensId = PrtManager::Instance()->getRun()->getLens();
  fRunType = PrtManager::Instance()->getRun()->getRunType();
  fEvType = PrtManager::Instance()->getRun()->getEv();
}

G4VParticleChange *PrtOpBoundaryProcess::PostStepDoIt(const G4Track &aTrack, const G4Step &aStep) {
  G4StepPoint *pPreStepPoint = aStep.GetPreStepPoint();
  G4StepPoint *pPostStepPoint = aStep.GetPostStepPoint();
  G4VParticleChange *particleChange = G4OpBoundaryProcess::PostStepDoIt(aTrack, aStep);

  // int parentId = aTrack.GetParentID();
  // std::cout<<"parentId   "<<parentId <<std::endl;
  // if(parentId==1) particleChange->ProposeTrackStatus(fStopAndKill);

  double endofbar = 0.5 * (4200 + 4 * 0.05); // 1250/2.;

  // ideal focusing
  if (fLensId == 10) {
    G4ThreeVector theGlobalPoint1 = pPostStepPoint->GetPosition();
    G4TouchableHistory *touchable = (G4TouchableHistory *)(pPostStepPoint->GetTouchable());
    G4ThreeVector lpoint = touchable->GetHistory()->GetTransform(1).TransformPoint(theGlobalPoint1);
    if (lpoint.getZ() < endofbar + 0.0001 && lpoint.getZ() > endofbar - 0.0001) {
      G4ThreeVector ww = pPreStepPoint->GetTouchableHandle()
                           ->GetHistory()
                           ->GetTopTransform()
                           .Inverse()
                           .TransformPoint(G4ThreeVector(0, 0, endofbar));

      if (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() != "wGlue")
        particleChange->ProposeTrackStatus(fStopAndKill);
      else aParticleChange.ProposePosition(ww.getX(), ww.getY(), lpoint.getZ() - 0.0005);
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->ComputeSafety(G4ThreeVector(ww.getX(), ww.getY(), lpoint.getZ() - 0.0005));
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
