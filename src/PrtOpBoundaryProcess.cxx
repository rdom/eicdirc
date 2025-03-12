#include "PrtOpBoundaryProcess.h"
#include "G4TouchableHistory.hh"
#include "PrtManager.h"
#include "G4TransportationManager.hh"

PrtOpBoundaryProcess::PrtOpBoundaryProcess() : G4OpBoundaryProcess() {
  fLensId = PrtManager::Instance()->getRun()->getLens();
  fRunType = PrtManager::Instance()->getRun()->getRunType();
  fEvType = PrtManager::Instance()->getRun()->getEv();
  fRadiatorL = PrtManager::Instance()->getRun()->getRadiatorL();
  fRadiatorW = PrtManager::Instance()->getRun()->getRadiatorW();
  fRadiatorH = PrtManager::Instance()->getRun()->getRadiatorH();
}

G4VParticleChange *PrtOpBoundaryProcess::PostStepDoIt(const G4Track &aTrack, const G4Step &aStep) {
  G4StepPoint *pPreStepPoint = aStep.GetPreStepPoint();
  G4StepPoint *pPostStepPoint = aStep.GetPostStepPoint();
  G4VParticleChange *particleChange = G4OpBoundaryProcess::PostStepDoIt(aTrack, aStep);

  auto prename = pPreStepPoint->GetPhysicalVolume()->GetName();
  auto posname = pPostStepPoint->GetPhysicalVolume()->GetName();

  // int parentId = aTrack.GetParentID();
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
      double newz = 2685 + 0.1; //endofbar + 630 + 0.1; // lpoint.getZ()

      if (prename != "wGlue") particleChange->ProposeTrackStatus(fStopAndKill);
      else aParticleChange.ProposePosition(ww.getX(), ww.getY(), newz);
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->ComputeSafety(
        G4ThreeVector(ww.getX(), ww.getY(), newz));
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
  }

  // ideal focusing
  if (fLensId == 10 && fEvType == 3) {
    endofbar = 1855.15;
    G4ThreeVector gpoint = pPostStepPoint->GetPosition();
    G4TouchableHistory *touchable = (G4TouchableHistory *)(pPostStepPoint->GetTouchable());
    G4ThreeVector lpoint = touchable->GetHistory()->GetTransform(1).TransformPoint(gpoint);

    if (gpoint.getZ() < endofbar + 0.1 && gpoint.getZ() > endofbar - 10.1) {

      // if(fEvType == 3) endofbar = endofbar - evprismlengh;
      G4ThreeVector ww = pPreStepPoint->GetTouchableHandle()
                           ->GetHistory()
                           ->GetTopTransform()
                           .Inverse()
                           .TransformPoint(G4ThreeVector(0, 0, endofbar));

      // in global CS
      endofbar = 1855.15;
      double newz = endofbar + 0.1; // lpoint.getZ()
      // if (prename != "wGlue")
      // particleChange->ProposeTrackStatus(fStopAndKill);
      // else
      aParticleChange.ProposePosition(ww.getX(), ww.getY(), newz);
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

  if (posname == "wExpVol" &&
      pPostStepPoint->GetPosition().z() < pPreStepPoint->GetPosition().z()) {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if (prename == "wLens3" && pPostStepPoint->GetPosition().z() < pPreStepPoint->GetPosition().z()) {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // kill photons outside bar and prizm
  if (GetStatus() == FresnelRefraction && posname == "wDirc") {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // black edge of the lens1
  if (prename == "wLens1" && posname == "wDirc") {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // black edge of the lens2
  if (prename == "wLens2" && posname == "wDirc") {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // // black edge of the lens3
  // if ((prename == "wLens3" && posname == "wDirc") || (prename == "wLens3" && posname == "wLens3")) {
  //   particleChange->ProposeTrackStatus(fStopAndKill);
  // }

  // no reflections inside lens
  if (prename == "wLens1" && posname == "wLens1") {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }
  if (prename == "wLens2" && posname == "wLens2") {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if(1){  
    if (prename == "wBar" && posname == "wDirc" && GetStatus()  == TotalInternalReflection) {      
      auto gpos = pPostStepPoint->GetPosition();
      auto touchable = (G4TouchableHistory *)(pPreStepPoint->GetTouchable());
      auto lpos = touchable->GetHistory()->GetTopTransform().TransformPoint(gpos);

      double a = PrtManager::Instance()->getRun()->getTest3();
      double as = atan(a / (0.5 * fRadiatorH)); // angle from sagita
      int study = PrtManager::Instance()->getRun()->getStudy();

      auto gmom = G4ThreeVector(*aParticleChange.GetMomentumDirection());
      if (study == 201) { // a
        if (lpos.y() > 0.5 * fRadiatorW - 0.001) gmom.rotateZ(-2 * a);
        if (lpos.y() < -0.5 * fRadiatorW + 0.001) gmom.rotateZ(-2 * a);
      }

      if (study == 202) { // b
        // if (lpos.y() > 0.5 * fRadiatorW - 0.001)
	  gmom.rotateZ(-2 * a);
      }

      if (study == 203) { // c
        if (lpos.y() > 0.5 * fRadiatorW - 0.001) gmom.rotateZ(-2 * a);
        if (lpos.y() < -0.5 * fRadiatorW + 0.001) gmom.rotateZ(2 * a);
      }

      if (study == 204) { // d
	if (lpos.y() > 0.5 * fRadiatorW - 0.001) {
          double h = 0.5 * fRadiatorH;
          double d = lpos.x();
          double r = (h * h + a * a) / (2 * a);
          double b = acos(sqrt(r * r - d * d) / r);
          gmom.rotateZ(2 * b);
          if (lpos.x() > 0) gmom.rotateZ(2 * b);
          else gmom.rotateZ(-2 * b);
        }

        if (lpos.y() < -0.5 * fRadiatorW + 0.001) {
	  double h = 0.5 * fRadiatorH;
          double d = lpos.x();
          double r = (h * h + a * a) / (2 * a);
          double b = acos(sqrt(r * r - d * d) / r);
          gmom.rotateZ(2 * b);
          if (lpos.x() > 0) gmom.rotateZ(-2 * b);
          else gmom.rotateZ(2 * b);
        }
      }

      if (study == 205) { // e
        if (lpos.y() > 0.5 * fRadiatorW - 0.001) {
          double h = 0.5 * fRadiatorH;
          double d = lpos.x();
          double r = (h * h + a * a) / (2 * a);
          double b = acos(sqrt(r * r - d * d) / r);
          gmom.rotateZ(2 * b);
          if (lpos.x() > 0) gmom.rotateZ(2 * b);
          else gmom.rotateZ(-2 * b);
        }
      }

      if (study == 206) { // f
        if (lpos.y() > 0.5 * fRadiatorW - 0.001) {
	  double h = fRadiatorH;
          double d = 0.5 * h + lpos.x();
          double r = (h * h + a * a) / (2 * a);
          double b = acos(sqrt(r * r - d * d) / r);
          gmom.rotateZ(2 * b);
        }
      }

      if (study == 207) { // g
        if (lpos.y() > 0.5 * fRadiatorW - 0.001 || lpos.y() < -0.5 * fRadiatorW + 0.001) {
          double b = G4RandGauss::shoot(0, a);
          gmom.rotateZ(2 * b);
        }
      }

      aParticleChange.ProposeMomentumDirection(gmom);
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);     
    }
  }

  return particleChange;
}
