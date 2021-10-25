// -----------------------------------------
// PrtPrizmSD.cxx
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#include "PrtPrizmSD.h"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"
#include <TVector3.h>

#include "PrtEvent.h"

#include "PrtRunAction.h"
#include "PrtManager.h"

PrtPrizmSD::PrtPrizmSD(const G4String &name, const G4String &hitsCollectionName)
  : G4VSensitiveDetector(name), fHitsCollection(NULL) {

  collectionName.insert(hitsCollectionName);
  fRun = PrtManager::Instance()->getRun();
}

PrtPrizmSD::~PrtPrizmSD() {}

void PrtPrizmSD::Initialize(G4HCofThisEvent *hce) {

  // Create hits collection
  fHitsCollection = new PrtPrizmHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

G4bool PrtPrizmSD::ProcessHits(G4Step *aStep, G4TouchableHistory *hist) {

  G4StepPoint *pPostStepPoint = aStep->GetPostStepPoint();
  if (pPostStepPoint->GetStepStatus() != fGeomBoundary) { // not at boundary
    return true;
  }

  PrtPrizmHit *newHit = new PrtPrizmHit();
  G4TouchableHistory *touchable = (G4TouchableHistory *)(aStep->GetPreStepPoint()->GetTouchable());
  G4String vname = touchable->GetVolume()->GetName();

  if (vname == "wDirc") {
    newHit->SetPrizmID(touchable->GetReplicaNumber());
  }

  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetPos(pPostStepPoint->GetPosition());

  // normal to the closest boundary
  G4Navigator *theNavigator =
    G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  int nid = 0;
  G4bool valid;
  G4ThreeVector normal = theNavigator->GetLocalExitNormal(&valid);
  G4ThreeVector gnormal = theNavigator->GetLocalToGlobalTransform().TransformAxis(-normal);
  normal = touchable->GetHistory()->GetTransform(1).TransformAxis(gnormal); // in lDirc

  if (valid) {
    if (fRun->getEv() == 1) {
      if (vname == "wWedge") {
        if (normal.y() > 0.99) nid = 1;                 // right
        if (normal.y() < -0.99) nid = 2;                // left
        if (normal.x() > 0.99) nid = 3;                 // bottom
        if (fabs(normal.x() + 0.866025) < 0.1) nid = 4; // top
      } else if (vname == "wSWedge" || vname == "wBlock") {
        if (normal.y() > 0.99) nid = 1;                                          // right
        if (normal.y() < -0.99) nid = 2;                                         // left
        if (normal.x() > 0.99) nid = 5;                                          // bottom
        if (fabs(normal.x() + 0.866025) < 0.1) nid = 6;                          // top
        if (normal.z() > 0.99) nid = 7;                                          // front
        if (normal.z() < -0.5 && normal.z() > -0.999 && normal.x() > 0) nid = 8; // mirror
        if (normal.z() < -0.5 && normal.z() > -0.999 && normal.x() < 0) nid = 9; // pd
      } else if (normal.z() < -0.999) nid = 0;
    } else {
      if (vname == "wPrizm") {
        if (normal.y() > 0.99) nid = 1;  // right
        if (normal.y() < -0.99) nid = 2; // left
        if (normal.x() > 0.99) nid = 3;  // bottom
        if (fabs(normal.x() + 0.866025) < 0.1) nid = 4;
      } else if (vname == "wLens3") {
        if (normal.y() > 0.99) nid = 5;  // right
        if (normal.y() < -0.99) nid = 6; // left
        if (normal.x() > 0.99) nid = 7;  // bottom
        if (normal.x() < -0.99) nid = 8; // top
      } else if (vname == "wLens2") {
        nid = 9;
      } else if (fRun->getRunType() == 5 && vname == "wLens1") {
	nid = -5;
	newHit->SetPos(aStep->GetPreStepPoint()->GetMomentum());
	newHit->SetEdep(aStep->GetPreStepPoint()->GetLocalTime());
      }
    }
  }

  // if (fRun->getRunType() == 5 && aStep->GetPreStepPoint()->GetPosition().z() < 2100.100001) {
  //   nid = -5;
  //   newHit->SetPos(aStep->GetPostStepPoint()->GetMomentum());
  //   newHit->SetEdep(aStep->GetPreStepPoint()->GetLocalTime());
  // }

  newHit->SetNormalId(nid);

  fHitsCollection->insert(newHit);

  return true;
}

void PrtPrizmSD::EndOfEvent(G4HCofThisEvent *) {
  if (verboseLevel > 1) {
    G4int nofHits = fHitsCollection->entries();
    G4cout << "\n-------->Prizm Hits Collection: in this event they are " << nofHits
           << " hits in the tracker chambers: " << G4endl;
    for (G4int i = 0; i < nofHits; i++) (*fHitsCollection)[i]->Print();
  }
}
