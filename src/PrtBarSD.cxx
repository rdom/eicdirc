#include "PrtBarSD.h"
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

G4Mutex PrtBarSD::fMutex = G4MUTEX_INITIALIZER;

PrtBarSD::PrtBarSD(const G4String &name, const G4String &hitsCollectionName, G4int nofCells)
  : G4VSensitiveDetector(name), fHitsCollection(NULL) {

  G4AutoLock tuberier(&fMutex);
  collectionName.insert(hitsCollectionName);
}

PrtBarSD::~PrtBarSD() {}

void PrtBarSD::Initialize(G4HCofThisEvent *hce) {

  // Create hits collection
  fHitsCollection = new PrtBarHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

G4bool PrtBarSD::ProcessHits(G4Step *step, G4TouchableHistory *hist) {

  // energy deposit
  G4Track *track = step->GetTrack();
  int parentId = track->GetParentID();
  if (parentId > 0) return true; // only primaries

  G4String ParticleName = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  if (ParticleName == "opticalphoton") return true;

  PrtBarHit *newHit = new PrtBarHit();
  newHit->SetTrackID(step->GetTrack()->GetTrackID());
  newHit->SetPos(step->GetPostStepPoint()->GetPosition());
  newHit->SetMom(track->GetMomentum());

  // store normal to the closest boundary
  // G4Navigator *theNavigator =
  //   G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  // double fP = track->GetDynamicParticle()->GetTotalMomentum();
  // double fEnergy = track->GetDynamicParticle()->GetTotalEnergy();
  // double cherenkov = acos(1 / (1.47125 * (fP / fEnergy)));
  // PrtManager::Instance()->SetCurrentCherenkov(cherenkov);
  // PrtManager::Instance()->Event()->SetTime();

  auto pstep = step->GetPostStepPoint();
  G4ThreeVector gpos = pstep->GetPosition();
  G4ThreeVector gmom = pstep->GetMomentum();
  G4TouchableHistory *touchable = (G4TouchableHistory *)(pstep->GetTouchable());
  G4ThreeVector lpos = touchable->GetHistory()->GetTransform(1).TransformPoint(gpos);

  if (fHitsCollection->entries() == 0) {
    PrtManager::Instance()->getEvent()->setMomentum(TVector3(gmom.x(), gmom.y(), gmom.z()));
    PrtManager::Instance()->getEvent()->setPosition(TVector3(lpos.x(), lpos.y(), lpos.z()));
  } else {    
    PrtManager::Instance()->getEvent()->setMomentumAfter(TVector3(gmom.x(), gmom.y(), gmom.z()));
    PrtManager::Instance()->getEvent()->setPositionAfter(TVector3(lpos.x(), lpos.y(), lpos.z()));
  }

  fHitsCollection->insert(newHit);

  return true;
}

void PrtBarSD::EndOfEvent(G4HCofThisEvent *) {

  if (verboseLevel > 1) {
    G4int nofHits = fHitsCollection->entries();
    G4cout << "\n-------->Bar Hits Collection: in this event they are " << nofHits
           << " hits in the tracker chambers: " << G4endl;
    for (G4int i = 0; i < nofHits; i++) (*fHitsCollection)[i]->Print();
  }
}
