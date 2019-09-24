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

PrtBarSD::PrtBarSD(const G4String& name, 
		   const G4String& hitsCollectionName,
		   G4int nofCells)
  : G4VSensitiveDetector(name), fHitsCollection(NULL){

  G4AutoLock tuberier(&fMutex);
  collectionName.insert(hitsCollectionName);
}

PrtBarSD::~PrtBarSD(){

}

void PrtBarSD::Initialize(G4HCofThisEvent* hce){
  
  // Create hits collection
  fHitsCollection = new PrtBarHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
}

G4bool PrtBarSD::ProcessHits(G4Step* aStep, G4TouchableHistory* hist){

  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4Track *track = aStep->GetTrack();
  G4String ParticleName = track->GetDynamicParticle()->
    GetParticleDefinition()->GetParticleName();
  if (ParticleName == "opticalphoton") return true;
  
  // if (edep==0.) return false;

  PrtBarHit* newHit = new PrtBarHit();
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  newHit->SetMom(track->GetMomentum());
  // // store normal to the closest boundary
  G4Navigator* theNavigator 
    = G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking();
  double fP = track->GetDynamicParticle()->GetTotalMomentum();
  
  double fEnergy = track->GetDynamicParticle()->GetTotalEnergy();
  double cherenkov = acos(1/(1.47125*(fP/fEnergy)));
  PrtManager::Instance()->SetCurrentCherenkov(cherenkov);

  fHitsCollection->insert(newHit);

  // newHit->Print();

  return true;
}

void PrtBarSD::EndOfEvent(G4HCofThisEvent*){

  if ( verboseLevel>1 ) { 
    G4int nofHits = fHitsCollection->entries();
    G4cout << "\n-------->Bar Hits Collection: in this event they are " << nofHits 
	   << " hits in the tracker chambers: " << G4endl;
    for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

