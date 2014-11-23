#ifndef PrtPrizmSD_h
#define PrtPrizmSD_h 1

#include <vector>
#include "G4VSensitiveDetector.hh"

#include "PrtPrizmHit.h"

class G4Step;
class G4HCofThisEvent;

/// Calorimeter sensitive detector class
///
/// In Initialize(), it creates one hit for each calorimeter layer and one more
/// hit for accounting the total quantities in all layers.
///
/// The values are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class PrtPrizmSD : public G4VSensitiveDetector
{
public:
  PrtPrizmSD(const G4String& name, 
	     const G4String& hitsCollectionName, 
	     G4int nofCells);
  virtual ~PrtPrizmSD();
  
  // methods from base class
  virtual void   Initialize(G4HCofThisEvent* hitCollection);
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
  virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

private: 
  PrtPrizmHitsCollection* fHitsCollection;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
