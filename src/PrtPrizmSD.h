#ifndef PrtPrizmSD_h
#define PrtPrizmSD_h 1

#include <vector>
#include "G4VSensitiveDetector.hh"

#include "PrtPrizmHit.h"

class G4Step;
class G4HCofThisEvent;


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

#endif
