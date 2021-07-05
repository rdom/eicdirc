// -----------------------------------------
// PrtBarSD.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtBarSD_h
#define PrtBarSD_h 1

#include <vector>
#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "PrtBarHit.h"

class G4Step;
class G4HCofThisEvent;

class PrtBarSD : public G4VSensitiveDetector {
 public:
  PrtBarSD(const G4String &name, const G4String &hitsCollectionName, G4int nofCells);
  virtual ~PrtBarSD();

  // methods from base class
  virtual void Initialize(G4HCofThisEvent *hitCollection);
  virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory *history);
  virtual void EndOfEvent(G4HCofThisEvent *hitCollection);

 private:
  PrtBarHitsCollection *fHitsCollection;
  static G4Mutex fMutex;
};

#endif
