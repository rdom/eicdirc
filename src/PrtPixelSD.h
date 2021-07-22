// -----------------------------------------
// PrtPixelSD.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtPixelSD_h
#define PrtPixelSD_h 1

#include <vector>

#include "G4VSensitiveDetector.hh"

class PrtPixelSD : public G4VSensitiveDetector {
 public:
  PrtPixelSD(const G4String &name, const G4String &hitsCollectionName);
  virtual ~PrtPixelSD();

  // methods from base class
  virtual void Initialize(G4HCofThisEvent *hitCollection);
  virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory *history);
  virtual void EndOfEvent(G4HCofThisEvent *hitCollection);

 private:
  int fEvType, fMcpLayout, fRunType;
  double fRadiatorL, fRadiatorW, fRadiatorH;
  int fMultHit[50][256];
  int fMap_Mpc[50][256];
};

#endif
