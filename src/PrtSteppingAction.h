#ifndef PrtSteppingAction_h
#define PrtSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

/// Stepping action class

class PrtSteppingAction : public G4UserSteppingAction {
 public:
  PrtSteppingAction();
  virtual ~PrtSteppingAction();

  // method from the base class
  virtual void UserSteppingAction(const G4Step *);

 private:
  G4int fScintillationCounter;
  G4int fCerenkovCounter;
  G4int fEventNumber;
};

#endif
