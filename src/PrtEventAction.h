// -----------------------------------------
// PrtEventAction.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtEventAction_h
#define PrtEventAction_h 1

#include "G4UserEventAction.hh"

#include "TGraph.h"
#include "TRandom.h"

class G4Track;
class PrtEventAction : public G4UserEventAction
{
public:
  PrtEventAction() {;}
  virtual ~PrtEventAction() {;}
  // void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
protected:
  G4EventManager* fpEventManager;
};
#endif
