#include "G4UserTrackingAction.hh"

#include "TGraph.h"
#include "TRandom.h"


class G4Track;
class PrtTrackingAction : public G4UserTrackingAction
{
public:
  PrtTrackingAction();
  void PreUserTrackingAction(const G4Track* aTrack);

private:
  TRandom* fRand;
  TGraph* fDetEff;
};
