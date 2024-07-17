#include "PrtEvent.h"

ClassImp(PrtEvent)

  PrtEvent::PrtEvent()
  : fPid(0), fTime(0), fTof(0), fTofPi(0), fTofP(0), fMomentum(TVector3(0, 0, 0)),
    fMomentumBefore(TVector3(0, 0, 0)), fMomentumAfter(TVector3(0, 0, 0)), fPosition(TVector3(0, 0, 0)),
    fPositionAfter(TVector3(0, 0, 0)) {}
