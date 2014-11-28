#include "PrtHit.h"

ClassImp(PrtHit)

// -----   Default constructor   -------------------------------------------
PrtHit::PrtHit(){ 
  fMcpId=-1;
  fPixelId=-1;
  fType=-1;
  fLocalPos = TVector3(0,0,0);
  fGlobalPos = TVector3(0,0,0);
  fDigiPos = TVector3(0,0,0);
  fMomentum = TVector3(0,0,0);
  fCherenkovMC = 0;
  fChannel= -1;
  fTdc = -1;
  fMultiplicity = 0;
  fLeadTime = -1; 
  fTotTime = -1; 
}
