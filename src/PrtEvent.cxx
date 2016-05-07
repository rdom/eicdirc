#include "PrtEvent.h"

ClassImp(PrtEvent)

// // -----   Default constructor   -------------------------------------------
PrtEvent::PrtEvent(): fDecoderId(-1), fId(-1),fType(0),fTime(-1),fPhysList(0),fAngle(0),
  fMomentum(TVector3(0,0,0)),fPosition(TVector3(0,0,0)),fHitSize(0),fGeometry(0),
  fLens(-1),fTrigger(0),fTest1(0),fTest2(0),fPrismStepX(0),fPrismStepY(0),fBeamX(0),fBeamZ(0),fTimeRes(0),fInfo("") { 
}

void PrtEvent::AddHit(PrtHit hit){
  fHitArray.push_back(hit);
  fHitSize++;
}

TString PrtEvent::PrintInfo(){
  TString info="Basic sim information: \n";
  info += fInfo + "\n";
  info += Form("Physics list %d \n",fPhysList);
  info += Form("Particle  id %d \n",fParticle);
  info += Form("Particle momentum %f \n", fMomentum.Mag());
  info += Form("Geometry id %d \n", fGeometry);
  info += Form("Lens  id %d \n",    fLens);
  info += Form("StemY %f \n",    fPrismStepX);
  info += Form("StemY %f \n",    fPrismStepY);
  info += Form("MCP's time resolution %f \n",    fTimeRes);
  return info;
}
