#include "PrtEvent.h"

ClassImp(PrtEvent)

// // -----   Default constructor   -------------------------------------------
PrtEvent::PrtEvent(){ 
  fDecoderId = -1;
  fId = -1;
  fType = 0;
  fTime = -1;

  fPhysList = 0;
  fAngle = 0;
  fMomentum = TVector3(0,0,0);
  fPosition = TVector3(0,0,0);
  fHitSize = 0;
  fGeometry = 0;
  fEvType = 0;
  fLens = -1;
  fTrigger = 0;
  fTest = 0;
}

void PrtEvent::AddHit(PrtHit hit){
  fHitArray.push_back(hit);
  fHitSize++;
}

TString PrtEvent::PrintInfo(){
  TString info="Basic sim information: \n";
  info += Form("Physics list %d \n",fPhysList);
  info += Form("Particle  id %d \n",fParticle);
  info += Form("Particle momentum %f \n", fMomentum.Mag());
  info += Form("Geometry id %d \n", fGeometry);
  info += Form("Ev type %d \n", fEvType);
  info += Form("Lens  id %d \n",    fLens);
  return info;
}
