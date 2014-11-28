// -----------------------------------------
// PrtTrackInfo.h
//
// Created on: 18.10.2013
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#include "PrtTrackInfo.h"

ClassImp(PrtTrackInfo)

// -----   Default constructor   -------------------------------------------
PrtTrackInfo::PrtTrackInfo(){ 
  fPhotonSize      = 0;
  fMcPdg           = 0;
  fMcMomentum      = TVector3(0,0,0);
  fMcMomentumInBar = TVector3(0,0,0);
  fMcCherenkov     = 0;
  fPdg             = 0;
  fMomentum        = TVector3(0,0,0);
  fCherenkov       = 0;
  fAngle           = 0;
  fInfo            = "";
}


PrtTrackInfo::~PrtTrackInfo(){ 
}

void PrtTrackInfo::AddPhoton(PrtPhotonInfo photon){
  fPhotonArray.push_back(photon);
  fPhotonSize++;
}
