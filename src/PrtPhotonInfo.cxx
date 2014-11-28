// -----------------------------------------
// PrtPhotonInfo.h
//
// Created on: 18.10.2013
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#include "PrtPhotonInfo.h"

ClassImp(PrtPhotonInfo)

// -----   Default constructor   -------------------------------------------
PrtPhotonInfo::PrtPhotonInfo(){ 
  fAmbiguitySize = 0;
  fHitTime = 0;
  fReflected = kFALSE;
  fEvReflections = 0;
  fSensorId = 0;
}

PrtPhotonInfo::~PrtPhotonInfo(){ 
}


void PrtPhotonInfo::AddAmbiguity(PrtAmbiguityInfo ambiguity){
  fAmbiguityArray.push_back(ambiguity);
  fAmbiguitySize++;
}
