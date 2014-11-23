#include "PrtEvent.h"

ClassImp(PrtEvent)

// // -----   Default constructor   -------------------------------------------
PrtEvent::PrtEvent(){ 
  fSize = -1;
  fDecoding = -1;
  fId = -1;
  fSeqNr = -1;
  fDate = -1;
  fTime = -1;
  fYear = -1;
  fMonth = -1;
  fDay = -1;
  fHour = -1;
  fMinute = -1;
  fSecond = -1;
  fPad = -1;
  fDataSize = -1; 
  fPaddedSize = -1; 

  fMaxMultiplicity = 0;
  fMaxChannel = 0;
  fSubEvtId = -1;
  fErrors = -1;
  fReferenceChannel = -1; 
  fReferenceTime = -1;

  fPhysList = 0;
  fAngle = 0;
  fMomentum = 0;
  fHitSize = 0;
  fGeometry = 0;
  fLens = -1;
}

void PrtEvent::AddHit(PrtHit hit){
  fHitArray.push_back(hit);
  fHitSize++;
}

