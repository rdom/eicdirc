// -----------------------------------------
// PrtHit.h
//
// author  : r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtHit_h
#define PrtHit_h 1

#include <vector>

#include "TObject.h"
#include "TVector3.h"

class PrtHit : public TObject {

 public:
  PrtHit();
  ~PrtHit(){};

  // accessors
  Int_t getChannel() { return fChannel; }
  Int_t getPrism() { return fPrism; }
  Int_t getPmt() { return fPmt; }
  Int_t getPixel() { return fPixel; }
  Double_t getLeadTime() { return fLeadTime; }
  Double_t getTotTime() { return fTotTime; }
  Long_t getPathInPrizm() { return fPathInPrizm; }
  TVector3 getPosition() { return fPosition; }
  TVector3 getMomentum() { return fMomentum; }

  // mutators
  void setChannel(Int_t val) { fChannel = val; }
  void setPrism(Int_t val) { fPrism = val; }
  void setPmt(Int_t val) { fPmt = val; }
  void setPixel(Int_t val) { fPixel = val; }
  void setLeadTime(Double_t val) { fLeadTime = val; }
  void setTotTime(Double_t val) { fTotTime = val; }
  void setPathInPrizm(Long_t val) { fPathInPrizm = val; }
  void setPosition(TVector3 val) { fPosition = val; }
  void setMomentum(TVector3 val) { fMomentum = val; }

 protected:
  Int_t fChannel;
  Int_t fPrism;  
  Int_t fPmt;
  Int_t fPixel;
  Double_t fLeadTime;
  Double_t fTotTime;
  Long_t fPathInPrizm;
  TVector3 fPosition;
  TVector3 fMomentum;

  ClassDef(PrtHit, 9)
};

#endif
