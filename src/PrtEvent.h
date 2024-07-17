// -----------------------------------------
// PrtEvent.h
//
// author  : r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtEvent_h
#define PrtEvent_h 1

#include "TObject.h"
#include "TString.h"

#include <vector>
#include "PrtHit.h"

class PrtEvent : public TObject {

 public:
  PrtEvent();
  ~PrtEvent(){};

  void addHit(PrtHit hit) { fHitArray.push_back(hit); }

  // accessors
  Int_t getPid() const { return fPid; }
  Double_t getTime() const { return fTime; }
  Double_t getTof() { return fTof; }
  Double_t getTofPi() { return fTofPi; }
  Double_t getTofP() { return fTofP; }
  TVector3 getMomentum() const { return fMomentum; }
  TVector3 getMomentumBefore() { return fMomentumBefore; }
  TVector3 getMomentumAfter() { return fMomentumAfter; }
  TVector3 getPosition() const { return fPosition; }
  TVector3 getPositionAfter() const { return fPositionAfter; }
  std::vector<PrtHit> getHits() { return fHitArray; }

  // mutators
  void setPid(Int_t v) { fPid = v; }
  void setTime(Double_t v) { fTime = v; }
  void setTof(Double_t v) { fTof = v; }
  void setTofPi(Double_t v) { fTofPi = v; }
  void setTofP(Double_t v) { fTofP = v; }
  void setMomentum(TVector3 v) { fMomentum = v; }
  void setMomentumBefore(TVector3 val) { fMomentumBefore = val; }
  void setMomentumAfter(TVector3 val) { fMomentumAfter = val; }
  void setPosition(TVector3 v) { fPosition = v; }
  void setPositionAfter(TVector3 v) { fPositionAfter = v; }

 private:
  Int_t fPid;
  Double_t fTime;
  Double_t fTof;
  Double_t fTofPi;
  Double_t fTofP;
  TVector3 fMomentum;
  TVector3 fMomentumBefore;
  TVector3 fMomentumAfter;
  TVector3 fPosition;
  TVector3 fPositionAfter;
  std::vector<PrtHit> fHitArray;

  ClassDef(PrtEvent, 8);
};
#endif
