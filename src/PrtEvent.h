// -----------------------------------------
// PrtEvent.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtEvent_h
#define PrtEvent_h 1

#include "TObject.h"
#include "TString.h"

#include <vector>
#include "PrtHit.h"

class PrtEvent: public TObject  {

public:

  PrtEvent(); 	//the default constructor
  ~PrtEvent(){}; 

  void AddHit(PrtHit hit);
  PrtHit GetHit(Int_t ind) { return fHitArray[ind]; }
  TString PrintInfo();

  // Accessors 
  Int_t GetDecoderId() const { return fDecoderId; }
  Int_t GetId() const { return fId; }
  Int_t GetType() const { return fType; }
  Double_t GetTime() const { return fTime; }

  Double_t GetAngle()      const { return fAngle; }
  Double_t GetPhi()        const { return fPhi; }
  Int_t GetPhysList()      const { return fPhysList; }
  Int_t GetParticle()      const { return fParticle; }
  TVector3 GetMomentum()   const { return fMomentum; }
  TVector3 GetPosition()   const { return fPosition; }
  Int_t GetHitSize()       const { return fHitSize; }
  Int_t GetGeometry()      const { return fGeometry; }
  Int_t GetLens()          const { return fLens; }
  Int_t GetTrigger()       const { return fTrigger; } 
  Double_t GetTest1()      const { return fTest1; }
  Double_t GetTest2()      const { return fTest2; }
  Double_t GetPrismStepX(){ return fPrismStepX; }
  Double_t GetPrismStepY(){ return fPrismStepY; }
  Double_t GetBeamX(){ return fBeamX; }
  Double_t GetBeamZ(){ return fBeamZ; }
  Double_t GetTimeRes(){ return fTimeRes; }
  TString  GetInfo() { return fInfo; }
  
  // Mutators
  void SetDecoderId(Int_t val)  { fDecoderId=val; }
  void SetId(Int_t val)        { fId=val; }
  void SetType(Int_t val)        { fType=val; }
  void SetTime(Double_t val)      { fTime=val; }

  void SetPhysList(Int_t val) { fPhysList = val; }
  void SetAngle(Double_t val) { fAngle = val; }
  void SetPhi(Double_t val) { fPhi = val; }
  void SetParticle(Int_t val) { fParticle = val; }
  void SetMomentum(TVector3 val) { fMomentum = val; }
  void SetPosition(TVector3 val) { fPosition = val; }
  void SetGeometry(Int_t val) { fGeometry = val; }
  void SetLens(Int_t val) { fLens = val; }
  void SetTrigger(Int_t val) { fTrigger = val; }
  void SetTest1(Double_t val) { fTest1 = val; }
  void SetTest2(Double_t val) { fTest2 = val; }
  void SetPrismStepX(Double_t val){ fPrismStepX = val; }
  void SetPrismStepY(Double_t val){ fPrismStepY = val; }
  void SetBeamX(Double_t val){ fBeamX = val; }
  void SetBeamZ(Double_t val){ fBeamZ = val; }
  void SetTimeRes(Double_t val){ fTimeRes = val; }
  void SetInfo(TString val){ fInfo = val; }

private: 
  Int_t fDecoderId;
  Int_t fId;
  Int_t fType;
  Double_t fTime;


  Int_t fHitSize;
  std::vector<PrtHit> fHitArray;

  Int_t fPhysList;
  Int_t fParticle;
  Double_t fAngle;
  Double_t fPhi;
  TVector3 fMomentum;
  TVector3 fPosition;
  Int_t fGeometry;
  Int_t fLens;
  Int_t fTrigger; 
  Double_t fTest1;
  Double_t fTest2;
  Double_t fPrismStepX;
  Double_t fPrismStepY;
  Double_t fBeamX;
  Double_t fBeamZ;
  Double_t fTimeRes;
  TString  fInfo;
  
  ClassDef(PrtEvent, 5);
};
#endif
