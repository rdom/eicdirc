#ifndef PrtHit_h
#define PrtHit_h 1

#include <vector>

#include "TObject.h"
#include "TVector3.h"

class PrtHit : public TObject {

public:   
 
  //Constructor
  PrtHit();

  ~PrtHit(){};
 
  // Accessors 
  Int_t GetMcpId()       { return fMcpId; }
  Int_t GetPixelId()     { return fPixelId; }
  Int_t GetType()        { return fType; }
  Int_t GetParticleId()  { return fParticleId; }   
  Int_t GetParentParticleId()  { return fParentParticleId; }  
  Int_t GetNreflectionsInPrizm()  { return fNreflectionsInPrizm; }
  Double_t GetPathInPrizm()  { return fPathInPrizm; }
  TVector3 GetLocalPos()     { return fLocalPos; }
  TVector3 GetGlobalPos()     { return fGlobalPos; }
  TVector3 GetDigiPos()     { return fDigiPos; }
  TVector3 GetMomentum()     { return fMomentum; }
  TVector3 GetPosition()     { return fPosition; }
  
  Int_t GetChannel() { return fChannel;}
  Int_t GetTdc() { return fTdc;}
  Int_t GetMultiplicity() { return fMultiplicity; }
  Int_t GetLeadTime(Int_t ind) { return fLeadTime[ind]; } 
  Int_t GetTrailTime(Int_t ind) { return fTrailTime[ind]; } 

  // Mutators
  void SetMcpId(Int_t val)   { fMcpId = val; }
  void SetPixelId(Int_t val) { fPixelId = val; }
  void SetType(Int_t val)    { fType = val; }
  void SetParticleId(Int_t val)  { fParticleId = val; }   
  void SetParentParticleId(Int_t val)  { fParentParticleId = val; }  
  void SetNreflectionsInPrizm(Int_t val)  { fNreflectionsInPrizm = val; }
  void SetPathInPrizm(Double_t val) { fPathInPrizm = val; }
  void SetLocalPos(TVector3 val)   { fLocalPos = val; }
  void SetGlobalPos(TVector3 val)  { fGlobalPos = val; }
  void SetDigiPos(TVector3 val)    { fDigiPos = val; }
  void SetMomentum(TVector3 val)    { fMomentum = val; }
  void SetPosition(TVector3 val)    { fPosition = val; }

  void SetChannel(Int_t val) { fChannel=val; }
  void SetTdc(Int_t val) { fTdc = val; }
  void SetMultiplicity(Int_t val) { fMultiplicity = val; }
 
  void SetLeadTime(Int_t ind, Int_t val)  { fLeadTime[ind]=val; } 
  void SetTrailTime(Int_t ind, Int_t val) { fTrailTime[ind]=val; } 

protected:

  Int_t fMcpId;
  Int_t fPixelId;
  Int_t fType;
  Int_t fParticleId; 
  Int_t fParentParticleId;
  Int_t fNreflectionsInPrizm;
  Double_t fPathInPrizm;
  TVector3 fLocalPos;
  TVector3 fGlobalPos;
  TVector3 fDigiPos;
  TVector3 fMomentum;
  TVector3 fPosition;

  Int_t fChannel;
  Int_t fTdc;
  Int_t fMultiplicity;
  Int_t fLeadTime[4];    
  Int_t fTrailTime[4];  

  ClassDef(PrtHit,1)
};

#endif
