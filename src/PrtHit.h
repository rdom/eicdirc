// -----------------------------------------
// PrtHit.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

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
  Double_t GetCherenkovMC()  { return fCherenkovMC;}
  
  Int_t GetMcpId()       { return fMcpId; }
  Int_t GetPixelId()     { return fPixelId; }
  Int_t GetChannel() { return fChannel;}
  Int_t GetTdc() { return fTdc;}
  Int_t GetTrb() { return fTrb;}
  Int_t GetMultiplicity() { return fMultiplicity; }
  Double_t GetLeadTime() { return fLeadTime; } 
  Double_t GetTotTime() { return fTotTime; } 
    
  // Mutators
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
  void SetCherenkovMC(Double_t val)    { fCherenkovMC = val; }

  
  void SetMcpId(Int_t val)   { fMcpId = val; }
  void SetPixelId(Int_t val) { fPixelId = val; }
  void SetChannel(Int_t val) { fChannel=val; }
  void SetTdc(Int_t val) { fTdc = val; }
  void SetTrb(Int_t val) { fTrb = val; }
  void SetMultiplicity(Int_t val) { fMultiplicity = val; }
  void SetLeadTime(Double_t val) { fLeadTime=val; } 
  void SetTotTime(Double_t val) { fTotTime=val; } 

protected:

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
  Double_t fCherenkovMC;
  
  Int_t fMcpId;
  Int_t fPixelId;
  Int_t fChannel;
  Int_t fTdc;
  Int_t fTrb;
  Int_t fMultiplicity;
  Double_t fLeadTime;    
  Double_t fTotTime;  

  ClassDef(PrtHit,4)
};

#endif
