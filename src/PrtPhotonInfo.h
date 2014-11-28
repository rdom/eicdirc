// -----------------------------------------
// PrtPhotonInfo.h
//
// Created on: 18.10.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------
// Container for DrcTrackInfo

#ifndef PrtPhotonInfo_h
#define PrtPhotonInfo_h 1

#include "PrtAmbiguityInfo.h"

#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include <vector>

class PrtPhotonInfo : public TObject {

public:
  
  // Default constructor
  PrtPhotonInfo ();

  // Default destructor
  ~PrtPhotonInfo ();

  // Copy constructor 
  PrtPhotonInfo (const PrtPhotonInfo& val):TObject() { *this = val; }  

  // Mutators
  void SetHitTime(Double_t val)               {fHitTime = val;}
  void SetReflected(Bool_t val)               {fReflected = val;}
  void SetEvReflections(Int_t val)            {fEvReflections = val;}
  void SetSensorId(Int_t val)                 {fSensorId = val;}
  void SetMcPrimeMomentumInBar(TVector3 val)  {fMcPrimeMomentumInBar = val;}
  void SetMcCherenkovInBar(Double_t val)      {fMcCherenkovInBar = val;}
  
  void AddAmbiguity(PrtAmbiguityInfo ambiguity);

  // Accessors  
  Double_t    GetHitTime()	              {return fHitTime;}    
  Bool_t GetReflected()	                      {return fReflected;}
  Int_t GetEvReflections()	              {return fEvReflections;}
  Int_t GetSensorId()	                      {return fSensorId;}

  Int_t GetAmbiguitySize()	              {return fAmbiguitySize;}
  PrtAmbiguityInfo GetAmbiguity(Int_t id)     {return fAmbiguityArray[id];}
  TVector3 GetMcPrimeMomentumInBar()	      {return fMcPrimeMomentumInBar;}
  Double_t  GetMcCherenkovInBar()	      {return fMcCherenkovInBar;}    

protected:

  std::vector<PrtAmbiguityInfo> fAmbiguityArray;
  Int_t fAmbiguitySize;

  Int_t    fSensorId;
  Double_t fHitTime;
  Bool_t   fReflected;
  Int_t    fEvReflections;
  TVector3 fMcPrimeMomentumInBar;
  Double_t fMcCherenkovInBar;
 
  ClassDef(PrtPhotonInfo,1)
};

#endif
