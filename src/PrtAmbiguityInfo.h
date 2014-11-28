// -----------------------------------------
// PndDrcLutInfo.h
//
// Created on: 18.10.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtAmbiguityInfo_h
#define PrtAmbiguityInfo_h 1

#include "PrtAmbiguityInfo.h"

#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include <vector>

class PrtAmbiguityInfo : public TObject {

public:    
  
  // Default constructor
  PrtAmbiguityInfo ();
  
  ~PrtAmbiguityInfo () { };
  
  // Copy constructor 
  PrtAmbiguityInfo (const PrtAmbiguityInfo& val):TObject() { *this = val; } 

  // Mutators
  void SetCherenkov(Double_t val)           {fCherenkov = val;}
  void SetBarTime(Double_t val)             {fBarTime = val;}
  void SetEvTime(Double_t val)              {fEvTime = val;}

  // Accessors  
  Double_t GetCherencov() 	            {return fCherenkov;}
  Double_t GetBarTime()	                    {return fBarTime;}
  Double_t GetEvTime()	                    {return fEvTime;}
  
protected:

  Double_t    fCherenkov;
  Double_t    fBarTime;
  Double_t    fEvTime;
    
  ClassDef(PrtAmbiguityInfo,1)
};

#endif
