// -----------------------------------------
// PrtLutNode.h
//
// Created on: 09.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------
// Container for look-up table

#ifndef PrtLutNode_h
#define PrtLutNode_h 1

#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include <vector>
#include <iostream>

class PrtLutNode : public TObject {

public:    
  
  // Default constructor
  PrtLutNode ();

  // Standard constructors
  PrtLutNode (Int_t detectorId);

  // Copy constructor 
  // PrtLutNode (PrtLutNode& node) { *this = node; }  

  // Modifiers
  void AddEntry(Int_t nodeId, TVector3 dir, Long_t pathid, Int_t nrefl, Double_t time, TVector3 hitpos, TVector3 digipos);
  void SetDigiPos(TVector3 pos){fDigiPos = pos;}

  // Accessors
  Int_t Entries() { return fSize; }
  Int_t GetDetectorId() { return fDetectorId; }

  TVector3 GetEntry(Int_t entry) { return fNodeArray[entry]; }
  Long_t GetPathId(Int_t entry){ return fPathIdArray[entry]; }
  Int_t GetNRefl(Int_t entry){ return fNRefl[entry]; }
  Double_t GetTime(Int_t entry){ return fTimeArray[entry]; }
  TVector3 GetHitPos(Int_t entry){ return fHitPos[entry]; }
  TVector3 GetDigiPos(){ return fDigiPos; }

private:

  Int_t fDetectorId;
  Int_t fSize;
  TVector3 fDigiPos;

  std::vector<TVector3> fNodeArray;
  std::vector<TVector3> fHitPos;
  std::vector<Long_t> fPathIdArray;
  std::vector<Int_t> fNRefl;
  std::vector<Double_t> fTimeArray;

  ClassDef(PrtLutNode, 1);
  
};

#endif
