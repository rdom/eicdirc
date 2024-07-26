// -----------------------------------------
// PrtRun.h
//
// author  : r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtRun_h
#define PrtRun_h 1

#include <iostream>
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

class PrtRun : public TObject {

 public:
  PrtRun();
  ~PrtRun(){};

  // accessors
  TString getInfo();
  TString getShortInfo() const { return fShortInfo; }
  TString getName() const { return fName; }
  Int_t getId() const { return fId; }
  Int_t getRunType() const { return fRunType; }
  Int_t getStudy() const { return fStudy; }
  Bool_t getMc() const { return fMc; }
  Double_t getTheta() const { return fTheta; }
  Double_t getPhi() const { return fPhi; }
  Int_t getPhysList() const { return fPhysList; }
  Int_t getPid() const { return fPid; }
  Double_t getMomentum() const { return fMomentum; }
  Int_t getField() const { return fField; }
  Int_t getGeometry() const { return fGeometry; }
  Int_t getRadiator() const { return fRadiator; }
  Int_t getPmtLayout() const { return fPmtLayout; }
  Int_t getEv() const { return fEv; }
  Int_t getLens() const { return fLens; }
  Int_t getTrigger() const { return fTrigger; }
  Int_t getNpmt() const { return fNpmt; }
  Int_t getNpix() const { return fNpix; }
  Int_t getCorrection() const { return fCorrection; }
  Double_t getPrismStepX() const { return fPrismStepX; }
  Double_t getPrismStepY() const { return fPrismStepY; }
  Double_t getBeamX() const { return fBeamX; }
  Double_t getBeamZ() const { return fBeamZ; }
  Double_t getBeamSize() const { return fBeamSize; }
  Double_t getTrackingResTheta() const { return fTrackingResTheta; }     	
  Double_t getTrackingResPhi() const { return fTrackingResPhi; }
  Double_t getTimeSigma() const { return fTimeSigma; }
  Double_t getTimeCut() const { return fTimeCut; }
  Double_t getSimOffset() const { return fSimOffset; }
  Double_t getRadiatorL() const { return fRadiatorL; }
  Double_t getRadiatorW() const { return fRadiatorW; }
  Double_t getRadiatorH() const { return fRadiatorH; }
  Double_t getDarkNoise() const { return fDarkNoise; }
  Double_t getTest1() const { return fTest1; }
  Double_t getTest2() const { return fTest2; }
  Double_t getTest3() const { return fTest3; }

  // mutators
  void setInfo(TString v) { fShortInfo = v; }
  void setName(TString v) { fName = v; }
  void setId(Int_t v) { fId = v; }
  void setRunType(Int_t v) { fRunType = v; }
  void setStudy(Int_t v) { fStudy = v; }
  void setMc(Bool_t v) { fMc = v; }
  void setTheta(Double_t v) { fTheta = v; }
  void setPhi(Double_t v) { fPhi = v; }
  void setPhysList(Int_t v) { fPhysList = v; }
  void setPid(Int_t v) { fPid = v; }
  void setMomentum(Double_t v) { fMomentum = v; }
  void setField(Int_t v) { fField = v; }
  void setGeometry(Int_t v) { fGeometry = v; }
  void setRadiator(Int_t v) { fRadiator = v; }
  void setPmtLayout(Int_t v);
  void setEv(Int_t v) { fEv = v; }
  void setLens(Int_t v) { fLens = v; }
  void setTrigger(Int_t v) { fTrigger = v; }
  void setNpmt(Int_t v) { fNpmt = v; }
  void setNpix(Int_t v) { fNpix = v; }  
  void setCorrection(Int_t v) { fCorrection = v; }
  void setPrismStepX(Double_t v) { fPrismStepX = v; }
  void setPrismStepY(Double_t v) { fPrismStepY = v; }
  void setBeamX(Double_t v) { fBeamX = v; }
  void setBeamZ(Double_t v) { fBeamZ = v; }
  void setBeamSize(Double_t v) { fBeamSize = v; }
  void setTrackingResTheta(Double_t v) { fTrackingResTheta = v; }
  void setTrackingResPhi(Double_t v) { fTrackingResPhi = v; }
  void setTimeSigma(Double_t v) { fTimeSigma = v; }
  void setTimeCut(Double_t v) { fTimeCut = v; }
  void setSimOffset(Double_t v) { fSimOffset = v; }
  void setRadiatorL(Double_t v) { fRadiatorL = v; }
  void setRadiatorW(Double_t v) { fRadiatorW = v; }
  void setRadiatorH(Double_t v) { fRadiatorH = v; }
  void setDarkNoise(Double_t v) { fDarkNoise = v; }
  void setTest1(Double_t v) { fTest1 = v; }
  void setTest2(Double_t v) { fTest2 = v; }
  void setTest3(Double_t v) { fTest3 = v; }

 private:
  TString fShortInfo;
  TString fName;
  Int_t fId;
  Int_t fRunType;
  Int_t fStudy;
  Bool_t fMc;
  Int_t fPhysList;
  Int_t fPid;
  Int_t fField;
  Int_t fGeometry;
  Int_t fEv;
  Int_t fLens;
  Int_t fRadiator;
  Int_t fPmtLayout;
  Int_t fTrigger;
  Int_t fNpmt;
  Int_t fNpix;
  Int_t fCorrection;

  Double_t fTheta;
  Double_t fPhi;
  Double_t fMomentum;

  Double_t fPrismStepX;
  Double_t fPrismStepY;
  Double_t fBeamX;
  Double_t fBeamZ;
  Double_t fBeamSize;
  Double_t fTrackingResTheta;
  Double_t fTrackingResPhi;
  Double_t fTimeSigma;
  Double_t fTimeCut;
  Double_t fSimOffset;

  Double_t fRadiatorL, fRadiatorW, fRadiatorH;
  Double_t fDarkNoise;
  
  Double_t fTest1;
  Double_t fTest2;
  Double_t fTest3;

  ClassDef(PrtRun, 3);
};
#endif
