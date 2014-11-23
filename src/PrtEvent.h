#ifndef PrtEvent_h
#define PrtEvent_h 1

#include "TObject.h"
#include "TClonesArray.h"

#include <vector>
#include "PrtHit.h"

class PrtEvent: public TObject  {

protected: 
  Int_t fSize;
  Int_t fDecoding;
  Int_t fId;
  Int_t fSeqNr;
  Int_t fDate;
  Int_t fTime;
  Int_t fYear;
  Int_t fMonth;
  Int_t fDay;
  Int_t fHour;
  Int_t fMinute;
  Int_t fSecond;
  Int_t fPad;
  Int_t fDataSize; 
  Int_t fPaddedSize; 

  Int_t fMaxMultiplicity;
  Int_t fMaxChannel;
  Int_t fSubEvtId;
  Int_t fErrors;
  Int_t fReferenceChannel;
  Int_t fReferenceTime;
  Int_t fHitSize;
  std::vector<PrtHit> fHitArray;

  Int_t fPhysList;
  Int_t fParticle;
  Double_t fAngle;
  Double_t fMomentum;
  Int_t fGeometry;
  Int_t fLens;

public:

  PrtEvent(); 	//the default constructor

  void AddHit(PrtHit hit);
  PrtHit GetHit(Int_t ind) { return fHitArray[ind]; }

  // Accessors 
  Int_t GetSize() const { return fSize; }
  Int_t GetDecoding() const { return fDecoding; }
  Int_t GetId() const { return fId; }
  Int_t GetDate() const { return fDate; }
  Int_t GetTime() const { return fTime; }
  Int_t GetYear() const { return fYear; }
  Int_t GetMonth() const { return fMonth; }
  Int_t GetDay() const { return fDay; }
  Int_t GetHour() const { return fHour; }
  Int_t GetMinute() const { return fMinute; }
  Int_t GetSecond() const { return fSecond; }
  Int_t GetPad() const { return fPad; }
  Int_t GetDataSize() const { return fDataSize; }
  Int_t GetPaddedSize() const { return fPaddedSize; }

  Int_t GetMaxMultiplicity()  const { return  fMaxMultiplicity; }
  Int_t GetMaxChannel()       const { return fMaxChannel; }
  Int_t GetSubEvtId()         const { return fSubEvtId; }
  Int_t GetErrors() const { return fErrors; }
  Int_t GetReferenceChannel() const { return fReferenceChannel; }
  Int_t GetReferenceTime()    const { return fReferenceTime; }

  Double_t GetAngle()         const { return fAngle; }
  Int_t GetPhysList()      const { return fPhysList; }
  Int_t GetParticle()      const { return fParticle; }
  Double_t GetMomentum()      const { return fMomentum; }
  Int_t GetHitSize()       const { return fHitSize; }
  Int_t GetGeometry()      const { return fGeometry; }
  Int_t GetLens()      const { return fLens; }

  // Mutators
  void SetSize(Int_t val)      { fSize=val; }
  void SetDecoding(Int_t val)  { fDecoding=val; }
  void SetId(Int_t val)        { fId=val; }
  void SetDate(Int_t val)      { fDate=val; }
  void SetTime(Int_t val)      { fTime=val; }
  void SetYear(Int_t val)      { fYear=val; }
  void SetMonth(Int_t val)     { fMonth=val; }
  void SetDay(Int_t val)       { fDay=val; }
  void SetHour(Int_t val)      { fHour=val; }
  void SetMinute(Int_t val)    { fMinute=val; }
  void SetSecond(Int_t val)    { fSecond=val; }
  void SetPad(Int_t val)       { fPad=val; }
  void SetDataSize(Int_t val)  { fDataSize=val; }
  void SetPaddedSize(Int_t val){ fPaddedSize=val; }

  void SetMaxMultInt(Int_t val) { fMaxMultiplicity = val; }
  void SetMaxChannelNr(Int_t val) { fMaxChannel = val; }
  void SetSubEvtId(Int_t val)  { fSubEvtId = val; }
  void SetErrors(Int_t val)    { fErrors = val; }
  void SetReferenceChannel(Int_t val) { fReferenceChannel = val; }

  void SetPhysList(Int_t val) { fPhysList = val; }
  void SetAngle(Double_t val) { fAngle = val; }
  void SetParticle(Int_t val) { fParticle = val; }
  void SetMomentum(Double_t val) { fMomentum = val; }
  void SetReferenceTime(Int_t val) { fReferenceTime = val; }
  void SetGeometry(Int_t val) { fGeometry = val; }
  void SetLens(Int_t val) { fLens = val; }

protected: 
  ClassDef(PrtEvent, 1);
};
#endif
