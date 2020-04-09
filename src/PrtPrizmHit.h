#ifndef PrtPrizmHit_h
#define PrtPrizmHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


class PrtPrizmHit : public G4VHit
{
  public:
    PrtPrizmHit();
    PrtPrizmHit(const PrtPrizmHit&);
    virtual ~PrtPrizmHit();

    // operators
    const PrtPrizmHit& operator=(const PrtPrizmHit&);
    G4int operator==(const PrtPrizmHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetPrizmID  (G4int id)         { fPrizmID = id; };
    void SetTrackID  (G4int track)      { fTrackID = track; };
    void SetNormalId (G4int chamb)      { fNormalId = chamb; };
    void SetEdep     (G4double de)      { fEdep = de; };
    void SetPos      (G4ThreeVector xyz){ fPos = xyz; };

    // Get methods
    G4int GetPrizmID() const     { return fPrizmID; };
    G4int GetTrackID() const     { return fTrackID; };
    G4int GetNormalId() const    { return fNormalId; };
    G4double GetEdep() const     { return fEdep; };
    G4ThreeVector GetPos() const { return fPos; };

  private:
      G4int         fPrizmID;
      G4int         fTrackID;
      G4int         fNormalId;
      G4double      fEdep;
      G4ThreeVector fPos;
};


typedef G4THitsCollection<PrtPrizmHit> PrtPrizmHitsCollection;

extern G4ThreadLocal G4Allocator<PrtPrizmHit>* PrtPrizmHitAllocator;

inline void* PrtPrizmHit::operator new(size_t)
{
  if(!PrtPrizmHitAllocator)
      PrtPrizmHitAllocator = new G4Allocator<PrtPrizmHit>;
  return (void *) PrtPrizmHitAllocator->MallocSingle();
}

inline void PrtPrizmHit::operator delete(void *hit)
{
  PrtPrizmHitAllocator->FreeSingle((PrtPrizmHit*) hit);
}

#endif
