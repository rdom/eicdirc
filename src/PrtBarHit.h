// -----------------------------------------
// PrtBarHit.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtBarHit_h
#define PrtBarHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"


class PrtBarHit : public G4VHit
{
  public:
    PrtBarHit();
    PrtBarHit(const PrtBarHit&);
    virtual ~PrtBarHit();

    // operators
    const PrtBarHit& operator=(const PrtBarHit&);
    G4int operator==(const PrtBarHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID  (G4int track)      { fTrackID = track; };
    void SetEdep     (G4double de)      { fEdep = de; };
    void SetPos      (G4ThreeVector xyz){ fPos = xyz; };
    void SetMom (G4ThreeVector mom){ fMom = mom; };

    // Get methods
    G4int GetTrackID() const     { return fTrackID; };
    G4double GetEdep() const     { return fEdep; };
    G4ThreeVector GetPos() const { return fPos; };
    G4ThreeVector GetMom() const { return fMom; };

  private:

      G4int         fTrackID;
      G4double      fEdep;
      G4ThreeVector fPos;
      G4ThreeVector fMom;
};


typedef G4THitsCollection<PrtBarHit> PrtBarHitsCollection;
extern G4ThreadLocal G4Allocator<PrtBarHit>* PrtBarHitAllocator;

inline void* PrtBarHit::operator new(size_t)
{
  if(!PrtBarHitAllocator)
      PrtBarHitAllocator = new G4Allocator<PrtBarHit>;
  return (void *) PrtBarHitAllocator->MallocSingle();
}

inline void PrtBarHit::operator delete(void *hit)
{
  PrtBarHitAllocator->FreeSingle((PrtBarHit*) hit);
}

#endif
