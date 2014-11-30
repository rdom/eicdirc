#include "PrtPrizmHit.h"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<PrtPrizmHit>* PrtPrizmHitAllocator=0;

PrtPrizmHit::PrtPrizmHit()
 : G4VHit(),
   fPrizmID(-1),
   fTrackID(-1),
   fNormalId(-1),
   fEdep(0.),
   fPos(G4ThreeVector())
{}

PrtPrizmHit::~PrtPrizmHit() {}

PrtPrizmHit::PrtPrizmHit(const PrtPrizmHit& right)
  : G4VHit()
{
  fPrizmID =  right.fPrizmID;
  fTrackID   = right.fTrackID;
  fNormalId = right.fNormalId;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
}

const PrtPrizmHit& PrtPrizmHit::operator=(const PrtPrizmHit& right)
{
  fPrizmID =  right.fPrizmID;
  fTrackID   = right.fTrackID;
  fNormalId = right.fNormalId;
  fEdep      = right.fEdep;
  fPos       = right.fPos;

  return *this;
}

G4int PrtPrizmHit::operator==(const PrtPrizmHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

void PrtPrizmHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void PrtPrizmHit::Print()
{
  G4cout
     << "  trackID: " << fTrackID << " chamberNb: " << fNormalId
     << "Edep: "
     << std::setw(7) << G4BestUnit(fEdep,"Energy")
     << " Position: "
     << std::setw(7) << G4BestUnit( fPos,"Length")
     << G4endl;
}

