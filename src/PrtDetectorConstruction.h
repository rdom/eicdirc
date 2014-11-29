#ifndef PrtDetectorConstruction_h
#define PrtDetectorConstruction_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"

#include "G4PVPlacement.hh"
#include "PrtDetectorConstructionMessenger.h"

class PrtDetectorConstructionMessenger;

class PrtDetectorConstruction : public G4VUserDetectorConstruction
{ 
public:
  PrtDetectorConstruction();
  virtual ~PrtDetectorConstruction();

public:
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  void DefineMaterials();
  void SetVisualization();
  void SetRotation(G4double angle);
  void DrawHitBox(G4int id);
  void SetLens(G4int id);
  void SetQuantumEfficiency(G4int id);
  

private:
  G4LogicalVolume* lExpHall;
  G4LogicalVolume* lFront;
  G4LogicalVolume* lDirc;
  G4LogicalVolume* lBar;
  G4LogicalVolume* lMirror;
  G4LogicalVolume* lLens1;
  G4LogicalVolume* lLens2;
  G4LogicalVolume* lPrizm;
  G4LogicalVolume* lMcp;
  G4LogicalVolume* lPixel;
  G4PVPlacement*   pPixel[1000];

  G4VPhysicalVolume* wBar;
  G4VPhysicalVolume* wMirror;
  G4VPhysicalVolume* wDirc;

  G4Material*        defaultMaterial; // material for bars
  G4Material*        BarMaterial; // material for bars
  G4Material*        OilMaterial;
  G4Material*        MirrorMaterial; // material of mirror
  G4Material*        epotekMaterial;  
  G4Material*        Nlak33aMaterial;  
  G4Material*        frontMaterial;
  
  G4int fNRow;
  G4int fNCol;
  G4int fGeomId;
  G4int fMcpLayout;
  G4double fNBar;
  G4double fHall[3];
  G4double fBar[3];
  G4double fMirror[3];
  G4double fPrizm[4];
  G4double fLens[4];
  G4double fMcpTotal[3];
  G4double fMcpActive[3];
  G4ThreeVector fPrismShift;


  G4double fRotAngle;
  G4RotationMatrix *fPrtRot;
  PrtDetectorConstructionMessenger* fGeomMessenger;
  G4double *fQuantumEfficiency;
};

#endif
