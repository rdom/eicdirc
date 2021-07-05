#ifndef PrtDetectorConstruction_h
#define PrtDetectorConstruction_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"


#include "TChain.h"

#include "PrtRun.h"
#include "PrtDetectorConstructionMessenger.h"

class PrtDetectorConstructionMessenger;

class PrtDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  PrtDetectorConstruction();
  virtual ~PrtDetectorConstruction();

 public:
  virtual G4VPhysicalVolume *Construct();
  virtual void ConstructSDandField();
  void DefineMaterials();
  void SetVisualization();
  void SetRotation(G4double angle);
  void DrawHitBox(G4int id);
  void SetLens(G4int id);
  void SetQuantumEfficiency(G4int id);

 private:
  PrtRun *fRun;

  G4LogicalVolume *lExpHall;
  G4LogicalVolume *lFront;
  G4LogicalVolume *lDirc;
  G4LogicalVolume *lFd;
  G4LogicalVolume *lBar;
  G4LogicalVolume *lGlue;
  G4LogicalVolume *lMirror;
  G4LogicalVolume *lLens1;
  G4LogicalVolume *lLens2;
  G4LogicalVolume *lLens3;
  G4LogicalVolume *lPrizm;
  G4LogicalVolume *lPrizmT1;
  G4LogicalVolume *lPrizmT2;
  G4LogicalVolume *lWedge;
  G4LogicalVolume *lSWedge;
  G4LogicalVolume *lBlock;
  G4LogicalVolume *lFmirror;
  G4LogicalVolume *lWindow;
  G4LogicalVolume *lMcp;
  G4LogicalVolume *lPixel;
  G4LogicalVolume *lExpVol;
  G4LogicalVolume *lGlueE;
  G4VPhysicalVolume *pDirc[100];

  G4VPhysicalVolume *wBar;
  G4VPhysicalVolume *wGlue;
  G4VPhysicalVolume *wMirror;
  G4VPhysicalVolume *wDirc;

  G4Material *defaultMaterial; // material for bars
  G4Material *BarMaterial;     // material for bars
  G4Material *OilMaterial;
  G4Material *MirrorMaterial; // material of mirror
  G4Material *epotekMaterial;
  G4Material *Nlak33aMaterial;
  G4Material *PbF2Material;
  G4Material *SapphireMaterial;
  G4Material *frontMaterial;

  int fNRow;
  int fNCol;
  int fNBoxes;
  double fRadius;
  double fNpix1;
  double fNpix2;
  double fBoxWidth;
  int fGeomType;
  int fEvType;
  int fMcpLayout;
  int fLensId;
  double fdTilt;
  double fNBar;
  double fHall[3];
  double fBar[3];
  double fMirror[3];
  double fFd[3];
  double fPrizm[4];
  double fPrizmT[6];
  double fLens[4];
  double fMcpTotal[3];
  double fMcpActive[3];
  double fBarsGap;
  double fRotAngle;
  double *fQuantumEfficiency;
  int fRunType, fTest1, fTest2, fTest3;

  G4ThreeVector fPrismShift;

  G4RotationMatrix *fPrtRot;
  PrtDetectorConstructionMessenger *fGeomMessenger;
};

#endif
