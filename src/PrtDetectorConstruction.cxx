
#include "PrtDetectorConstruction.h"

#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4AutoDelete.hh"
#include "G4ProductionCuts.hh"


#include "PrtManager.h"
#include "PrtBarSD.h"
#include "PrtPrizmSD.h"
#include "PrtPixelSD.h"
#include "PrtField.h"
#include "PrtFastSimModelTracker.h"


PrtDetectorConstruction::PrtDetectorConstruction() : G4VUserDetectorConstruction() {

  fRun = PrtManager::Instance()->getRun();
  fRunType = fRun->getRunType();
  fGeomType = fRun->getGeometry();
  fEvType = fRun->getEv();
  fMcpLayout = fRun->getPmtLayout();
  fLensId = fRun->getLens();
  fNBar = fRun->getRadiator();
  fTest1 = fRun->getTest1();
  fTest2 = fRun->getTest2();
  fTest3 = fRun->getTest3();

  fNRow = 6;
  fNCol = 4;
  if (fNBar < 0) fNBar = 1;

  fHall[0] = 1500;
  fHall[1] = 1500;
  fHall[2] = 3500;
  
  fBar[0] = 17; fBar[1] = 32; fBar[2] = 1050;

  // fPrizm[0] = 170; fPrizm[1] = 300; fPrizm[2] = 30+300*tan(37*deg); fPrizm[3] = 30;
  // fPrizm[0] = 390; fPrizm[1] = 300; fPrizm[3] = 50; fPrizm[2]= fPrizm[3]+300*tan(32*deg);

  fPrizm[0] = 350;
  fPrizm[1] = 300;
  fPrizm[3] = 50;
  fPrizm[2] = fPrizm[3] + 300 * tan(32 * deg);
  fBarsGap = 0.15;

  fdTilt = 80 * deg;
  fPrizmT[0] = 390;
  fPrizmT[1] = 400 - 290 * cos(fdTilt); //
  fPrizmT[2] = 290 * sin(fdTilt);       // hight
  fPrizmT[3] = 50;                      // face
  fPrizmT[4] = 290;
  fPrizmT[5] = 290 * cos(fdTilt);

  fBBWindow[0] = 30;
  fBBWindow[1] = fPrizm[0];
  fBBWindow[2] = 0; // fTest3;

  fMcpTotal[0] = fMcpTotal[1] = 53 + 4;
  fMcpTotal[2] = 1;
  fMcpActive[0] = fMcpActive[1] = 53;
  fMcpActive[2] = 1;
  fLens[0] = fLens[1] = fPrizm[3];
  fLens[2] = 12;
  fRadius = 970;
  fNBoxes = 16;
  
  if (fGeomType == 0 || fGeomType == 10) { // ATHENA
    fNBoxes = 16;
    fRadius = 972.8; // middle of the barbox at 90 degree
    fNBar = 10;
    fBar[2] = 1100;
  }

  if (fGeomType == 1 || fGeomType == 11) { // ePIC ECCE
    fNBoxes = 12;
    // fRadius = 700 + 0.5 * fBar[0]; // old = 729.6;
    fRadius = 770.5;
    fNBar = 10;
    fBar[2] = 1225; // BaBar bars
    fBar[1] = 35;
    fPrizm[0] = 350 + 1.35;
  }

  if (fGeomType == 2 || fGeomType == 12) { // CORE
    fNBoxes = 16;
    fRadius = 501.7;
    fNBar = 5;
    fBar[2] = 700;

    fNRow = 3;
    fNCol = 4;
    fPrizm[0] = 175;

    // fNBoxes = 12; // alternative
    // fNBar = 7;
    // fPrizm[0] = 245;
  }

  fBoxWidth = fPrizm[0];
  std::cout << "fBoxWidth/Prism  " << fBoxWidth<<" x "<< fPrizm[2] << std::endl;

  // fBar[1] = (fPrizm[0] - (fNBar - 1) * fBarsGap) / fNBar;
  // std::cout << "N bars " << fNBar << " bar width " << fBar[1] << std::endl;

  fMirror[0] = fBar[0] + 1;
  fMirror[1] = fPrizm[0];
  fMirror[2] = 1;

  if (fEvType == 1) { // BaBar
    fBar[0] = 17.25;
    fBar[1] = 35;
    fBar[2] = 1050;
    fNBar = 12;
    fRadius = 900;
    fLensId = 0;
    fBoxWidth = fNBar * (fBar[1] + fBarsGap);
    fNRow = 7;
    fNBoxes = 12;
  }

  if (fGeomType < 10) fNBoxes = 1;

  fFd[0] = fBoxWidth;
  fFd[1] = fPrizm[2];
  fFd[2] = 1;

  if (fEvType == 4) {
    fFd[1] = fPrizmT[4];
  }
  if (fEvType == 3) {
    fLens[0] =  fBar[0];
    fLens[2] = 9;
  }
  
  if (fLensId == 0 || fLensId == 10) {
    fLens[2] = 0;
  }
  if (fLensId == 2) {
    fLens[1] = 175;
    fLens[2] = 14.4;
  }
  if (fLensId == 3) {
    fLens[1] = fPrizm[0] / fNBar;
    if (fNBar == 1) fLens[1] = fPrizm[0] / 11.;
  }

  if (fLensId == 6) {
    fLens[1] = fPrizm[0];
  }


  if (fMcpLayout == 4) {
    fMcpTotal[0] = fMcpTotal[1] = 120;
    fMcpActive[0] = fMcpActive[1] = 108;
    fNRow = 3;
    fNCol = 2;
  }

  fRun->setRadiatorW(fBar[1]);
  fRun->setRadiatorH(fBar[0]);

  fPrtRot = new G4RotationMatrix();
  // create a messenger for this class
  fGeomMessenger = new PrtDetectorConstructionMessenger(this);
}

PrtDetectorConstruction::~PrtDetectorConstruction() {}

G4VPhysicalVolume *PrtDetectorConstruction::Construct() {
  DefineMaterials();

  // ------------- Volumes --------------

  // The experimental Hall
  G4Box *gExpHall = new G4Box("gExpHall", fHall[0], fHall[1], fHall[2]);
  lExpHall = new G4LogicalVolume(gExpHall, defaultMaterial, "lExpHall", 0, 0, 0);
  G4VPhysicalVolume *wExpHall =
    new G4PVPlacement(0, G4ThreeVector(), lExpHall, "wExpHall", 0, false, 0);

  double gluethickness = 0.05;
  double dirclength = fBar[2] * 4 + gluethickness * 4;
  fRun->setRadiatorL(dirclength);

  // The DIRC bar box
  G4Box *gDirc = new G4Box("gDirc", 130, 180, 0.5 * dirclength + 500);
  lDirc = new G4LogicalVolume(gDirc, defaultMaterial, "lDirc", 0, 0, 0);
  G4Box *gFd = new G4Box("gFd", 0.5 * fFd[1], 0.5 * fFd[0], 0.5 * fFd[2]);
  lFd = new G4LogicalVolume(gFd, defaultMaterial, "lFd", 0, 0, 0);

  double tphi, dphi = 360 * deg / (double)fNBoxes;

  double center_shift = 630; // makes end at -182
  double rshift = -100; // shift x-center of the bar box to avoid overlaps
  if (fRunType == 1) {
    // LUT
    new G4PVPlacement(0, G4ThreeVector(rshift, 0, center_shift), lDirc, "wDirc", lExpHall, false, 0);
  } else {
    for (int i = 0; i < fNBoxes; i++) {
      tphi = dphi * i;
      
      double dx = (fRadius - rshift) * cos(tphi);
      double dy = (fRadius - rshift) * sin(tphi);

      G4RotationMatrix *tRot = new G4RotationMatrix();
      tRot->rotateZ(-tphi);
      new G4PVPlacement(tRot, G4ThreeVector(dx, dy, center_shift), lDirc, "wDirc", lExpHall, false,
                        i);
    }
  }

  // The Bar
  G4Box *gBar = new G4Box("gBar", 0.5 * fBar[0], 0.5 * fBar[1], 0.5 * fBar[2]);
  lBar = new G4LogicalVolume(gBar, BarMaterial, "lBar", 0, 0, 0);
  double evprismhight = fBar[0]; //fPrizm[3]; // fBar[0]
  double evprismlengh = 893; // Bar[2] // fTest1

  G4Box *gExpVol = new G4Box("gExpVol", 0.5 * evprismhight, 0.5 * fBoxWidth, 0.5 * evprismlengh); 
  lExpVol = new G4LogicalVolume(gExpVol, BarMaterial, "lExpVol", 0, 0, 0);
  
  // Glue
  G4Box *gGlue = new G4Box("gGlue", 0.5 * fBar[0], 0.5 * fBar[1], 0.5 * gluethickness);
  lGlue = new G4LogicalVolume(gGlue, epotekMaterial, "lGlue", 0, 0, 0);
  G4Box *gGlueE = new G4Box("gGlueE", 0.5 * evprismhight, 0.5 * fBoxWidth, 0.5 * gluethickness);
  lGlueE = new G4LogicalVolume(gGlueE, epotekMaterial, "lGlueE", 0, 0, 0);

  // Window
  if (fBBWindow[2] > 0.1) {
    G4Box *gBBWindow = new G4Box("gBBWindow", 0.5 * fBBWindow[0], 0.5 * fBBWindow[1], 0.5 * fBBWindow[2]);
    lBBWindow = new G4LogicalVolume(gBBWindow, BarMaterial, "lBBWindow", 0, 0, 0);    
    new G4PVPlacement(0, G4ThreeVector(rshift, 0, 0.5 * (dirclength + fBBWindow[2])), lBBWindow,
		      "wBBWindow", lDirc, false, 0);
  }
    
  // Tracker
  G4Box *gTracker = new G4Box("gTracker", 0.5, fNBar * 0.5 * fBar[1], 4 * 0.5 * fBar[2]);
  lTracker = new G4LogicalVolume(gTracker, defaultMaterial, "lTracker", 0, 0, 0);

  if (fNBar == 1) {
    for (int j = 0; j < 4; j++) {
      double z = -0.5 * dirclength + 0.5 * fBar[2] + (fBar[2] + gluethickness) * j;
      wBar = new G4PVPlacement(0, G4ThreeVector(rshift, 0, z), lBar, "wBar", lDirc, false, 0);
      wGlue = new G4PVPlacement(0, G4ThreeVector(rshift, 0, z + 0.5 * (fBar[2] + gluethickness)),
                                lGlue, "wGlue", lDirc, false, 0);
    }
    wTracker = new G4PVPlacement(0, G4ThreeVector(rshift - 0.5 * fBar[0] - 1, 0, 0), lTracker,
                                 "wTracker", lDirc, false, 0);
  } else {
    int id = 0, nparts = 4;
    if (fEvType == 3 || fEvType == 5) {
      dirclength = fBar[2] * 3 + evprismlengh + gluethickness * 4;
      fRun->setRadiatorL(dirclength - 2 * evprismlengh);
      double sh = 0;
      if (fEvType == 3) sh = fLens[2];
      nparts = 3;
      double z = -0.5 * dirclength + 0.5 * evprismlengh + (fBar[2] + gluethickness) * 3;
      new G4PVPlacement(0, G4ThreeVector(rshift, 0, z + sh), lExpVol, "wExpVol", lDirc, false, id);
      new G4PVPlacement(0, G4ThreeVector(rshift, 0, z + 0.5 * (evprismlengh + gluethickness) + sh),
                        lGlueE, "wGlue", lDirc, false, id);
    }

    for (int i = 0; i < fNBar; i++) {
      double shifty = i * (fBar[1] + fBarsGap) - 0.5 * fBoxWidth + fBar[1] / 2.;
      for (int j = 0; j < nparts; j++) {
        double z = -0.5 * dirclength + 0.5 * fBar[2] + (fBar[2] + gluethickness) * j;
	new G4PVPlacement(0, G4ThreeVector(rshift, shifty, z), lBar, "wBar", lDirc, false, id);
        wGlue = new G4PVPlacement(0, G4ThreeVector(rshift, shifty, z + 0.5 * (fBar[2] + gluethickness)),
                                  lGlue, "wGlue", lDirc, false, id);
        id++;
      }
    }
    wTracker = new G4PVPlacement(0, G4ThreeVector(rshift - 0.5 * fBar[0] - 1, 0, 0), lTracker,
                                 "wTracker", lDirc, false, 0);
  }

  // The Mirror
  G4Box *gMirror = new G4Box("gMirror", fMirror[0] / 2., fMirror[1] / 2., fMirror[2] / 2.);
  lMirror = new G4LogicalVolume(gMirror, MirrorMaterial, "lMirror", 0, 0, 0);
  wMirror = new G4PVPlacement(0, G4ThreeVector(rshift, 0, -0.5 * dirclength - fMirror[2] / 2.),
                              lMirror, "wMirror", lDirc, false, 0);

  // The Lenses
  if (fLensId == 2) { // 2-layer lens
    double lensrad = 70;
    double lensMinThikness = 2;
    G4Box *gfbox = new G4Box("Fbox", fLens[0] / 2., fLens[1] / 2., fLens[2] / 2.);
    G4ThreeVector zTrans(0, 0, -lensrad + fLens[2] / 2. - lensMinThikness);

    G4Sphere *gsphere = new G4Sphere("Sphere", 0, 70, 0. * deg, 360. * deg, 0. * deg, 380. * deg);
    G4IntersectionSolid *gLens1 =
      new G4IntersectionSolid("Fbox*Sphere", gfbox, gsphere, new G4RotationMatrix(), zTrans);
    G4SubtractionSolid *gLens2 =
      new G4SubtractionSolid("Fbox-Sphere", gfbox, gsphere, new G4RotationMatrix(), zTrans);

    lLens1 = new G4LogicalVolume(gLens1, Nlak33aMaterial, "lLens1", 0, 0, 0); // Nlak33aMaterial
    lLens2 = new G4LogicalVolume(gLens2, BarMaterial, "lLens2", 0, 0, 0);
  }

  if (fLensId == 3) { // 3-component spherical lens
    double lensMinThikness = 2;

    double r1 = 0;// = fTest1;
    double r2 = 0;// = fTest2;

    r1 = (r1 == 0) ? 47.8 : r1;
    r2 = (r2 == 0) ? 29.1 : r2;
    r1 = 62;
    r2 = 36;

    if (fEvType == 3) {
      r1 = 150;
      r2 = 90;
    }

    double thight = fBar[0];
    if (thight < 12) thight = 12;
    double cr2 = sqrt(fLens[1] * fLens[1] / 4. + thight * thight / 4.);
    if (cr2 > r2) {
      std::cout << "bad lens" << std::endl;
      cr2 = r2;
    }
    std::cout << "cr2 " << cr2 << std::endl;
    
    // fLens[2]
    double optimallt= (2 * lensMinThikness + r2 - sqrt(r2 * r2 - cr2 * cr2) + lensMinThikness);
    std::cout << "Used lens thickness = " << fLens[2] << " optimal = "<< optimallt << std::endl;

    G4ThreeVector zTrans1(0, 0,
                          -r1 - fLens[2] / 2. + r1 - sqrt(r1 * r1 - cr2 * cr2) + lensMinThikness);
    G4ThreeVector zTrans2(
      0, 0, -r2 - fLens[2] / 2. + r2 - sqrt(r2 * r2 - cr2 * cr2) + lensMinThikness * 2);

    G4Box *gfbox = new G4Box("Fbox", fLens[0] / 2., fLens[1] / 2., fLens[2] / 2.);
    G4Tubs *gfstube = new G4Tubs("ftube", 0, cr2, fLens[2] / 2., 0, 360 * deg);
    // G4Tubs *gfstube = new G4Tubs("ftube", 0, cr2, 2*2, 0, 360 * deg);

    G4Sphere *gsphere1 = new G4Sphere("Sphere1", 0, r1, 0, 360 * deg, 0, 360 * deg);
    G4Sphere *gsphere2 = new G4Sphere("Sphere2", 0, r2, 0, 360 * deg, 0, 360 * deg);

    G4IntersectionSolid *gbbox = new G4IntersectionSolid(
      "bbox", gfbox, gfbox, new G4RotationMatrix(), G4ThreeVector(0, 0, lensMinThikness * 2));
    G4IntersectionSolid *gsbox = new G4IntersectionSolid(
      "sbox", gfstube, gfbox, new G4RotationMatrix(), G4ThreeVector(0, 0, -lensMinThikness * 2));

    G4UnionSolid *gubox =
      new G4UnionSolid("unionbox", gbbox, gsbox, new G4RotationMatrix(), G4ThreeVector(0, 0, 0));

    G4IntersectionSolid *gLens1 =
      new G4IntersectionSolid("Lens1", gubox, gsphere1, new G4RotationMatrix(), zTrans1);
    G4SubtractionSolid *gLenst =
      new G4SubtractionSolid("temp", gubox, gsphere1, new G4RotationMatrix(), zTrans1);

    G4IntersectionSolid *gLens2 =
      new G4IntersectionSolid("Lens2", gLenst, gsphere2, new G4RotationMatrix(), zTrans2);
    G4SubtractionSolid *gLens3 =
      new G4SubtractionSolid("Lens3", gbbox, gsphere2, new G4RotationMatrix(), zTrans2);

    lLens1 = new G4LogicalVolume(gLens1, BarMaterial, "lLens1", 0, 0, 0);
    // Nlak33aMaterial //PbF2Material //SapphireMaterial
    lLens2 = new G4LogicalVolume(gLens2, SapphireMaterial, "lLens2", 0, 0, 0);
    lLens3 = new G4LogicalVolume(gLens3, BarMaterial, "lLens3", 0, 0, 0);
  }
 
  if (fLensId == 6) { // 3-component cylindrical lens
    double lensMinThikness = 2.0;

    double r1 = 0; // PrtManager::Instance()->GetTest1();
    double r2 = 0; // PrtManager::Instance()->GetTest2();

    lensMinThikness = 2;
    double layer12 = lensMinThikness * 2;

    // r1 = (r1==0)? 27.45: r1;
    // r2 = (r2==0)? 20.02: r2;

    r1 = (r1 == 0) ? 33 : r1;
    r2 = (r2 == 0) ? 24 : r2;

    r1 = 62;
    r2 = 36;
    
    double shight = fBar[0];//25;

    if (fEvType == 3) {
      // fLens[0] = fBar[0];
      fLens[2] = 9;
      r1 = 150;
      r2 = 90;
    }

    G4ThreeVector zTrans1(
      0, 0, -r1 - fLens[2] / 2. + r1 - sqrt(r1 * r1 - shight / 2. * shight / 2.) + lensMinThikness);
    G4ThreeVector zTrans2(
      0, 0, -r2 - fLens[2] / 2. + r2 - sqrt(r2 * r2 - shight / 2. * shight / 2.) + layer12);

    G4Box *gfbox = new G4Box("fbox", 0.5 * fLens[0], 0.5 * fLens[1], 0.5 * fLens[2]);
    G4Box *gcbox = new G4Box("cbox", 0.5 * fLens[0], 0.5 * fLens[1] + 1, 0.5 * fLens[2]);
    G4ThreeVector tTrans1(0.5 * (fLens[0] + shight), 0, -fLens[2] + layer12);
    G4ThreeVector tTrans0(-0.5 * (fLens[0] + shight), 0, -fLens[2] + layer12);
    G4SubtractionSolid *tubox =
      new G4SubtractionSolid("tubox", gfbox, gcbox, new G4RotationMatrix(), tTrans1);
    G4SubtractionSolid *gubox =
      new G4SubtractionSolid("gubox", tubox, gcbox, new G4RotationMatrix(), tTrans0);

    G4Tubs *gcylinder1 = new G4Tubs("Cylinder1", 0, r1, 0.5 * fLens[1], 0 * deg, 360 * deg);
    G4Tubs *gcylinder2 = new G4Tubs("Cylinder2", 0, r2, 0.5 * fLens[1] - 0.5, 0 * deg, 360 * deg);
    G4Tubs *gcylinder1c = new G4Tubs("Cylinder1c", 0, r1, 0.5 * fLens[1] + 0.5, 0 * deg, 360 * deg);
    G4Tubs *gcylinder2c = new G4Tubs("Cylinder2c", 0, r2, 0.5 * fLens[1] + 0.5, 0 * deg, 360 * deg);
    G4RotationMatrix *xRot = new G4RotationMatrix();
    xRot->rotateX(M_PI / 2. * rad);

    G4IntersectionSolid *gLens1 =
      new G4IntersectionSolid("Lens1", gubox, gcylinder1, xRot, zTrans1);
    G4SubtractionSolid *gLenst = new G4SubtractionSolid("temp", gubox, gcylinder1c, xRot, zTrans1);

    G4IntersectionSolid *gLens2 =
      new G4IntersectionSolid("Lens2", gLenst, gcylinder2, xRot, zTrans2);
    G4SubtractionSolid *gLens3 =
      new G4SubtractionSolid("Lens3", gLenst, gcylinder2c, xRot, zTrans2);

    lLens1 = new G4LogicalVolume(gLens1, BarMaterial, "lLens1", 0, 0, 0);
    lLens2 = new G4LogicalVolume(gLens2, Nlak33aMaterial, "lLens2", 0, 0, 0);
    lLens3 = new G4LogicalVolume(gLens3, BarMaterial, "lLens3", 0, 0, 0);
  }

  if (fLensId == 100) {
    fLens[2] = 200;
    G4Box *gLens3 = new G4Box("gLens1", fBar[0] / 2., 0.5 * fBoxWidth, 100);
    lLens3 = new G4LogicalVolume(gLens3, BarMaterial, "lLens3", 0, 0, 0);
  }

  if (fLensId != 0 && fLensId != 10) {
    double shifth = 0.5 * (dirclength + fLens[2]) + fBBWindow[2];
    if (fLensId == 100) {
      new G4PVPlacement(0, G4ThreeVector(rshift, 0, shifth), lLens3, "wLens3", lDirc, false, 0);
    } else if (fNBar == 1 && fLensId == 3) {
      for (int i = 0; i < 11; i++) {
        double shifty = i * fLens[1] - 0.5 * (fBoxWidth - fLens[1]);
        auto pos = G4ThreeVector(rshift, shifty, shifth);
        new G4PVPlacement(0, pos, lLens1, "wLens1", lDirc, false, i);
        new G4PVPlacement(0, pos, lLens2, "wLens2", lDirc, false, i);
        new G4PVPlacement(0, pos, lLens3, "wLens3", lDirc, false, i);
      }
    } else {
      double sh = 0;
      if (fEvType == 3) sh = evprismlengh + gluethickness;
      for (int i = 0; i < fNBar; i++) {
        double shifty = i * fLens[1] - 0.5 * (fBoxWidth - fLens[1]);
        if (fLensId != 6) {
          auto pos = G4ThreeVector(rshift, shifty, shifth - sh);
          new G4PVPlacement(0, pos, lLens1, "wLens1", lDirc, false, i);
          new G4PVPlacement(0, pos, lLens2, "wLens2", lDirc, false, i);
          if (fLensId == 3) new G4PVPlacement(0, pos, lLens3, "wLens3", lDirc, false, i);
        }
      }
      if (fLensId == 6) {
        auto pos = G4ThreeVector(rshift, 0, shifth - sh);
        new G4PVPlacement(0, pos, lLens1, "wLens1", lDirc, false, 0);
        new G4PVPlacement(0, pos, lLens2, "wLens2", lDirc, false, 0);
        new G4PVPlacement(0, pos, lLens3, "wLens3", lDirc, false, 0);
      }
    }
  }

  // The Prizm
  G4Trap *gPrizmT1 = new G4Trap("gPrizmT1", fPrizmT[0], fPrizmT[1], fPrizmT[2], fPrizmT[3]);
  lPrizmT1 = new G4LogicalVolume(gPrizmT1, BarMaterial, "lPrizmT1", 0, 0, 0);

  G4Trap *gPrizmT2 = new G4Trap("gPrizmT2", fPrizmT[0], fPrizmT[5], fPrizmT[2], 0.0001);
  lPrizmT2 = new G4LogicalVolume(gPrizmT2, BarMaterial, "lPrizmT2", 0, 0, 0);

  G4RotationMatrix *xRot = new G4RotationMatrix();
  xRot->rotateX(M_PI / 2. * rad);

  G4RotationMatrix *fdRot = new G4RotationMatrix();
  G4RotationMatrix *fdrot = new G4RotationMatrix();
  double evshiftz = 0.5 * dirclength + fPrizm[1] + fMcpActive[2] / 2. + fLens[2] + fBBWindow[2];
  double evshiftx = 0;

  if (fEvType == 1) { // focusing prism

    double fWindow[] = {124, 437, 14.0}; // 9.6
    double fWedge[] = {33.25, 91, 79.537, 27.0};
    double fSWedge[] = {fBoxWidth, 60, 123, 89};
    double fBlock[] = {fBoxWidth, 184.5, 495, 269};

    double currentz = 0.5 * dirclength;

    double fdH = fWedge[1] * tan(0.006); // 6 mrad tilt
    fWedge[2] = tan(30 / 180. * M_PI) * fWedge[1] + fWedge[3] -
                fdH; // update the side of the prizm assuming bottom tilt
    double theta = atan((fWedge[2] + 2 * fdH - fWedge[3]) / 2. / fWedge[1]);
    G4Trap *gWedge =
      new G4Trap("gWedge", fWedge[1] / 2., theta, 0., fWedge[0] / 2., fWedge[2] / 2.,
                 fWedge[2] / 2., 0., fWedge[0] / 2., fWedge[3] / 2., fWedge[3] / 2., 0.);
    lWedge = new G4LogicalVolume(gWedge, BarMaterial, "lWedge", 0, 0, 0);
    G4RotationMatrix *tRot = new G4RotationMatrix();
    tRot->rotateY(M_PI);

    // window
    G4Box *gWindow = new G4Box("gWindow", 0.5 * fWindow[0], 0.5 * fWindow[1], 0.5 * fWindow[2]);
    lWindow = new G4LogicalVolume(gWindow, BarMaterial, "lWindow", 0, 0, 0);

    double xstep = 0.5 * (fWedge[3] - fBar[0]);
    for (int i = 0; i < fNBar; i++) {
      double yshift = 0.5 * fBoxWidth - 0.5 * fBar[1] - (fBar[1] + fBarsGap) * i;
      fPrismShift = G4ThreeVector(rshift + (fWedge[2] + fWedge[3]) / 4. - 0.5 * fWedge[3] + xstep,
                                  yshift, currentz + 0.5 * fWedge[1]);
      new G4PVPlacement(tRot, fPrismShift, lWedge, "wWedge", lDirc, false, i);
    }
    currentz += fWedge[1];

    new G4PVPlacement(0, G4ThreeVector(rshift + 28 + xstep, 0, currentz + 0.5 * fWindow[2]),
                      lWindow, "wWindow", lDirc, false, 0);
    currentz += fWindow[2];

    // second wedge
    G4Trap *gSWedge = new G4Trap("gSWedge", fSWedge[0], fSWedge[1], fSWedge[2], fSWedge[3]);
    lSWedge = new G4LogicalVolume(gSWedge, BarMaterial, "lSWedge", 0, 0, 0);
    G4RotationMatrix *xwRot = new G4RotationMatrix();
    xwRot->rotateX(0.5 * M_PI);
    new G4PVPlacement(xwRot,
                      G4ThreeVector(rshift + (fSWedge[2] + fSWedge[3]) / 4. - 0.5 * fWedge[3] + xstep, 0,
                                    currentz + 0.5 * fSWedge[1]),
                      lSWedge, "wSWedge", lDirc, false, 0);
    currentz += fSWedge[1];

    // focusing block
    auto gBlockA = new G4Trap("gBlockA", fBlock[0], fBlock[1], fBlock[2], fBlock[3]);

    double cradius = 900; // PrtManager::Instance()->GetTest1();
    double cangle = 16;   // PrtManager::Instance()->GetTest2();

    double cpos = 240;
    double b = atan((cradius - sqrt(cradius * cradius - cpos * cpos)) / cpos) -
               cangle * TMath::DegToRad(); // atan(apos/cpos);
    double mxcor = cradius * (1 - cos(b));
    double mycor = cradius * sin(b);

    std::cout << b << " mxcor " << mxcor << " " << mycor << std::endl;

    G4Tubs *gCyl = new G4Tubs("gCyl", 0, cradius, fBoxWidth, 0 * deg, 360 * deg);
    G4RotationMatrix *cRot = new G4RotationMatrix();
    auto cTrans = G4ThreeVector(-mycor + cpos - ((fBlock[2] + fBlock[3]) / 4.),
                                mxcor - cradius + 0.5 * fBlock[1], 0);
    auto gBlock = new G4IntersectionSolid("gBlock", gBlockA, gCyl, cRot, cTrans);

    double mirrorthickness = 1;
    cTrans.setY(cTrans.getY() - mirrorthickness);

    auto gFmirror =
      new G4SubtractionSolid("gFmirror", gBlock, gCyl, new G4RotationMatrix(), cTrans);
    lFmirror = new G4LogicalVolume(gFmirror, MirrorMaterial, "lFmirror", 0, 0, 0);

    fPrismShift = G4ThreeVector(rshift + (fBlock[2] + fBlock[3]) / 4. - 0.5 * fWedge[3] + xstep, 0,
                                currentz + 0.5 * fBlock[1]);
    xRot->rotateX(M_PI);

    lBlock = new G4LogicalVolume(gBlock, BarMaterial, "lBlock", 0, 0, 0);
    new G4PVPlacement(xRot, fPrismShift, lBlock, "wBlock", lDirc, false, 0);

    fPrismShift.setZ(fPrismShift.getZ() + mirrorthickness);
    new G4PVPlacement(xRot, fPrismShift, lFmirror, "wFmirror", lDirc, false, 0);

    // PMT plane
    double pmtrot = atan((fBlock[3] - fBlock[2]) / fBlock[1]);
    fdrot->rotateY(0.5 * M_PI - pmtrot);
    new G4PVPlacement(
      fdrot,
      G4ThreeVector(rshift + 0.5 * (fBlock[2] + fBlock[3]) - 0.5 * fWedge[3] + xstep, 0,
                    currentz + 0.5 * fBlock[1]),
      lFd, "wFd", lDirc, false, 0);
  } else if (fEvType == 4) { // tilted prism

    fPrismShift = G4ThreeVector(rshift + (fPrizmT[2] + fPrizmT[3]) / 4. - fPrizmT[3] / 2., 0,
                                0.5 * dirclength + fPrizmT[1] / 2. + fLens[2]);
    new G4PVPlacement(xRot, fPrismShift, lPrizmT1, "wPrizm", lDirc, false, 0);
    fPrismShift = G4ThreeVector((fPrizmT[2] + 0.0001) / 4. - 0.0001 / 2. - 0.5 * fPrizm[3], 0,
                                fPrizmT[1] + 0.5 * dirclength + fPrizmT[5] / 2. + fLens[2]);
    G4RotationMatrix *yRot = new G4RotationMatrix();
    yRot->rotateX(-0.5 * M_PI * rad);
    fdRot->rotateY(fdTilt);
    fdrot->rotateY(fdTilt - 90 * deg);
    evshiftz =
      0.5 * dirclength + fLens[2] + fPrizmT[1] + 0.5 * fPrizmT[5] + 0.5 * fFd[2] * sin(fdTilt);
    evshiftx = 0.5 * fPrizmT[4] * (1 - sin(fdTilt)) - 0.5 * fFd[2] * cos(fdTilt);
    new G4PVPlacement(yRot, fPrismShift, lPrizmT2, "wPrizm", lDirc, false, 0);
    new G4PVPlacement(
      fdrot, G4ThreeVector(rshift + 0.5 * fFd[1] - 0.5 * fPrizm[3] - evshiftx, 0, evshiftz), lFd,
      "wFd", lDirc, false, 0);
  } else if (fEvType == 6){ // split prism
    
    G4Trap *gPrizm = new G4Trap("gPrizm", 0.5 * fPrizm[0] - 0.1, fPrizm[1], fPrizm[2], fPrizm[3]);
    lPrizm = new G4LogicalVolume(gPrizm, BarMaterial, "lPrizm", 0, 0, 0);
    fPrismShift = G4ThreeVector(rshift + (fPrizm[2] + fPrizm[3]) / 4. - 0.5 * fPrizm[3], -0.25 * fPrizm[0],
                                0.5 * (dirclength + fPrizm[1]) + fLens[2] + fBBWindow[2]);
    new G4PVPlacement(xRot, fPrismShift, lPrizm, "wPrizm", lDirc, false, 0);
    
    fPrismShift = G4ThreeVector(rshift + (fPrizm[2] + fPrizm[3]) / 4. - 0.5 * fPrizm[3], 0.25 * fPrizm[0],
                                0.5 * (dirclength + fPrizm[1]) + fLens[2] + fBBWindow[2]);
    new G4PVPlacement(xRot, fPrismShift, lPrizm, "wPrizm", lDirc, false, 1);

  } else { // normal prism

    G4Trap *gPrizm = new G4Trap("gPrizm", fPrizm[0], fPrizm[1], fPrizm[2], fPrizm[3]);
    lPrizm = new G4LogicalVolume(gPrizm, BarMaterial, "lPrizm", 0, 0, 0);
    
    fPrismShift = G4ThreeVector(rshift + (fPrizm[2] + fPrizm[3]) / 4. - 0.5 * fPrizm[3], 0,
                                0.5 * (dirclength + fPrizm[1]) + fLens[2] + fBBWindow[2]);
    new G4PVPlacement(xRot, fPrismShift, lPrizm, "wPrizm", lDirc, false, 0);
  }

  new G4PVPlacement(fdrot,
                    G4ThreeVector(rshift + 0.5 * fFd[1] - 0.5 * fPrizm[3] - evshiftx, 0, evshiftz),
                    lFd, "wFd", lDirc, false, 0);
    
  if (fMcpLayout == 1) {
    // standard mcp pmt layout
    // The MCP
    G4Box *gMcp = new G4Box("gMcp", fMcpTotal[0] / 2., fMcpTotal[1] / 2., fMcpTotal[2] / 2.);
    lMcp = new G4LogicalVolume(gMcp, BarMaterial, "lMcp", 0, 0, 0);

    // The MCP Pixel
    G4Box *gPixel = new G4Box("gPixel", fMcpActive[0] / 16., fMcpActive[1] / 16., 0.01);
    lPixel = new G4LogicalVolume(gPixel, BarMaterial, "lPixel", 0, 0, 0);

    for (int i = 0; i < 8; i++) {
      for (int j = 0; j < 8; j++) {
        double shiftx = i * (fMcpActive[0] / 8.) - fMcpActive[0] / 2. + fMcpActive[0] / 16.;
        double shifty = j * (fMcpActive[0] / 8.) - fMcpActive[0] / 2. + fMcpActive[0] / 16.;
        new G4PVPlacement(0, G4ThreeVector(shiftx, shifty, 0), lPixel, "wPixel", lMcp, false,
                          8 * i + j);
      }
    }

    int mcpId = 0;
    double gapx = (fPrizm[2] - 4 * fMcpTotal[0]) / 5.;
    double gapy = (fPrizm[0] - 6 * fMcpTotal[0]) / 7.;
    for (int i = 0; i < fNCol; i++) {
      for (int j = 0; j < fNRow; j++) {
        double shiftx = i * (fMcpTotal[0] + gapx) - fFd[1] / 2. + fMcpTotal[0] / 2. + gapx;
        double shifty = j * (fMcpTotal[0] + gapy) - fFd[0] / 2. + fMcpTotal[0] / 2. + gapy;
        new G4PVPlacement(0, G4ThreeVector(shiftx, shifty, 0), lMcp, "wMcp", lFd, false, mcpId++);
      }
    }
  }
  if (fMcpLayout == 0) {
    // The MCP
    G4Box *gMcp = new G4Box("gMcp", fPrizm[2] / 2., fPrizm[0] / 2., fMcpTotal[2] / 2.);
    lMcp = new G4LogicalVolume(gMcp, BarMaterial, "lMcp", 0, 0, 0);

    // The MCP Pixel
    fNpix1 = 4;
    fNpix2 = 6;
    G4Box *gPixel = new G4Box("gPixel", fMcpTotal[0] / 2., fMcpTotal[1] / 2., 0.01);
    lPixel = new G4LogicalVolume(gPixel, BarMaterial, "lPixel", 0, 0, 0);

    // double disx = (fPrizm[0]-fNpix2*fMcpTotal[0])/(double)fNpix2;
    double disy = (fPrizm[1] - fNpix1 * fMcpTotal[0]) / (double)(fNpix1 + 1);
    if (true) {
      for (int i = 0; i < fNpix1; i++) {
        for (int j = 0; j < fNpix2; j++) {
          double shiftx =
            i * (fMcpTotal[0] + disy / 2.) - fPrizm[2] / 2. + fMcpTotal[0] / 2. + disy / 2.;
          double shifty =
            j * (fMcpTotal[0] + disy / 2.) - fPrizm[0] / 2. + fMcpTotal[0] / 2. + disy / 2.;
          new G4PVPlacement(0, G4ThreeVector(shiftx, shifty, 0), lPixel, "wPixel", lMcp, false,
                            fNpix1 * i + j);
        }
      }
    } else {
      new G4PVPlacement(fdRot, G4ThreeVector(0, 0, 0), lPixel, "wPixel", lMcp, false, 1);
    }
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), lMcp, "wMcp", lFd, false, 1);
  }
  if (fMcpLayout == 3) {
    // one mcp layout
    G4Box *gMcp = new G4Box("gMcp", 0.5 * fPrizm[2], 0.5 * fPrizm[0], 0.5 * fMcpTotal[2]);
    lMcp = new G4LogicalVolume(gMcp, BarMaterial, "lMcp", 0, 0, 0);

    double pixSize = 3 * mm;

    fNpix1 = fPrizm[2] / pixSize - 1;
    fNpix2 = fPrizm[0] / pixSize - 1;

    G4Box *gPixel = new G4Box("gPixel", 0.5 * fPrizm[2] / fNpix1, 0.5 * fPrizm[0] / fNpix2, 0.01);
    lPixel = new G4LogicalVolume(gPixel, BarMaterial, "lPixel", 0, 0, 0);
    for (int i = 0; i < fNpix1; i++) {
      for (int j = 0; j < fNpix2; j++) {
        double shiftx = i * (fPrizm[2] / fNpix1) - 0.5 * fPrizm[2] + 0.5 * fPrizm[2] / fNpix1;
        double shifty = j * (fPrizm[0] / fNpix2) - 0.5 * fPrizm[0] + 0.5 * fPrizm[0] / fNpix2;
        new G4PVPlacement(fdRot, G4ThreeVector(shiftx, shifty, 0), lPixel, "wPixel", lMcp, false,
                          fNpix1 * i + j);
      }
    }
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), lMcp, "wMcp", lFd, false, 1);
  }
  if (fMcpLayout == 4) {
    // LAPPD/HRPPD pmt

    G4Box *gMcp = new G4Box("gMcp", fMcpTotal[0] / 2., fMcpTotal[1] / 2., fMcpTotal[2] / 2.);
    lMcp = new G4LogicalVolume(gMcp, BarMaterial, "lMcp", 0, 0, 0);

    fNpix1 = 40;
    fNpix2 = 40;
    fRun->setNpix(fNpix1 * fNpix1);

    std::cout << "fNpix1=" << fNpix1 << " fNpix2=" << fNpix2 << " size = " << fMcpActive[0] / fNpix1
              << std::endl;

    // The MCP Pixel
    G4Box *gPixel =
      new G4Box("gPixel", 0.5 * fMcpActive[0] / fNpix1, 0.5 * fMcpActive[1] / fNpix2, 0.01);
    lPixel = new G4LogicalVolume(gPixel, BarMaterial, "lPixel", 0, 0, 0);

    for (int i = 0; i < fNpix1; i++) {
      for (int j = 0; j < fNpix2; j++) {
        double shiftx =
          i * (fMcpActive[0] / fNpix1) - fMcpActive[0] / 2. + 0.5 * fMcpActive[0] / fNpix1;
        double shifty =
          j * (fMcpActive[1] / fNpix2) - fMcpActive[1] / 2. + 0.5 * fMcpActive[1] / fNpix2;
        new G4PVPlacement(0, G4ThreeVector(shiftx, shifty, 0), lPixel, "wPixel", lMcp, false,
                          fNpix1 * i + j);
      }
    }

    int mcpId = 0;
    double gapx = 0;
    double gapy = 0;

    for (int i = 0; i < fNCol; i++) {
      for (int j = 0; j < fNRow; j++) {
        double shiftx = i * (fMcpTotal[0] + gapx) - 0.5 * (fFd[1] - fMcpTotal[0]) + gapx +
                        0.5 * (fPrizm[2] - fNCol * fMcpTotal[1]);
        ;
        double shifty = j * (fMcpTotal[1] + gapy) - 0.5 * (fBoxWidth - fMcpTotal[1]) + gapy +
                        0.5 * (fBoxWidth - fNRow * fMcpTotal[1]);
        new G4PVPlacement(0, G4ThreeVector(shiftx, shifty, 0), lMcp, "wMcp", lFd, false, mcpId++);
      }
    }
  }
  if (fMcpLayout == 2031) {
    // alternative mcp pmt
    // The MCP

    // fMcpTotal[1]=fPrizm[0]/6.-1;
    // fMcpTotal[0]=fPrizm[2]/4.-5;
    // fMcpActive[1]= fMcpTotal[1];
    // fMcpActive[0]= fMcpTotal[0];

    G4Box *gMcp = new G4Box("gMcp", fMcpTotal[0] / 2., fMcpTotal[1] / 2., fMcpTotal[2] / 2.);
    lMcp = new G4LogicalVolume(gMcp, BarMaterial, "lMcp", 0, 0, 0);

    // double pixSize = 6*mm;
    // fNpix1 = 32; // fMcpActive[1]/pixSize-1;
    // fNpix2 = 32; // fMcpActive[1]/pixSize-1;

    fNpix1 = 16;
    fNpix2 = 16;
    fRun->setNpix(fNpix1 * fNpix1);
    // fRun->setTest2(fMcpActive[0] / fNpix1);
    // fRun->setNpix(fNpix1 * fNpix1);

    std::cout << "fNpix1=" << fNpix1 << " fNpix2=" << fNpix2 << " size = " << fMcpActive[0] / fNpix1
              << std::endl;

    // The MCP Pixel
    G4Box *gPixel =
      new G4Box("gPixel", 0.5 * fMcpActive[0] / fNpix1, 0.5 * fMcpActive[1] / fNpix2, 0.01);
    lPixel = new G4LogicalVolume(gPixel, BarMaterial, "lPixel", 0, 0, 0);

    for (int i = 0; i < fNpix1; i++) {
      for (int j = 0; j < fNpix2; j++) {
        double shiftx =
          i * (fMcpActive[0] / fNpix1) - fMcpActive[0] / 2. + 0.5 * fMcpActive[0] / fNpix1;
        double shifty =
          j * (fMcpActive[1] / fNpix2) - fMcpActive[1] / 2. + 0.5 * fMcpActive[1] / fNpix2;
        new G4PVPlacement(0, G4ThreeVector(shiftx, shifty, 0), lPixel, "wPixel", lMcp, false,
                          fNpix1 * i + j);
      }
    }

    int mcpId = 0;
    double gapx = (fPrizm[2] - fNCol * fMcpTotal[0]) / (double)(fNCol + 1);
    double gapy = (fBoxWidth - fNRow * fMcpTotal[1]) / (double)(fNRow + 1);


    fNRow = 6;    
    for (int i = 0; i < fNCol; i++) {
      for (int j = 0; j < fNRow; j++) {
        // double shiftx = i*(fMcpTotal[0]+gapx)-0.5*(fPrizm[3]-fMcpTotal[0])+gapx;
        // double shifty = j*(fMcpTotal[1]+gapy)-0.5*(fBoxWidth-fMcpTotal[1])+gapy;
        // new
        // G4PVPlacement(0,G4ThreeVector(shiftx,shifty,0.5*dirclength+fPrizm[1]+fMcpActive[2]/2.+fLens[2]),lMcp,"wMcp",
        // lDirc,false,mcpId++);

        double shiftx = i * (fMcpTotal[0] + gapx) - 0.5 * (fFd[1] - fMcpTotal[0]) + gapx;
        double shifty = j * (fMcpTotal[1] + gapy) - 0.5 * (fBoxWidth - fMcpTotal[1]) + gapy;
	// shifty += 0.5*fMcpTotal[0];
        new G4PVPlacement(0, G4ThreeVector(shiftx, shifty, 0), lMcp, "wMcp", lFd, false, mcpId++);
      }
    }
  }

  {
    const int num = 36;
    double WaveLength[num];
    double PhotonEnergy[num];
    double PMTReflectivity[num];
    double EfficiencyMirrors[num];
    const double LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;
    for (int i = 0; i < num; i++) {
      WaveLength[i] = (300 + i * 10) * nanometer;
      PhotonEnergy[num - (i + 1)] = LambdaE / WaveLength[i];
      PMTReflectivity[i] = 0.;
      EfficiencyMirrors[i] = 0;
    }

    /***************** QUANTUM EFFICIENCY OF BURLE AND HAMAMTSU PMT'S *****/

    // ideal pmt quantum efficiency
    double QuantumEfficiencyIdial[num] = {
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    // Burle PMT's
    double QuantumEfficiencyB[num] = {0.,   0.001, 0.002, 0.005, 0.01,  0.015, 0.02,  0.03, 0.04,
                                      0.05, 0.06,  0.07,  0.09,  0.1,   0.13,  0.15,  0.17, 0.2,
                                      0.24, 0.26,  0.28,  0.282, 0.284, 0.286, 0.288, 0.29, 0.28,
                                      0.26, 0.24,  0.22,  0.20,  0.18,  0.15,  0.13,  0.12, 0.10};

    // hamamatsu pmt quantum efficiency
    double QuantumEfficiencyPMT[num] = {
      0.001, 0.002, 0.004, 0.007, 0.011, 0.015, 0.020, 0.026, 0.033, 0.040, 0.045, 0.056,
      0.067, 0.085, 0.109, 0.129, 0.138, 0.147, 0.158, 0.170, 0.181, 0.188, 0.196, 0.203,
      0.206, 0.212, 0.218, 0.219, 0.225, 0.230, 0.228, 0.222, 0.217, 0.210, 0.199, 0.177};

    // these quantum efficiencies have to be multiplied by geometry
    //   efficiency of given PMT's
    //   for Hamamatsu by factor 0.7
    //   for Burle by factor 0.45
    for (int k = 0; k < 36; k++) {
      QuantumEfficiencyB[k] = QuantumEfficiencyB[k] * 0.45;
      QuantumEfficiencyPMT[k] = QuantumEfficiencyPMT[k] * .7;
    }

    // double QuantumEfficiency[num]=
    //    { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
    //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
    //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

    //  double QuantumEfficiencyPMT[num] =
    //    { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
    //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
    //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

    /* define quantum efficiency for burle PMT's - the same efficiency is
       assigned to pads and also to slots !!!! */

    // burle pmt - bigger slots => logicPad
    G4MaterialPropertiesTable *PhotocatodBurleMPT = new G4MaterialPropertiesTable();
    PhotocatodBurleMPT->AddProperty("EFFICIENCY", PhotonEnergy, QuantumEfficiencyB, num);
    PhotocatodBurleMPT->AddProperty("REFLECTIVITY", PhotonEnergy, PMTReflectivity, num);

    G4OpticalSurface *BurlePMTOpSurface =
      new G4OpticalSurface("BurlePMTOpSurface", glisur, polished, dielectric_metal);
    BurlePMTOpSurface->SetMaterialPropertiesTable(PhotocatodBurleMPT);

    // // assignment for pad
    // if(burle)
    //   new G4LogicalSkinSurface("BurlePMTSurface",logicBurPad,BurlePMTOpSurface);

    // if(burle1)
    //   new G4LogicalSkinSurface("Burle1PMTSurface",logicBur1Pad,BurlePMTOpSurface);

    /* hamamatsu pmt's - smaller slots => quantum efficiency again
       assign to slot and pad */

    fQuantumEfficiency = QuantumEfficiencyIdial;
    G4MaterialPropertiesTable *PhotocatodHamamatsuMPT = new G4MaterialPropertiesTable();
    PhotocatodHamamatsuMPT->AddProperty("EFFICIENCY", PhotonEnergy, fQuantumEfficiency, num);
    PhotocatodHamamatsuMPT->AddProperty("REFLECTIVITY", PhotonEnergy, PMTReflectivity, num);

    G4OpticalSurface *HamamatsuPMTOpSurface =
      new G4OpticalSurface("HamamatsuPMTOpSurface", glisur, polished, dielectric_dielectric);
    HamamatsuPMTOpSurface->SetMaterialPropertiesTable(PhotocatodHamamatsuMPT);

    // // assignment to pad
    // if(hamamatsu8500)
    new G4LogicalSkinSurface("HamamatsuPMTSurface", lPixel, HamamatsuPMTOpSurface);

    // Mirror
    G4OpticalSurface *MirrorOpSurface =
      new G4OpticalSurface("MirrorOpSurface", glisur, polished, dielectric_metal);

    double ReflectivityMirrorBar[num] = {
      0.8755, 0.882,  0.889, 0.895, 0.9,   0.9025, 0.91,   0.913, 0.9165, 0.92,   0.923,  0.9245,
      0.9285, 0.932,  0.934, 0.935, 0.936, 0.9385, 0.9395, 0.94,  0.9405, 0.9405, 0.9405, 0.9405,
      0.94,   0.9385, 0.936, 0.934, 0.931, 0.9295, 0.928,  0.928, 0.921,  0.92,   0.927,  0.9215};

    G4MaterialPropertiesTable *MirrorMPT = new G4MaterialPropertiesTable();
    MirrorMPT->AddProperty("REFLECTIVITY", PhotonEnergy, ReflectivityMirrorBar, num);
    MirrorMPT->AddProperty("EFFICIENCY", PhotonEnergy, EfficiencyMirrors, num);

    MirrorOpSurface->SetMaterialPropertiesTable(MirrorMPT);
    new G4LogicalSkinSurface("MirrorSurface", lMirror, MirrorOpSurface);
    new G4LogicalSkinSurface("MirrorSurfaceF", lFmirror, MirrorOpSurface);
  }

  SetVisualization();

  return wExpHall;
}

void PrtDetectorConstruction::DefineMaterials() {
  G4String symbol;      // a=mass of a mole;
  double a, z, density; // z=mean number of protons;

  int ncomponents, natoms;
  double fractionmass;

  // define Elements
  G4Element *H = new G4Element("Hydrogen", symbol = "H", z = 1., a = 1.01 * g / mole);
  G4Element *C = new G4Element("Carbon", symbol = "C", z = 6., a = 12.01 * g / mole);
  G4Element *N = new G4Element("Nitrogen", symbol = "N", z = 7., a = 14.01 * g / mole);
  G4Element *O = new G4Element("Oxygen", symbol = "O", z = 8., a = 16.00 * g / mole);
  G4Element *Si = new G4Element("Silicon", symbol = "Si", z = 14., a = 28.09 * g / mole);
  G4Element *Al = new G4Element("Aluminum", symbol = "Al", z = 13., a = 26.98 * g / mole);

  // quartz material = SiO2
  G4Material *SiO2 = new G4Material("quartz", density = 2.200 * g / cm3, ncomponents = 2);
  SiO2->AddElement(Si, natoms = 1);
  SiO2->AddElement(O, natoms = 2);

  Nlak33aMaterial = new G4Material("Nlak33a", density = 4.220 * g / cm3, ncomponents = 2);
  Nlak33aMaterial->AddElement(Si, natoms = 1);
  Nlak33aMaterial->AddElement(O, natoms = 2);

  PbF2Material = new G4Material("PbF2", density = 3.97 * g / cm3, ncomponents = 2);
  PbF2Material->AddElement(Si, natoms = 1);
  PbF2Material->AddElement(O, natoms = 2);

  SapphireMaterial = new G4Material("Sapphire", density = 3.98 * g / cm3, ncomponents = 2);
  SapphireMaterial->AddElement(Al, natoms = 2);
  SapphireMaterial->AddElement(O, natoms = 3);

  G4Material *Vacuum = new G4Material("interGalactic", 1., 1.008 * g / mole, 1.e-25 * g / cm3,
                                      kStateGas, 2.73 * kelvin, 3.e-18 * pascal);
  G4Material *Air = new G4Material("Air", density = 1.290 * mg / cm3, ncomponents = 2);
  Air->AddElement(N, fractionmass = 0.7);
  Air->AddElement(O, fractionmass = 0.3);

  G4Material *Aluminum = new G4Material("Aluminum", density = 2.7 * g / cm3, ncomponents = 1);
  Aluminum->AddElement(Al, fractionmass = 1.0);

  G4Material *KamLandOil = new G4Material("KamLandOil", density = 0.914 * g / cm3, ncomponents = 2);
  KamLandOil->AddElement(C, natoms = 12);
  KamLandOil->AddElement(H, natoms = 26);

  G4Material *CarbonFiber =
    new G4Material("CarbonFiber", density = 0.145 * g / cm3, ncomponents = 1);
  CarbonFiber->AddElement(C, natoms = 1);

  /* as I don't know the exact material composition,
     I will use Epoxyd material composition and add
     the optical property of Epotek to this material */

  G4Material *Epotek = new G4Material("Epotek", density = 1.2 * g / cm3, ncomponents = 3);

  Epotek->AddElement(C, natoms = 3);
  Epotek->AddElement(H, natoms = 5);
  Epotek->AddElement(O, natoms = 2);

  // assign main materials
  if (fGeomType < 20) defaultMaterial = Vacuum;
  else defaultMaterial = Air; // Vacuum // material of world

  frontMaterial = CarbonFiber;
  BarMaterial = SiO2;        // material of all Bars, Quartz and Window
  OilMaterial = KamLandOil;  // material of volume 1,2,3,4
  MirrorMaterial = Aluminum; // mirror material
  epotekMaterial = Epotek;   // Epotek material - glue between bars

  // ------------ Generate & Add Material Properties Table ------------

  static const double LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;
  const int num = 36;
  double WaveLength[num];
  double AirAbsorption[num];      // absorption value for air
  double AirRefractiveIndex[num]; // air refractive index
  double PhotonEnergy[num];       // energy of photons which correspond to the given

  // refractive or absoprtion values
  double PhotonEnergyNlak33a[76] = {
    1 * eV,       1.2511 * eV,  1.26386 * eV, 1.27687 * eV, 1.29016 * eV, 1.30372 * eV,
    1.31758 * eV, 1.33173 * eV, 1.34619 * eV, 1.36097 * eV, 1.37607 * eV, 1.39152 * eV,
    1.40731 * eV, 1.42347 * eV, 1.44 * eV,    1.45692 * eV, 1.47425 * eV, 1.49199 * eV,
    1.51016 * eV, 1.52878 * eV, 1.54787 * eV, 1.56744 * eV, 1.58751 * eV, 1.6081 * eV,
    1.62923 * eV, 1.65092 * eV, 1.6732 * eV,  1.69609 * eV, 1.71961 * eV, 1.7438 * eV,
    1.76868 * eV, 1.79427 * eV, 1.82062 * eV, 1.84775 * eV, 1.87571 * eV, 1.90452 * eV,
    1.93423 * eV, 1.96488 * eV, 1.99652 * eV, 2.0292 * eV,  2.06296 * eV, 2.09787 * eV,
    2.13398 * eV, 2.17135 * eV, 2.21006 * eV, 2.25017 * eV, 2.29176 * eV, 2.33492 * eV,
    2.37973 * eV, 2.42631 * eV, 2.47473 * eV, 2.52514 * eV, 2.57763 * eV, 2.63236 * eV,
    2.68946 * eV, 2.7491 * eV,  2.81143 * eV, 2.87666 * eV, 2.94499 * eV, 3.01665 * eV,
    3.09187 * eV, 3.17095 * eV, 3.25418 * eV, 3.34189 * eV, 3.43446 * eV, 3.53231 * eV,
    3.6359 * eV,  3.74575 * eV, 3.86244 * eV, 3.98663 * eV, 4.11908 * eV, 4.26062 * eV,
    4.41225 * eV, 4.57506 * eV, 4.75035 * eV, 4.93961 * eV};

  int n_PbF2 = 56;
  double en_PbF2[] = {
    1.55 * eV,  1.569 * eV, 1.59 * eV,  1.61 * eV,  1.631 * eV, 1.653 * eV, 1.675 * eV, 1.698 * eV,
    1.722 * eV, 1.746 * eV, 1.771 * eV, 1.797 * eV, 1.823 * eV, 1.851 * eV, 1.879 * eV, 1.907 * eV,
    1.937 * eV, 1.968 * eV, 2 * eV,     2.033 * eV, 2.066 * eV, 2.101 * eV, 2.138 * eV, 2.175 * eV,
    2.214 * eV, 2.254 * eV, 2.296 * eV, 2.339 * eV, 2.384 * eV, 2.431 * eV, 2.48 * eV,  2.53 * eV,
    2.583 * eV, 2.638 * eV, 2.695 * eV, 2.755 * eV, 2.818 * eV, 2.883 * eV, 2.952 * eV, 3.024 * eV,
    3.1 * eV,   3.179 * eV, 3.263 * eV, 3.351 * eV, 3.444 * eV, 3.542 * eV, 3.647 * eV, 3.757 * eV,
    3.875 * eV, 3.999 * eV, 4.133 * eV, 4.275 * eV, 4.428 * eV, 4.592 * eV, 4.769 * eV, 4.959 * eV};
  double ab_PbF2[] = {407,   403.3, 379.1, 406.3, 409.7, 408.9, 406.7, 404.7, 391.7, 397.7,
                      409.6, 403.7, 403.8, 409.7, 404.9, 404.2, 407.1, 411.1, 403.1, 406.1,
                      415.4, 399.1, 405.8, 408.2, 385.7, 405.6, 405.2, 401.6, 402.6, 407.1,
                      417.7, 401.1, 389.9, 411.9, 400.9, 398.3, 402.1, 408.7, 384.8, 415.8,
                      413.1, 385.7, 353.7, 319.1, 293.6, 261.9, 233.6, 204.4, 178.3, 147.6,
                      118.2, 78.7,  51.6,  41.5,  24.3,  8.8};
  for (int i = 0; i < n_PbF2; i++) {
    ab_PbF2[i] *= cm; // cm to mm
  }
  double ref_PbF2[] = {1.749, 1.749, 1.75,  1.75,  1.751, 1.752, 1.752, 1.753, 1.754, 1.754,
                       1.755, 1.756, 1.757, 1.757, 1.758, 1.759, 1.76,  1.761, 1.762, 1.764,
                       1.765, 1.766, 1.768, 1.769, 1.771, 1.772, 1.774, 1.776, 1.778, 1.78,
                       1.782, 1.785, 1.787, 1.79,  1.793, 1.796, 1.8,   1.804, 1.808, 1.813,
                       1.818, 1.824, 1.83,  1.837, 1.845, 1.854, 1.865, 1.877, 1.892, 1.91,
                       1.937, 1.991, 1.38,  1.915, 1.971, 2.019};

  const int n_Sapphire = 57;
  double en_Sapphire[] = {
    1.02212 * eV, 1.05518 * eV, 1.09045 * eV, 1.12610 * eV, 1.16307 * eV, 1.20023 * eV,
    1.23984 * eV, 1.28043 * eV, 1.32221 * eV, 1.36561 * eV, 1.41019 * eV, 1.45641 * eV,
    1.50393 * eV, 1.55310 * eV, 1.60393 * eV, 1.65643 * eV, 1.71059 * eV, 1.76665 * eV,
    1.82437 * eV, 1.88397 * eV, 1.94576 * eV, 2.00946 * eV, 2.07504 * eV, 2.14283 * eV,
    2.21321 * eV, 2.28542 * eV, 2.36025 * eV, 2.43727 * eV, 2.51693 * eV, 2.59924 * eV,
    2.68480 * eV, 2.77245 * eV, 2.86337 * eV, 2.95693 * eV, 3.05379 * eV, 3.15320 * eV,
    3.25674 * eV, 3.36273 * eV, 3.47294 * eV, 3.58646 * eV, 3.70433 * eV, 3.82549 * eV,
    3.94979 * eV, 4.07976 * eV, 4.21285 * eV, 4.35032 * eV, 4.49380 * eV, 4.64012 * eV,
    4.79258 * eV, 4.94946 * eV, 5.11064 * eV, 5.27816 * eV, 5.44985 * eV, 5.62797 * eV,
    5.81266 * eV, 6.00407 * eV, 6.19920 * eV};
  double ref_Sapphire[] = {
    1.75188, 1.75253, 1.75319, 1.75382, 1.75444, 1.75505, 1.75567, 1.75629, 1.75691, 1.75754,
    1.75818, 1.75883, 1.75949, 1.76017, 1.76088, 1.76160, 1.76235, 1.76314, 1.76395, 1.76480,
    1.76570, 1.76664, 1.76763, 1.76867, 1.76978, 1.77095, 1.77219, 1.77350, 1.77490, 1.77639,
    1.77799, 1.77968, 1.78150, 1.78343, 1.78551, 1.78772, 1.79011, 1.79265, 1.79540, 1.79835,
    1.80154, 1.80497, 1.80864, 1.81266, 1.81696, 1.82163, 1.82674, 1.83223, 1.83825, 1.84480,
    1.85191, 1.85975, 1.86829, 1.87774, 1.88822, 1.89988, 1.91270};

  double ab_Sapphire[n_Sapphire];
  // transmittance for 2 mm thickness (crystan standard UV grade)
  // https://www.crystran.co.uk/optical-materials/sapphire-al2o3
  double transmittance_Sapphire[] = {
    0.845, 0.844, 0.843, 0.843, 0.843, 0.843, 0.843, 0.843, 0.843, 0.843, 0.843, 0.843,
    0.843, 0.843, 0.843, 0.843, 0.843, 0.843, 0.843, 0.843, 0.842, 0.842, 0.842, 0.842,
    0.838, 0.838, 0.838, 0.838, 0.838, 0.838, 0.836, 0.836, 0.836, 0.836, 0.834, 0.834,
    0.833, 0.832, 0.832, 0.831, 0.831, 0.83,  0.828, 0.828, 0.828, 0.828, 0.828, 0.828,
    0.828, 0.828, 0.828, 0.825, 0.80,  0.67,  0.47,  0.23,  0.22};

  // conversion to mean free path
  for (int i = 0; i < n_Sapphire; i++) {
    double fresnel = (ref_Sapphire[i] - 1) * (ref_Sapphire[i] - 1) /
                     ((ref_Sapphire[i] + 1) * (ref_Sapphire[i] + 1));
    ab_Sapphire[i] = (-1) / log(transmittance_Sapphire[i] + 2 * fresnel) * 2 * mm;
  }
  /*************************** ABSORPTION COEFFICIENTS *****************************/

  // absorption of KamLandOil per 50 cm - from jjv
  double KamLandOilAbsorption[num] = {
    0.97469022,  0.976603956, 0.978511548, 0.980400538, 0.982258449, 0.984072792,
    0.985831062, 0.987520743, 0.989129303, 0.990644203, 0.992052894, 0.993342822,
    0.994501428, 0.995516151, 0.996374433, 0.997063719, 0.997571464, 0.997885132,
    0.997992205, 0.997880183, 0.997536591, 0.99,        0.98,        0.97,
    0.96,        0.94,        0.93,        0.924507,    0.89982,     0.883299,
    0.85657,     0.842637,    0.77020213,  0.65727,     0.324022,    0.019192};

  // absorption of quartz per 1 m - from jjv
  double QuartzAbsorption[num] = {
    0.999572036, 0.999544661, 0.999515062, 0.999483019, 0.999448285, 0.999410586,
    0.999369611, 0.999325013, 0.999276402, 0.999223336, 0.999165317, 0.999101778,
    0.999032079, 0.998955488, 0.998871172, 0.998778177, 0.99867541,  0.998561611,
    0.998435332, 0.998294892, 0.998138345, 0.997963425, 0.997767484, 0.997547418,
    0.99729958,  0.99701966,  0.99670255,  0.996342167, 0.995931242, 0.995461041,
    0.994921022, 0.994298396, 0.993577567, 0.992739402, 0.991760297, 0.990610945};

  // absorption of epotek per one layer - thicknes 0.001'' - from jjv
  double EpotekAbsorption[num] = {
    0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999,
    0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999,
    0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999,
    0.99999999, 0.9999,     0.9998,     0.9995,     0.999,      0.998,      0.997,      0.996,
    0.9955,     0.993,      0.9871,     0.9745};

  // N-Lak 33a
  double Nlak33aAbsorption[76] = {
    371813,    352095,    331021,    310814,    291458,  272937,  255238,  238342,   222234,
    206897,    192313,    178463,    165331,    152896,  141140,  130043,  119585,   109747,
    100507,    91846.3,   83743.1,   76176.7,   69126.1, 62570.2, 56488,   50858.3,  45660.1,
    40872.4,   36474.6,   32445.8,   28765.9,   25414.6, 22372.2, 19619.3, 17136.9,  14906.5,
    12910.2,   11130.3,   9550.13,   8153.3,    6924.25, 5848.04, 4910.46, 4098.04,  3398.06,
    2798.54,   2288.32,   1856.99,   1494.92,   1193.28, 943.973, 739.657, 573.715,  440.228,
    333.94,    250.229,   185.064,   134.967,   96.9664, 68.5529, 47.6343, 32.4882,  21.7174,
    14.2056,   9.07612,   5.65267,   3.4241,    2.01226, 1.14403, 0.62722, 0.330414, 0.166558,
    0.0799649, 0.0363677, 0.0155708, 0.00623089};
  
  for (int i = 0; i < 76; i++) {
    Nlak33aAbsorption[i] *= cm; // cm to mm
  }

  double EpotekThickness = 0.001 * 2.54 * cm;
  for (int i = 0; i < num; i++) {
    WaveLength[i] = (300 + i * 10) * nanometer;
    AirAbsorption[i] = 4. * cm; // if photon in the air -> kill it immediately
    AirRefractiveIndex[i] = 1.;
    PhotonEnergy[num - (i + 1)] = LambdaE / WaveLength[i];

    /* as the absorption is given per length and G4 needs
       mean free path length, calculate it here
       mean free path length - taken as probability equal 1/e
       that the photon will be absorbed */

    EpotekAbsorption[i] = (-1) / log(EpotekAbsorption[i]) * EpotekThickness;
    QuartzAbsorption[i] = (-1) / log(QuartzAbsorption[i]) * 100 * cm;
    KamLandOilAbsorption[i] = (-1) / log(KamLandOilAbsorption[i]) * 50 * cm;
  }

  /**************************** REFRACTIVE INDEXES ****************************/

  // only phase refractive indexes are necessary -> g4 calculates group itself !!

  double QuartzRefractiveIndex[num] = {
    1.456535, 1.456812, 1.4571,   1.457399, 1.457712, 1.458038, 1.458378, 1.458735, 1.459108,
    1.4595,   1.459911, 1.460344, 1.460799, 1.46128,  1.461789, 1.462326, 1.462897, 1.463502,
    1.464146, 1.464833, 1.465566, 1.46635,  1.46719,  1.468094, 1.469066, 1.470116, 1.471252,
    1.472485, 1.473826, 1.475289, 1.476891, 1.478651, 1.480592, 1.482739, 1.485127, 1.487793};

  double EpotekRefractiveIndex[num] = {
    1.554034, 1.555575, 1.55698,  1.558266, 1.559454, 1.56056,  1.561604, 1.562604, 1.563579,
    1.564547, 1.565526, 1.566536, 1.567595, 1.568721, 1.569933, 1.57125,  1.57269,  1.574271,
    1.576012, 1.577932, 1.580049, 1.582381, 1.584948, 1.587768, 1.590859, 1.59424,  1.597929,
    1.601946, 1.606307, 1.611033, 1.616141, 1.621651, 1.62758,  1.633947, 1.640771, 1.64807};

  double KamLandOilRefractiveIndex[num] = {
    1.433055, 1.433369, 1.433698, 1.434045, 1.434409, 1.434793, 1.435198, 1.435626, 1.436077,
    1.436555, 1.4371,   1.4376,   1.4382,   1.4388,   1.4395,   1.4402,   1.4409,   1.4415,
    1.4425,   1.4434,   1.4444,   1.4455,   1.4464,   1.4479,   1.4501,   1.450428, 1.451976,
    1.453666, 1.455513, 1.45754,  1.45977,  1.462231, 1.464958, 1.467991, 1.471377, 1.475174};

  double Nlak33aRefractiveIndex[76] = {
    1.73816, 1.73836, 1.73858, 1.73881, 1.73904, 1.73928, 1.73952, 1.73976, 1.74001, 1.74026,
    1.74052, 1.74078, 1.74105, 1.74132, 1.7416,  1.74189, 1.74218, 1.74249, 1.74279, 1.74311,
    1.74344, 1.74378, 1.74412, 1.74448, 1.74485, 1.74522, 1.74562, 1.74602, 1.74644, 1.74687,
    1.74732, 1.74779, 1.74827, 1.74878, 1.7493,  1.74985, 1.75042, 1.75101, 1.75163, 1.75228,
    1.75296, 1.75368, 1.75443, 1.75521, 1.75604, 1.75692, 1.75784, 1.75882, 1.75985, 1.76095,
    1.76211, 1.76335, 1.76467, 1.76608, 1.76758, 1.7692,  1.77093, 1.77279, 1.7748,  1.77698,
    1.77934, 1.7819,  1.7847,  1.78775, 1.79111, 1.79481, 1.79889, 1.80343, 1.8085,  1.81419,
    1.82061, 1.8279,  1.83625, 1.84589, 1.85713, 1.87039};

  /* ASSIGNING REFRACTIVE AND ABSORPTION PROPERTIES TO THE GIVEN MATERIALS */

  // Quartz material => Si02
  G4MaterialPropertiesTable *QuartzMPT = new G4MaterialPropertiesTable();
  QuartzMPT->AddProperty("RINDEX", PhotonEnergy, QuartzRefractiveIndex, num);
  QuartzMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);

  // assign this parameter table to BAR material
  BarMaterial->SetMaterialPropertiesTable(QuartzMPT);

  // Air
  G4MaterialPropertiesTable *AirMPT = new G4MaterialPropertiesTable();
  AirMPT->AddProperty("RINDEX", PhotonEnergy, AirRefractiveIndex, num);
  AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAbsorption, num);
  //  assign this parameter table to the air
  defaultMaterial->SetMaterialPropertiesTable(AirMPT);

  // KamLandOil
  G4MaterialPropertiesTable *KamLandOilMPT = new G4MaterialPropertiesTable();
  KamLandOilMPT->AddProperty("RINDEX", PhotonEnergy, KamLandOilRefractiveIndex, num);
  KamLandOilMPT->AddProperty("ABSLENGTH", PhotonEnergy, KamLandOilAbsorption, num);
  // assing this parameter table  to the KamLandOil
  OilMaterial->SetMaterialPropertiesTable(KamLandOilMPT);

  // N-Lak 33a
  G4MaterialPropertiesTable *Nlak33aMPT = new G4MaterialPropertiesTable();
  Nlak33aMPT->AddProperty("RINDEX", PhotonEnergyNlak33a, Nlak33aRefractiveIndex, 76);
  Nlak33aMPT->AddProperty("ABSLENGTH", PhotonEnergyNlak33a, Nlak33aAbsorption, 76);
  Nlak33aMaterial->SetMaterialPropertiesTable(Nlak33aMPT);

  // PbF2
  G4MaterialPropertiesTable *PbF2MPT = new G4MaterialPropertiesTable();
  PbF2MPT->AddProperty("RINDEX", en_PbF2, ref_PbF2, n_PbF2);
  PbF2MPT->AddProperty("ABSLENGTH", en_PbF2, ab_PbF2, n_PbF2);
  PbF2Material->SetMaterialPropertiesTable(PbF2MPT);

  // Sapphire
  G4MaterialPropertiesTable *SapphireMPT = new G4MaterialPropertiesTable();
  SapphireMPT->AddProperty("RINDEX", en_Sapphire, ref_Sapphire, n_Sapphire);
  SapphireMPT->AddProperty("ABSLENGTH", en_Sapphire, ab_Sapphire, n_Sapphire);
  SapphireMaterial->SetMaterialPropertiesTable(SapphireMPT);

  // Epotek Glue
  G4MaterialPropertiesTable *EpotekMPT = new G4MaterialPropertiesTable();
  EpotekMPT->AddProperty("RINDEX", PhotonEnergy, EpotekRefractiveIndex, num);
  EpotekMPT->AddProperty("ABSLENGTH", PhotonEnergy, EpotekAbsorption, num);
  // assign this parameter table to the epotek
  epotekMaterial->SetMaterialPropertiesTable(EpotekMPT);
}

void PrtDetectorConstruction::SetVisualization() {

  G4Colour DircColour = G4Colour(1., 1.0, 0.);

  G4VisAttributes *waExpHall = new G4VisAttributes(DircColour);
  waExpHall->SetVisibility(false);
  waExpHall->SetForceWireframe(true);
  lExpHall->SetVisAttributes(waExpHall);

  G4VisAttributes *waDirc = new G4VisAttributes(DircColour);
  waDirc->SetForceWireframe(true);
  waDirc->SetVisibility(false);
  lDirc->SetVisAttributes(waDirc);

  G4VisAttributes *waFd = new G4VisAttributes(DircColour);
  waFd->SetForceWireframe(true);
  lFd->SetVisAttributes(waFd);

  G4VisAttributes *waBar = new G4VisAttributes(G4Colour(0., 1., 0.9, 0.05)); // 0.05
  waBar->SetVisibility(true);
  if(fBBWindow[2] > 0.1) lBBWindow->SetVisAttributes(waBar);
  lBar->SetVisAttributes(waBar);
  lExpVol->SetVisAttributes(waBar);
  
  G4VisAttributes *waGlue = new G4VisAttributes(G4Colour(0., 0.4, 0.9, 0.1));
  waGlue->SetVisibility(true);
  lGlue->SetVisAttributes(waGlue);
  lGlueE->SetVisAttributes(waGlue);

  G4VisAttributes *waMirror = new G4VisAttributes(G4Colour(1., 1., 0.9, 0.2));
  waMirror->SetVisibility(true);
  lMirror->SetVisAttributes(waMirror);

  G4VisAttributes *waTracker = new G4VisAttributes(G4Colour(1., 1., 0.7, 0.2));
  waTracker->SetVisibility(false);
  lTracker->SetVisAttributes(waTracker);
  
  double transp = 0.4;
  G4VisAttributes *vaLens = new G4VisAttributes(G4Colour(0., 1., 1., transp));
  vaLens->SetForceWireframe(true);
  // vaLens->SetForceAuxEdgeVisible(true);
  // vaLens->SetForceLineSegmentsPerCircle(10);
  // vaLens->SetLineWidth(4);

  if (fLensId == 100) lLens3->SetVisAttributes(vaLens);

  if (fLensId == 2 || fLensId == 3 || fLensId == 6) {
    lLens1->SetVisAttributes(vaLens);
    G4VisAttributes *vaLens2 = new G4VisAttributes(G4Colour(0., 1., 1., transp));
    vaLens2->SetColour(G4Colour(0., 0.5, 1., transp));
    vaLens2->SetForceWireframe(true);
    lLens2->SetVisAttributes(vaLens2);
    if (fLensId == 3 || fLensId == 6) lLens3->SetVisAttributes(vaLens);
  }

  G4VisAttributes *waPrizm = new G4VisAttributes(G4Colour(0., 0.9, 0.9, 0.4)); // 0.4
  // waPrizm->SetForceAuxEdgeVisible(true);
  // waPrizm->SetForceSolid(true);
  lPrizm->SetVisAttributes(waPrizm);
  lPrizmT1->SetVisAttributes(waPrizm);
  lPrizmT2->SetVisAttributes(waPrizm);

  if (fEvType == 1) {
    lBlock->SetVisAttributes(waPrizm);
    lWindow->SetVisAttributes(waPrizm);
    lWedge->SetVisAttributes(waBar);
    lSWedge->SetVisAttributes(waPrizm);
    lFmirror->SetVisAttributes(waMirror);
  }

  G4VisAttributes *waMcp = new G4VisAttributes(G4Colour(0.3,0.0,0.3,0.1));
  //G4VisAttributes *waMcp = new G4VisAttributes(G4Colour(1.0, 0., 0.1, 0.4));
  // waMcp->SetForceWireframe(true);
  waMcp->SetForceSolid(true);
  lMcp->SetVisAttributes(waMcp);

  // G4VisAttributes *waPixel = new G4VisAttributes(G4Colour(0.7, 0.0, 0.1, 0.5));
 G4VisAttributes *waPixel = new G4VisAttributes(G4Colour(0.1, 0.4, 0.6, 0.3));
  waPixel->SetForceWireframe(true);
  if (fMcpLayout == 3) waPixel->SetVisibility(false);
  else waPixel->SetVisibility(true);
  lPixel->SetVisAttributes(waPixel);
}

void PrtDetectorConstruction::ConstructSDandField() {

  // smearing at the tracker
  G4Region *r = new G4Region("BarRegion");
  r->AddRootLogicalVolume(lTracker);
  auto fastSimModelTracker = new PrtFastSimModelTracker("tracker", r);

  // Register the fast simulation model for deleting
  G4AutoDelete::Register(fastSimModelTracker);

  // Sensitive detectors
  PrtPixelSD *pixelSD = new PrtPixelSD("PixelSD", "PixelHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(pixelSD);
  SetSensitiveDetector("lPixel", pixelSD);
  PrtPrizmSD *prizmSD = new PrtPrizmSD("PrizmSD", "PrizmHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(prizmSD);

  if (fEvType == 0) {
    SetSensitiveDetector("lPrizm", prizmSD);
    if (fLensId == 2 || fLensId == 3 || fLensId == 6) {
      SetSensitiveDetector("lLens1", prizmSD);
      SetSensitiveDetector("lLens2", prizmSD);
      SetSensitiveDetector("lLens3", prizmSD);
    }
  } else if (fEvType == 1) {
    SetSensitiveDetector("lWedge", prizmSD);
    SetSensitiveDetector("lSWedge", prizmSD);
    SetSensitiveDetector("lBlock", prizmSD);
  } else if (fEvType == 4) {
    SetSensitiveDetector("lPrizmT1", prizmSD);
  }

  PrtBarSD *barSD = new PrtBarSD("BarSD", "BarHitsCollection", 0);
  G4SDManager::GetSDMpointer()->AddNewDetector(barSD);
  SetSensitiveDetector("lBar", barSD);

  // Magnetic Fields
  int fid = fRun->getField();
  if (fid == 1 || fid == 2 || fid == 3 || fid == 4) {

    double ZOffset = 0.;
    G4MagneticField *mField;
    if (fid == 1) {
      // CORE field
      mField = new PrtField("../data/field.tab", ZOffset);
    } else if (fid == 2) {
      // MARCO 1.7T
      mField = new PrtField(
        "../data/MARCO_v.6.4.1.1.3_1.7T_Magnetic_Field_Map_2022_11_14_rad_coords_cm_T.CART.txt",
        ZOffset);
    } else if (fid == 3) {
      // MARCO 2.0T
      mField = new PrtField(
        "../data/MARCO_v.6.4.1.1.3_2T_Magnetic_Field_Map_2022_11_14_rad_coords_cm_T.CART.txt",
        ZOffset);
    } else if (fid == 4) {
      // simple solenoidal field
      mField = new G4UniformMagField(G4ThreeVector(3. * tesla, 0., 0.));
    }

    G4FieldManager *globalFieldMgr =
      G4TransportationManager::GetTransportationManager()->GetFieldManager();
    globalFieldMgr->SetDetectorField(mField);
    globalFieldMgr->CreateChordFinder(mField);
  }
}

void PrtDetectorConstruction::SetRotation(double angle) {
  fPrtRot->rotateY(-fRotAngle);
  fPrtRot->rotateY(angle);
  fRotAngle = angle;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void PrtDetectorConstruction::DrawHitBox(int id) {

  // Int_t dispopt = PrtManager::Instance()->getDisplayOpt();
  // if (dispopt != 1 && dispopt != 2) return;

  std::cout << "reading occuhits.root " << std::endl;
  
  PrtEvent *event = new PrtEvent();
  TChain *tree;
  if (id == 1) tree = (TChain *)PrtManager::Instance()->getTree();
  else {
    tree = new TChain("data");
    tree->Add("occuhits.root");
  }

  TBranch *branch = tree->GetBranch("PrtEvent");
  if (!branch) return;
  branch->SetAddress(&event);

  Int_t nevent = tree->GetEntries();
  std::cout << id << " nevent  " << nevent << std::endl;

  int pix(0), mcp(0), prism(0);
  int npmt = fNRow * fNCol;
  int hload[12][6][6150] = {{{0}}};

  for (int i = 0; i < nevent; i++) {
    tree->GetEvent(i);
    for (auto hit : event->getHits()) {      
      pix = hit.getPixel();
      mcp = hit.getPmt();
      prism = hit.getPrism();
      if (prism < 0 || pix < 0 || pix > 6150 || mcp < 0) continue;
      hload[prism][mcp][pix]++;
    }
  }

  int maxload[16] = {0};
  for (int p = 0; p < fNBoxes; p++) {
    for (int im = 0; im < npmt; im++) {
      for (int i = 0; i < 6150; i++) {
        if (maxload[p] < hload[p][im][i]) maxload[p] = hload[p][im][i];
      }
    }
  }

  int rgbcolors[256][3] = {
    {0, 0, 130},     {0, 2, 131},     {0, 4, 132},     {0, 7, 134},     {0, 9, 135},
    {0, 11, 137},    {0, 14, 138},    {0, 16, 140},    {0, 18, 141},    {0, 21, 142},
    {0, 23, 144},    {0, 26, 145},    {0, 28, 147},    {0, 30, 148},    {0, 33, 150},
    {0, 35, 151},    {0, 37, 153},    {0, 40, 154},    {0, 42, 155},    {0, 45, 157},
    {0, 47, 158},    {0, 49, 160},    {0, 52, 161},    {0, 54, 163},    {0, 56, 164},
    {0, 59, 165},    {0, 61, 167},    {0, 64, 168},    {0, 66, 170},    {0, 68, 171},
    {0, 71, 173},    {0, 73, 174},    {0, 75, 176},    {0, 78, 177},    {0, 80, 178},
    {0, 83, 180},    {0, 85, 181},    {0, 87, 183},    {0, 90, 184},    {0, 92, 186},
    {0, 94, 187},    {0, 97, 188},    {0, 99, 190},    {0, 102, 191},   {0, 104, 193},
    {0, 106, 194},   {0, 109, 196},   {0, 111, 197},   {0, 113, 198},   {0, 116, 200},
    {0, 118, 201},   {0, 121, 203},   {0, 123, 204},   {0, 125, 206},   {0, 128, 207},
    {0, 130, 209},   {0, 132, 210},   {0, 135, 211},   {0, 137, 213},   {0, 140, 214},
    {0, 142, 216},   {0, 144, 217},   {0, 147, 219},   {0, 149, 220},   {0, 151, 221},
    {0, 154, 223},   {0, 156, 224},   {0, 159, 226},   {0, 161, 227},   {0, 163, 229},
    {0, 166, 230},   {0, 168, 232},   {0, 170, 233},   {0, 173, 234},   {0, 175, 236},
    {0, 178, 237},   {0, 180, 239},   {0, 182, 240},   {0, 185, 242},   {0, 187, 243},
    {0, 189, 244},   {0, 192, 246},   {0, 194, 247},   {0, 197, 249},   {0, 199, 250},
    {0, 201, 252},   {0, 204, 253},   {0, 206, 255},   {3, 207, 251},   {6, 207, 248},
    {9, 208, 245},   {12, 209, 241},  {16, 210, 238},  {19, 210, 235},  {22, 211, 232},
    {25, 212, 228},  {28, 212, 225},  {32, 213, 222},  {35, 214, 219},  {38, 214, 215},
    {41, 215, 212},  {45, 216, 209},  {48, 217, 206},  {51, 217, 202},  {54, 218, 199},
    {57, 219, 196},  {61, 219, 193},  {64, 220, 189},  {67, 221, 186},  {70, 221, 183},
    {73, 222, 180},  {77, 223, 176},  {80, 224, 173},  {83, 224, 170},  {86, 225, 167},
    {90, 226, 163},  {93, 226, 160},  {96, 227, 157},  {99, 228, 154},  {102, 229, 150},
    {106, 229, 147}, {109, 230, 144}, {112, 231, 141}, {115, 231, 137}, {118, 232, 134},
    {122, 233, 131}, {125, 233, 128}, {128, 234, 124}, {131, 235, 121}, {135, 236, 118},
    {138, 236, 115}, {141, 237, 111}, {144, 238, 108}, {147, 238, 105}, {151, 239, 102},
    {154, 240, 98},  {157, 240, 95},  {160, 241, 92},  {163, 242, 89},  {167, 243, 85},
    {170, 243, 82},  {173, 244, 79},  {176, 245, 76},  {180, 245, 72},  {183, 246, 69},
    {186, 247, 66},  {189, 247, 63},  {192, 248, 59},  {196, 249, 56},  {199, 250, 53},
    {202, 250, 50},  {205, 251, 46},  {208, 252, 43},  {212, 252, 40},  {215, 253, 37},
    {218, 254, 33},  {221, 255, 30},  {222, 251, 30},  {222, 248, 29},  {223, 244, 29},
    {224, 241, 28},  {224, 237, 28},  {225, 234, 27},  {225, 230, 26},  {226, 227, 26},
    {226, 223, 25},  {227, 220, 25},  {228, 216, 24},  {228, 213, 24},  {229, 210, 23},
    {229, 206, 23},  {230, 203, 22},  {230, 199, 22},  {231, 196, 21},  {231, 192, 21},
    {232, 189, 20},  {233, 185, 20},  {233, 182, 19},  {234, 178, 19},  {234, 175, 18},
    {235, 172, 18},  {235, 168, 17},  {236, 165, 17},  {237, 161, 16},  {237, 158, 16},
    {238, 154, 15},  {238, 151, 15},  {239, 147, 14},  {239, 144, 14},  {240, 140, 13},
    {240, 137, 12},  {241, 133, 12},  {242, 130, 11},  {242, 127, 11},  {243, 123, 10},
    {243, 120, 10},  {244, 116, 9},   {244, 113, 9},   {245, 109, 8},   {246, 106, 8},
    {246, 102, 7},   {247, 99, 7},    {247, 95, 6},    {248, 92, 6},    {248, 89, 5},
    {249, 85, 5},    {249, 82, 4},    {250, 78, 4},    {251, 75, 3},    {251, 71, 3},
    {252, 68, 2},    {252, 64, 2},    {253, 61, 1},    {253, 57, 1},    {254, 54, 0},
    {255, 51, 0},    {251, 49, 0},    {248, 48, 0},    {245, 47, 0},    {242, 46, 0},
    {239, 44, 0},    {236, 43, 0},    {233, 42, 0},    {230, 41, 0},    {227, 39, 0},
    {224, 38, 0},    {221, 37, 0},    {218, 36, 0},    {215, 34, 0},    {212, 33, 0},
    {209, 32, 0},    {206, 31, 0},    {203, 29, 0},    {200, 28, 0},    {197, 27, 0},
    {194, 26, 0},    {191, 24, 0},    {187, 23, 0},    {184, 22, 0},    {181, 21, 0},
    {178, 19, 0},    {175, 18, 0},    {172, 17, 0},    {169, 16, 0},    {166, 14, 0},
    {163, 13, 0},    {160, 12, 0},    {157, 11, 0},    {154, 9, 0},     {151, 8, 0},
    {148, 7, 0},     {145, 6, 0},     {142, 4, 0},     {139, 3, 0},     {136, 2, 0},
    {133, 1, 0}};

  G4VisAttributes *waDircHit = new G4VisAttributes(G4Colour(1., 1., 0.9, 0.2));
  waDircHit->SetVisibility(false);

  double tphi, dphi = 360 / fNBoxes * deg;
  G4LogicalVolume *lDircHit[fNBoxes];
  for (int i = 0; i < fNBoxes; i++) {
    tphi = dphi * i;
    double dx = fRadius * cos(tphi);
    double dy = fRadius  * sin(tphi);
    G4ThreeVector dirc =
      G4ThreeVector(dx, dy, 0.5 * fBar[2] * 4 + fPrizm[1] + fMcpActive[2] / 2. + fLens[2] + 630);

    G4Box *gDircHit = new G4Box("gDircHit", 400., 300., 10);
    lDircHit[i] = new G4LogicalVolume(gDircHit, defaultMaterial, Form("lDircHit%d", i), 0, 0, 0);
    lDircHit[i]->SetVisAttributes(waDircHit);

    G4RotationMatrix *tRot = new G4RotationMatrix();
    tRot->rotateZ(-tphi);
    new G4PVPlacement(tRot, dirc, lDircHit[i], "wDircHit", lExpHall, false, i);
  }

  G4LogicalVolume *lHit;
  G4VisAttributes *waHit;
  G4Box *gHit;
  for (int p = 0; p < 12; p++) {
    mcp = 0;
    double gapx = (fPrizm[2] - fNCol * fMcpTotal[0]) / (double)(fNCol + 1);
    double gapy = (fPrizm[0] - fNRow * fMcpTotal[1]) / (double)(fNRow + 1);
    for (int i = 0; i < fNCol; i++) {
      for (int j = 0; j < fNRow; j++) {
        double shiftx = i * (fMcpTotal[0] + gapx) - fPrizm[3] / 2. + fMcpTotal[0] / 2 + gapx;
        double shifty = j * (fMcpTotal[1] + gapy) - fPrizm[0] / 2. + fMcpTotal[1] / 2 + gapy;
        int pixelId = 0;

	for (int i1 = 0; i1 < fNpix1; i1++) {
	  for (int j1 = 0; j1 < fNpix2; j1++) {
	    
            double colorid = hload[p][mcp][pixelId] / (double)maxload[p];
            if (colorid > 1) colorid = 1;
            double hight = colorid * 5;
            int cid = colorid * 255;
            if (cid > 1 && maxload[p] > 10) {
              gHit = new G4Box("gHit", 0.5 * fMcpActive[0] / fNpix1, 0.5 * fMcpActive[1] / fNpix2,
                               hight + 0.2);
              lHit = new G4LogicalVolume(gHit, BarMaterial, "lHit", 0, 0, 0);
              waHit = new G4VisAttributes(G4Colour(
                rgbcolors[cid][0] / 255., rgbcolors[cid][1] / 255., rgbcolors[cid][2] / 255., 1.0));
              waHit->SetForceSolid(true);
              lHit->SetVisAttributes(waHit);

              double shiftx1 = shiftx + i1 * (fMcpActive[0] / fNpix1) - fMcpActive[0] / 2. +
                               0.5 * fMcpActive[0] / fNpix1;
              double shifty1 = shifty + j1 * (fMcpActive[1] / fNpix2) - fMcpActive[1] / 2. +
                               0.5 * fMcpActive[1] / fNpix2;
              auto vpix = G4ThreeVector(shiftx1, shifty1, hight);
              new G4PVPlacement(0, vpix, lHit, "wDircHit", lDircHit[p], false, pixelId);
            }
            pixelId++;
          }
        }
        mcp++;
      }
    }
  }

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
}

void PrtDetectorConstruction::SetLens(int id) {
  
  // if(id==0){
  //   fPrismShift.setZ(fPrismShift.z()-fLens[2]);
  // }
  // G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void PrtDetectorConstruction::SetQuantumEfficiency(int id) {
  const int num = 36;
  // ideal pmt quantum efficiency
  double QuantumEfficiencyIdial[num] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  // Burle PMT's
  double QuantumEfficiencyB[num] = {0.,   0.001, 0.002, 0.005, 0.01,  0.015, 0.02,  0.03, 0.04,
                                    0.05, 0.06,  0.07,  0.09,  0.1,   0.13,  0.15,  0.17, 0.2,
                                    0.24, 0.26,  0.28,  0.282, 0.284, 0.286, 0.288, 0.29, 0.28,
                                    0.26, 0.24,  0.22,  0.20,  0.18,  0.15,  0.13,  0.12, 0.10};

  // hamamatsu pmt quantum efficiency
  double QuantumEfficiencyPMT[num] = {
    0.001, 0.002, 0.004, 0.007, 0.011, 0.015, 0.020, 0.026, 0.033, 0.040, 0.045, 0.056,
    0.067, 0.085, 0.109, 0.129, 0.138, 0.147, 0.158, 0.170, 0.181, 0.188, 0.196, 0.203,
    0.206, 0.212, 0.218, 0.219, 0.225, 0.230, 0.228, 0.222, 0.217, 0.210, 0.199, 0.177};

  if (id == 0) fQuantumEfficiency = QuantumEfficiencyIdial;
  if (id == 1) fQuantumEfficiency = QuantumEfficiencyPMT;
  if (id == 2) fQuantumEfficiency = QuantumEfficiencyB;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
