
#include "PrtDetectorConstruction.h"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"


#include "PrtManager.h"
#include "PrtPrizmSD.h"
#include "PrtPixelSD.h"

PrtDetectorConstruction::PrtDetectorConstruction()
  : G4VUserDetectorConstruction(){

  fGeomId = PrtManager::Instance()->GetGeometry();
  fMcpLayout = PrtManager::Instance()->GetMcpLayout();

  //fGeomId=2014; // 2012

  fNRow = 6;
  fNCol = 4;
  fNBar = PrtManager::Instance()->GetRadiator();
  if(fNBar<0) fNBar=1;

  fHall[0] = 2000; fHall[1] = 2000; fHall[2] = 4000;
  //fBar[0] = 17; fBar[1] = 32; fBar[2] = 4200;
  //fMirror[0] = 20; fMirror[1] = 40; fMirror[2] =1;
  //fPrizm[0] = 170; fPrizm[1] = 300; fPrizm[2] = 30+300*tan(37*deg); fPrizm[3] = 30;
  fPrizm[0] = 390; fPrizm[1] = 300; fPrizm[3] = 50; fPrizm[2]= fPrizm[3]+300*tan(38*deg); 
  fMirror[0] = 20; fMirror[1] = fPrizm[0]; fMirror[2] =1;
  //  fPrizm[0] = 170; fPrizm[1] = 300; fPrizm[2] = 50+300*tan(45*deg); fPrizm[3] = 50;
  fBar[0] = 17; fBar[1] = fPrizm[0]/fNBar; fBar[2] = 4200; //4200

  fMcpTotal[0] = fMcpTotal[1] = 53+4; fMcpTotal[2]=1;
  fMcpActive[0] = fMcpActive[1] = 53; fMcpActive[2]=1;
  fLens[0] = fLens[1] = 40; fLens[2]=10;

  if(PrtManager::Instance()->GetLens() == 2){
    fLens[0] = 50; fLens[1] = 175; fLens[2]=14.4;
  }
  if(PrtManager::Instance()->GetLens() == 3){
    fLens[0] = 30; fLens[1] = fBar[1]; fLens[2]=15;
  }

  PrtManager::Instance()->SetRadiatorL(fBar[2]);
  PrtManager::Instance()->SetRadiatorW(fBar[1]);
  PrtManager::Instance()->SetRadiatorH(fBar[0]);			  
  fPrtRot = new G4RotationMatrix();
  //create a messenger for this class
  fGeomMessenger = new PrtDetectorConstructionMessenger(this);
}

PrtDetectorConstruction::~PrtDetectorConstruction(){}

G4VPhysicalVolume* PrtDetectorConstruction::Construct(){
  DefineMaterials();

  // ------------- Volumes --------------

  // The experimental Hall
  G4Box* gExpHall = new G4Box("gExpHall",fHall[0],fHall[1],fHall[2]);
  lExpHall = new G4LogicalVolume(gExpHall,defaultMaterial,"lExpHall",0,0,0);
  G4VPhysicalVolume* wExpHall  = new G4PVPlacement(0,G4ThreeVector(),lExpHall,"gExpHall",0,false,0);

  // The front material
  G4Box* gFront = new G4Box("gFront",200.,200.,5);
  lFront = new G4LogicalVolume(gFront,frontMaterial,"lFront",0,0,0);
  
  if(PrtManager::Instance()->GetGeometry() == 3){
    // new G4PVPlacement(0,G4ThreeVector(0,0,-1200),lFront,"wFront",lExpHall,false,0);
  }
  
  // The DIRC
  G4Box* gDirc = new G4Box("gDirc",400.,300.,fBar[2]/2.+350);
  lDirc = new G4LogicalVolume(gDirc,defaultMaterial,"lDirc",0,0,0);

  G4int fNBoxes = (PrtManager::Instance()->GetGeometry() == 0)? 16 :1;
  G4double radius = 1000;
  G4double tphi, dphi = 22.5*deg; //22.5*deg;
  G4int BoxId = 0;

  if(PrtManager::Instance()->GetRunType()==1){ 
    // LUT
    new G4PVPlacement(0,G4ThreeVector(0,0,0),lDirc,"wDirc",lExpHall,false,0);
  }else{ 
    for(int i=0; i<fNBoxes; i++){
      tphi = dphi*i; 
      G4double dx = radius * cos(tphi);
      G4double dy = radius * sin(tphi);

      G4RotationMatrix *tRot = new G4RotationMatrix();
      tRot->rotateZ(-tphi);     
      new G4PVPlacement(tRot,G4ThreeVector(dx,dy,0),lDirc,"wDirc",lExpHall,false,i); //G4ThreeVector(dx,dy,-500)
    }
  }


  // The Bar
  G4Box *gBar = new G4Box("gBar",fBar[0]/2.,fBar[1]/2.,fBar[2]/2.);
  lBar = new G4LogicalVolume(gBar,BarMaterial,"lBar",0,0,0);
  //wBar =  new G4PVPlacement(0,G4ThreeVector(0,0,0),lBar,"wBar", lDirc,false,0);

  bool plate = false;
  if(fNBar==1){
    wBar =  new G4PVPlacement(0,G4ThreeVector(0,0,0),lBar,"wBar", lDirc,false,0);
  }else{
    int pixelId = 0;
    for(int i=0; i<fNBar; i++){
      double shifty = i*(fBar[1]+0.01)-fPrizm[0]/2. + fBar[1]/2.;
      pDirc[i] = new G4PVPlacement(0,G4ThreeVector(0,shifty,0),lBar,"wBar",lDirc,false,i);      
    }
  }
  

  // The Mirror
  G4Box* gMirror = new G4Box("gMirror",fMirror[0]/2.,fMirror[1]/2.,fMirror[2]/2.);
  lMirror = new G4LogicalVolume(gMirror,MirrorMaterial,"lMirror",0,0,0);
  wMirror =new G4PVPlacement(0,G4ThreeVector(0,0,-fBar[2]/2.-fMirror[2]/2.),lMirror,"wMirror", lDirc,false,0);

  // The Lens
  G4double lensrad = 70;
  G4double lensMinThikness = 2;
  G4Box* gfbox = new G4Box("Fbox",fLens[0]/2.,fLens[1]/2.,fLens[2]/2.);
  G4ThreeVector zTrans(0, 0, -lensrad+fLens[2]/2.-lensMinThikness);

  G4Sphere* gsphere = new G4Sphere("Sphere",0,70,0.*deg,360.*deg,0.*deg,380.*deg);
  G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Fbox*Sphere", gfbox, gsphere,new G4RotationMatrix(),zTrans); 
  G4SubtractionSolid* gLens2 = new G4SubtractionSolid("Fbox-Sphere", gfbox, gsphere,new G4RotationMatrix(),zTrans);

  lLens1 = new G4LogicalVolume(gLens1,Nlak33aMaterial,"lLens1",0,0,0); //Nlak33aMaterial  
  lLens2 = new G4LogicalVolume(gLens2,BarMaterial,"lLens2",0,0,0);

  if(PrtManager::Instance()->GetLens() == 3){ // 3-component spherical lens
    G4double lensMinThikness = 2; 
  
    G4double r1 = 0; //PrtManager::Instance()->GetTest1();
    G4double r2 = 0; //PrtManager::Instance()->GetTest2();
  
    G4double lensrad1 = (r1==0)? 47.8: r1;
    G4double lensrad2 = (r2==0)? 29.1: r2;
    
 
    G4ThreeVector zTrans1(0, 0, -lensrad1-fLens[2]/2.+lensrad1-sqrt(lensrad1*lensrad1-fLens[0]/2.*fLens[0]/2.)+lensMinThikness);
    G4ThreeVector zTrans2(0, 0, -lensrad2+fLens[2]/2.-lensMinThikness);
    
    G4Sphere* gsphere1 = new G4Sphere("Sphere1",0,lensrad1,0,360*deg,0,360*deg);
    G4Sphere* gsphere2 = new G4Sphere("Sphere2",0,lensrad2,0,360*deg,0,360*deg);


    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Fbox*Sphere1", gfbox, gsphere1,new G4RotationMatrix(),zTrans1); 
    G4SubtractionSolid* gLenst = new G4SubtractionSolid("Fbox-Sphere1", gfbox, gsphere1, new G4RotationMatrix(),zTrans1);

    G4IntersectionSolid* gLens2 = new G4IntersectionSolid("gLenst*Sphere2", gLenst, gsphere2, new G4RotationMatrix(),zTrans2);
    G4SubtractionSolid* gLens3 = new G4SubtractionSolid("gLenst-Sphere2", gLenst, gsphere2,new G4RotationMatrix(),zTrans2);
    
    lLens1 = new G4LogicalVolume(gLens1,BarMaterial,"lLens1",0,0,0);
    lLens2 = new G4LogicalVolume(gLens2,Nlak33aMaterial,"lLens2",0,0,0);
    lLens3 = new G4LogicalVolume(gLens3,BarMaterial,"lLens3",0,0,0);
  }

  if(PrtManager::Instance()->GetLens() != 0 && PrtManager::Instance()->GetLens() != 10){
    new G4PVPlacement(0,G4ThreeVector(0,0,fBar[2]/2.+fLens[2]/2.),lLens1,"wLens1", lDirc,false,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,fBar[2]/2.+fLens[2]/2.),lLens2,"wLens2", lDirc,false,0);
    if(PrtManager::Instance()->GetLens() == 3)  new G4PVPlacement(0,G4ThreeVector(0,0,fBar[2]/2.+fLens[2]/2.),lLens3,"wLens3", lDirc,false,0);
  }else{
    fLens[2]=0; 
  }

  // The Prizm
  G4Trap* gPrizm = new G4Trap("gPrizm",fPrizm[0],fPrizm[1],fPrizm[2],fPrizm[3]);
  lPrizm = new G4LogicalVolume(gPrizm, BarMaterial,"lPrizm",0,0,0);
  G4RotationMatrix* xRot = new G4RotationMatrix();
  xRot->rotateX(M_PI/2.*rad);
  fPrismShift = G4ThreeVector((fPrizm[2]+fPrizm[3])/4.-fPrizm[3]/2.,0,fBar[2]/2.+fPrizm[1]/2.+fLens[2]);
  new G4PVPlacement(xRot,fPrismShift,lPrizm,"wPrizm", lDirc,false,0);

  if(fMcpLayout==1){
    // standard mcp pmt layout
    // The MCP
    G4Box* gMcp = new G4Box("gMcp",fMcpTotal[0]/2.,fMcpTotal[1]/2.,fMcpTotal[2]/2.);
    lMcp = new G4LogicalVolume(gMcp,BarMaterial,"lMcp",0,0,0);
    
    // The MCP Pixel
    G4Box* gPixel = new G4Box("gPixel",fMcpActive[0]/16.,fMcpActive[1]/16.,fMcpActive[2]/16.);
    lPixel = new G4LogicalVolume(gPixel,BarMaterial,"lPixel",0,0,0);
    
    int pixelId = 0;
    for(int i=0; i<8; i++){
      for(int j=0; j<8; j++){
	double shiftx = i*(fMcpActive[0]/8.)-fMcpActive[0]/2.+fMcpActive[0]/16.;
	double shifty = j*(fMcpActive[0]/8.)-fMcpActive[0]/2.+fMcpActive[0]/16.;
	new G4PVPlacement(0,G4ThreeVector(shiftx,shifty,0),lPixel,"wPixel", lMcp,false,pixelId++);      
      }
    }
 
    int mcpId = 0;
    G4double gapx = (fPrizm[2]-4*fMcpTotal[0])/5.;
    G4double gapy = (fPrizm[0]-6*fMcpTotal[0])/7.;
    for(int i=0; i<fNCol; i++){
      for(int j=0; j<fNRow; j++){
	double shiftx = i*(fMcpTotal[0]+gapx)-fPrizm[3]/2.+fMcpTotal[0]/2+gapx;
	double shifty = j*(fMcpTotal[0]+gapy)-fPrizm[0]/2.+fMcpTotal[0]/2+gapy;
	new G4PVPlacement(0,G4ThreeVector(shiftx,shifty,fBar[2]/2.+fPrizm[1]+fMcpActive[2]/2.+fLens[2]),lMcp,"wMcp", lDirc,false,mcpId++);
      }
    }

  }
  if(fMcpLayout==0){
    // only mcps
    // The MCP
    G4Box* gMcp = new G4Box("gMcp",fPrizm[2]/2.,fPrizm[0]/2.,fMcpTotal[2]/2.);
    lMcp = new G4LogicalVolume(gMcp,BarMaterial,"lMcp",0,0,0);
  
    // The MCP Pixel
    fNpix1=4;
    fNpix2=6;
    G4Box* gPixel = new G4Box("gPixel",fMcpTotal[0]/2.,fMcpTotal[1]/2.,fMcpTotal[2]/16.);
    lPixel = new G4LogicalVolume(gPixel,BarMaterial,"lPixel",0,0,0);

    double disx = (fPrizm[0]-fNpix2*fMcpTotal[0])/(double)fNpix2;
    double disy = (fPrizm[1]-fNpix1*fMcpTotal[0])/(double)(fNpix1+1);
    if(true){
      int pixelId=0;
      for(int i=0; i<fNpix1; i++){
	for(int j=0; j<fNpix2; j++){
	  double shiftx = i*(fMcpTotal[0] +disy/2.) - fPrizm[2]/2. + fMcpTotal[0]/2.+disy/2.;
	  double shifty = j*(fMcpTotal[0] +disy/2.) - fPrizm[0]/2. + fMcpTotal[0]/2.+disy/2.;
	  new G4PVPlacement(0,G4ThreeVector(shiftx,shifty,0),lPixel,"wPixel", lMcp,false,pixelId++);      
	}
      }
    }else{
      new G4PVPlacement(0,G4ThreeVector(0,0,0),lPixel,"wPixel", lMcp,false,1);
    } 
    new G4PVPlacement(0,G4ThreeVector(fPrizm[2]/2.-fPrizm[3]/2.,0,fBar[2]/2.+fPrizm[1]+fMcpActive[2]/2.+fLens[2]),lMcp,"wMcp", lDirc,false,1);
  }
  if(fMcpLayout==3){
    // one mcp layout
    G4Box* gMcp = new G4Box("gMcp",0.5*fPrizm[2],0.5*fPrizm[0],0.5*fMcpTotal[2]);
    lMcp = new G4LogicalVolume(gMcp,BarMaterial,"lMcp",0,0,0);
  

    G4double pixSize = 3*mm;
    fNpix1 = fPrizm[2]/pixSize-1;
    fNpix2 = fPrizm[0]/pixSize-1;
    std::cout<<"fNpix1  "<< fNpix1<<"   "<<fPrizm[0]<<std::endl;
    std::cout<<"fNpix2  "<< fNpix2 << "  "<< fPrizm[2]<<std::endl;
 
    // fNpix1 = 90;  
    // fNpix2 = 120;
    // fNpix1 = 4;
    // fNpix2 = 7;
    G4Box* gPixel = new G4Box("gPixel",0.5*fPrizm[2]/fNpix1,0.5*fPrizm[0]/fNpix2,0.2);
    lPixel = new G4LogicalVolume(gPixel,BarMaterial,"lPixel",0,0,0);
    int pixelId=0;
    for(int i=0; i<fNpix1; i++){
      for(int j=0; j<fNpix2; j++){
	//std::cout<<"i  "<<i <<"   j  "<< j <<std::endl;
	
	double shiftx = i*(fPrizm[2]/fNpix1) - 0.5*fPrizm[2]+0.5*fPrizm[2]/fNpix1;
	double shifty = j*(fPrizm[0]/fNpix2) - 0.5*fPrizm[0]+0.5*fPrizm[0]/fNpix2;
	new G4PVPlacement(0,G4ThreeVector(shiftx,shifty,0),lPixel,"wPixel", lMcp,false,pixelId++);      
      }
    }
    new G4PVPlacement(0,G4ThreeVector(fPrizm[2]/2.-fPrizm[3]/2.,0,fBar[2]/2.+fPrizm[1]+fMcpActive[2]/2.+fLens[2]),lMcp,"wMcp", lDirc,false,1);
  }

  const G4int num = 36; 
  G4double WaveLength[num];
  G4double PhotonEnergy[num];
  G4double PMTReflectivity[num];
  G4double EfficiencyMirrors[num];
  const G4double LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;
  for(int i=0;i<num;i++){
    WaveLength[i]= (300 +i*10)*nanometer;
    PhotonEnergy[num-(i+1)]= LambdaE/WaveLength[i];
    PMTReflectivity[i]=0.;
    EfficiencyMirrors[i]=0; 
  }

  /***************** QUANTUM EFFICIENCY OF BURLE AND HAMAMTSU PMT'S *****/

  //ideal pmt quantum efficiency
  G4double QuantumEfficiencyIdial[num]=
    {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     1.0,1.0,1.0,1.0,1.0,1.0};

  // Burle PMT's 
  G4double QuantumEfficiencyB[num] =
    {0.,0.001,0.002,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.06,
     0.07,0.09,0.1,0.13,0.15,0.17,0.2,0.24,0.26,0.28,0.282,0.284,0.286,
     0.288,0.29,0.28,0.26,0.24,0.22,0.20,0.18,0.15,0.13,0.12,0.10};
  
  //hamamatsu pmt quantum efficiency
  G4double QuantumEfficiencyPMT[num]=
    {0.001,0.002,0.004,0.007,0.011,0.015,0.020,0.026,0.033,0.040,0.045,
     0.056,0.067,0.085,0.109,0.129,0.138,0.147,0.158,0.170,
     0.181,0.188,0.196,0.203,0.206,0.212,0.218,0.219,0.225,0.230,
     0.228,0.222,0.217,0.210,0.199,0.177};

  // these quantum efficiencies have to be multiplied by geometry
  //   efficiency of given PMT's
  //   for Hamamatsu by factor 0.7
  //   for Burle by factor 0.45 
  for(G4int k=0;k<36;k++){
      QuantumEfficiencyB[k] =  QuantumEfficiencyB[k] * 0.45 ;
      QuantumEfficiencyPMT[k] =  QuantumEfficiencyPMT[k] *.7;
    }
 
  // G4double QuantumEfficiency[num]= 
  //    { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
  //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
  //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

  //  G4double QuantumEfficiencyPMT[num] =
  //    { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
  //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
  //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
  
  /* define quantum efficiency for burle PMT's - the same efficiency is 
     assigned to pads and also to slots !!!! */
    
  //burle pmt - bigger slots => logicPad
  G4MaterialPropertiesTable* PhotocatodBurleMPT = new G4MaterialPropertiesTable();
  PhotocatodBurleMPT->AddProperty("EFFICIENCY",  PhotonEnergy,QuantumEfficiencyB,num);
  PhotocatodBurleMPT->AddProperty("REFLECTIVITY",PhotonEnergy,PMTReflectivity,num);

 
  G4OpticalSurface* BurlePMTOpSurface= 
    new G4OpticalSurface("BurlePMTOpSurface",glisur,polished,
			 dielectric_metal);
  BurlePMTOpSurface->SetMaterialPropertiesTable(PhotocatodBurleMPT);

  // // assignment for pad
  // if(burle)
  //   new G4LogicalSkinSurface("BurlePMTSurface",logicBurPad,BurlePMTOpSurface); 

  // if(burle1)
  //   new G4LogicalSkinSurface("Burle1PMTSurface",logicBur1Pad,BurlePMTOpSurface); 

  /* hamamatsu pmt's - smaller slots => quantum efficiency again
     assign to slot and pad */
  
  fQuantumEfficiency = QuantumEfficiencyIdial;
  G4MaterialPropertiesTable* PhotocatodHamamatsuMPT = new G4MaterialPropertiesTable();
  PhotocatodHamamatsuMPT->AddProperty("EFFICIENCY",  PhotonEnergy,fQuantumEfficiency,num);
  PhotocatodHamamatsuMPT->AddProperty("REFLECTIVITY",PhotonEnergy,PMTReflectivity,num);

  G4OpticalSurface* HamamatsuPMTOpSurface= 
    new G4OpticalSurface("HamamatsuPMTOpSurface",glisur,polished,dielectric_metal);
  HamamatsuPMTOpSurface->SetMaterialPropertiesTable(PhotocatodHamamatsuMPT);

  // // assignment to pad
  // if(hamamatsu8500)
  new G4LogicalSkinSurface("HamamatsuPMTSurface",lPixel,HamamatsuPMTOpSurface);

  // Mirror
  G4OpticalSurface* MirrorOpSurface= 
    new G4OpticalSurface("MirrorOpSurface",glisur,polished,dielectric_metal);
  
  G4double ReflectivityMirrorBar[num]={
    0.8755,0.882,0.889,0.895,0.9,0.9025,0.91,0.913,0.9165,0.92,0.923,
    0.9245,0.9285,0.932,0.934,0.935,0.936,0.9385,0.9395,0.94,
    0.9405,0.9405,0.9405,0.9405,0.94,0.9385,0.936,0.934,
    0.931,0.9295,0.928,0.928,0.921,0.92,0.927,0.9215};

  G4MaterialPropertiesTable *MirrorMPT = new G4MaterialPropertiesTable();
  MirrorMPT->AddProperty("REFLECTIVITY", PhotonEnergy, ReflectivityMirrorBar, num);
  MirrorMPT->AddProperty("EFFICIENCY", PhotonEnergy, EfficiencyMirrors,   num);
  
  MirrorOpSurface->SetMaterialPropertiesTable(MirrorMPT);
  new G4LogicalSkinSurface("MirrorSurface", lMirror,MirrorOpSurface);

  SetVisualization();

  return wExpHall;
}

void PrtDetectorConstruction::DefineMaterials(){
  G4String symbol;             //a=mass of a mole;
  G4double a, z,  density;      //z=mean number of protons;  

  G4int ncomponents, natoms;
  G4double fractionmass;

  // define Elements
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* Si = new G4Element("Silicon" ,symbol="Si", z= 14., a= 28.09*g/mole);
  G4Element* Al = new G4Element("Aluminum",symbol="Al", z=13., a=26.98*g/mole);


  // quartz material = SiO2
  G4Material* SiO2 = new G4Material("quartz",density= 2.200*g/cm3, ncomponents=2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O , natoms=2);

  Nlak33aMaterial  = new G4Material("Nlak33a",density= 4.220*g/cm3, ncomponents=2);
  Nlak33aMaterial->AddElement(Si, natoms=1);
  Nlak33aMaterial->AddElement(O , natoms=2);

  G4Material* Vacuum = new G4Material("interGalactic", 1., 1.008*g/mole, 
				      1.e-25*g/cm3, kStateGas, 
				      2.73*kelvin, 3.e-18*pascal);
  G4Material* Air = new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material* Aluminum = new G4Material("Aluminum",density=2.7*g/cm3,ncomponents=1);
  Aluminum->AddElement(Al,fractionmass=1.0);

  G4Material* KamLandOil = new G4Material("KamLandOil",density=0.914*g/cm3,ncomponents=2);
  KamLandOil->AddElement(C,natoms=12);
  KamLandOil->AddElement(H,natoms=26);

  G4Material* CarbonFiber = new G4Material("CarbonFiber", density=0.145*g/cm3, ncomponents=1);
  CarbonFiber->AddElement(C,natoms=1);
			

  /* as I don't know the exact material composition,
     I will use Epoxyd material composition and add
     the optical property of Epotek to this material */

  G4Material* Epotek = new G4Material("Epotek",density=1.2*g/cm3,ncomponents=3);

  Epotek->AddElement(C,natoms=3);
  Epotek->AddElement(H,natoms=5);
  Epotek->AddElement(O,natoms=2);


  // assign main materials
  if(PrtManager::Instance()->GetGeometry() < 10) defaultMaterial = Vacuum;
  else defaultMaterial = Air; //Vacuum // material of world
  frontMaterial =  CarbonFiber; 
  BarMaterial = SiO2; // material of all Bars, Quartz and Window
  OilMaterial = KamLandOil; // material of volume 1,2,3,4
  MirrorMaterial = Aluminum; // mirror material
  epotekMaterial = Epotek; // Epotek material - glue between bars


  // ------------ Generate & Add Material Properties Table ------------

  static const G4double LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;
  const G4int num = 36;
  G4double WaveLength[num];
  G4double Absorption[num]; // default value for absorption
  G4double AirAbsorption[num]; // absorption value for air
  G4double AirRefractiveIndex[num]; // air refractive index
  G4double PhotonEnergy[num]; // energy of photons which correspond to the given 
  // refractive or absoprtion values

  G4double PhotonEnergyNlak33a[76] = {1,1.2511,1.26386,1.27687,1.29016,1.30372,1.31758,1.33173,1.34619,1.36097,1.37607,1.39152,1.40731,1.42347,1.44,1.45692,1.47425,1.49199,1.51016,1.52878,1.54787,1.56744,1.58751,1.6081,1.62923,1.65092,1.6732,1.69609,1.71961,1.7438,1.76868,1.79427,1.82062,1.84775,1.87571,1.90452,1.93423,1.96488,1.99652,2.0292,2.06296,2.09787,2.13398,2.17135,2.21006,2.25017,2.29176,2.33492,2.37973,2.42631,2.47473,2.52514,2.57763,2.63236,2.68946,2.7491,2.81143,2.87666,2.94499,3.01665,3.09187,3.17095,3.25418,3.34189,3.43446,3.53231,3.6359,3.74575,3.86244,3.98663,4.11908,4.26062,4.41225,4.57506,4.75035,4.93961};

  /*************************** ABSORPTION COEFFICIENTS *****************************/

  // absorption of KamLandOil per 50 cm - from jjv
  G4double KamLandOilAbsorption[num]=
    {0.97469022,0.976603956,0.978511548,0.980400538,0.982258449,0.984072792,
     0.985831062,0.987520743,0.989129303,0.990644203,0.992052894,
     0.993342822,0.994501428,0.995516151,0.996374433,0.997063719,
     0.997571464,0.997885132,0.997992205,0.997880183,0.997536591,
     0.99,0.98,0.97,0.96,0.94,0.93,0.924507,0.89982,0.883299,
     0.85657,0.842637,0.77020213,0.65727,0.324022,0.019192};

  // absorption of quartz per 1 m - from jjv
  G4double QuartzAbsorption[num] = 
    {0.999572036,0.999544661,0.999515062,0.999483019,0.999448285,
     0.999410586,0.999369611,0.999325013,0.999276402,0.999223336,
     0.999165317,0.999101778,0.999032079,0.998955488,0.998871172,
     0.998778177,0.99867541,0.998561611,0.998435332,0.998294892,0.998138345,
     0.997963425,0.997767484,0.997547418,
     0.99729958,0.99701966,0.99670255,0.996342167,0.995931242,0.995461041,
     0.994921022,0.994298396,0.993577567,0.992739402,0.991760297,0.990610945};
  
  // absorption of epotek per one layer - thicknes 0.001'' - from jjv
  G4double EpotekAbsorption[num] = 
    {0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.9999,0.9998,0.9995,0.999,0.998,0.997,0.996,0.9955,0.993,
     0.9871,0.9745};

  //N-Lak 33a
  G4double Nlak33aAbsorption[76]={371813,352095,331021,310814,291458,272937,255238,238342,222234,206897,192313,178463,165331,152896,141140,130043,119585,109747,100507,91846.3,83743.1,76176.7,69126.1,62570.2,56488,50858.3,45660.1,40872.4,36474.6,32445.8,28765.9,25414.6,22372.2,19619.3,17136.9,14906.5,12910.2,11130.3,9550.13,8153.3,6924.25,5848.04,4910.46,4098.04,3398.06,2798.54,2288.32,1856.99,1494.92,1193.28,943.973,739.657,573.715,440.228,333.94,250.229,185.064,134.967,96.9664,68.5529,47.6343,32.4882,21.7174,14.2056,9.07612,5.65267,3.4241,2.01226,1.14403,0.62722,0.330414,0.166558,0.0799649,0.0363677,0.0155708,0.00623089};
  

  G4double EpotekThickness = 0.001*2.54*cm;
  for(int i=0;i<num;i++){
    WaveLength[i]= (300 +i*10)*nanometer;
    Absorption[i]= 100*m; // not true, just due to definiton -> not absorb any
    AirAbsorption[i] = 4.*cm; // if photon in the air -> kill it immediately
    AirRefractiveIndex[i] = 1.; 
    PhotonEnergy[num-(i+1)]= LambdaE/WaveLength[i];

    /* as the absorption is given per length and G4 needs 
       mean free path length, calculate it here
       mean free path length - taken as probability equal 1/e
       that the photon will be absorbed */
      
    EpotekAbsorption[i] = (-1)/log(EpotekAbsorption[i])*EpotekThickness;
    QuartzAbsorption[i] = (-1)/log(QuartzAbsorption[i])*100*cm;
    KamLandOilAbsorption[i] = (-1)/log(KamLandOilAbsorption[i])*50*cm;
  }

  /**************************** REFRACTIVE INDEXES ****************************/
  
  // only phase refractive indexes are necessary -> g4 calculates group itself !!
  
  G4double QuartzRefractiveIndex[num]={
    1.456535,1.456812,1.4571,1.457399,1.457712,1.458038,1.458378,
    1.458735,1.459108,1.4595,1.459911,1.460344,1.460799,1.46128,
    1.461789,1.462326,1.462897,1.463502,1.464146,1.464833,
    1.465566,1.46635,1.46719,1.468094,1.469066,1.470116,1.471252,1.472485,
    1.473826,1.475289,1.476891,1.478651,1.480592,1.482739,1.485127,1.487793};

  G4double EpotekRefractiveIndex[num]={
    1.554034,1.555575,1.55698,1.558266,1.559454,1.56056,1.561604,
     1.562604,1.563579,1.564547,1.565526,1.566536,1.567595,
     1.568721,1.569933,1.57125,1.57269,1.574271,1.576012,
     1.577932,1.580049,1.582381,1.584948,1.587768,1.590859,
     1.59424,1.597929,1.601946,1.606307,1.611033,1.616141,1.621651,1.62758,
     1.633947,1.640771,1.64807};

  G4double KamLandOilRefractiveIndex[num]={
    1.433055,1.433369,1.433698,1.434045,1.434409,1.434793,1.435198,
    1.435626,1.436077,1.436555,1.4371,1.4376,1.4382,1.4388,1.4395,
    1.4402,1.4409,1.4415,1.4425,1.4434,1.4444,1.4455,1.4464,1.4479,1.4501,
    1.450428,1.451976,1.453666,1.455513,1.45754,1.45977,1.462231,1.464958,
    1.467991,1.471377,1.475174};

  double Nlak33aRefractiveIndex[76]={1.73816,1.73836,1.73858,1.73881,1.73904,1.73928,1.73952,1.73976,1.74001,1.74026,1.74052,1.74078,1.74105,1.74132,1.7416,1.74189,1.74218,1.74249,1.74279,1.74311,1.74344,1.74378,1.74412,1.74448,1.74485,1.74522,1.74562,1.74602,1.74644,1.74687,1.74732,1.74779,1.74827,1.74878,1.7493,1.74985,1.75042,1.75101,1.75163,1.75228,1.75296,1.75368,1.75443,1.75521,1.75604,1.75692,1.75784,1.75882,1.75985,1.76095,1.76211,1.76335,1.76467,1.76608,1.76758,1.7692,1.77093,1.77279,1.7748,1.77698,1.77934,1.7819,1.7847,1.78775,1.79111,1.79481,1.79889,1.80343,1.8085,1.81419,1.82061,1.8279,1.83625,1.84589,1.85713,1.87039};

  /* ASSIGNING REFRACTIVE AND ABSORPTION PROPERTIES TO THE GIVEN MATERIALS */

  // Quartz material => Si02
  G4MaterialPropertiesTable* QuartzMPT = new G4MaterialPropertiesTable();
  QuartzMPT->AddProperty("RINDEX",       PhotonEnergy, QuartzRefractiveIndex,num);
  QuartzMPT->AddProperty("ABSLENGTH",    PhotonEnergy, QuartzAbsorption,           num);

  // assign this parameter table to BAR material
  BarMaterial->SetMaterialPropertiesTable(QuartzMPT);

  // Air
  G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
  AirMPT->AddProperty("RINDEX",    PhotonEnergy, AirRefractiveIndex, num);
  AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAbsorption,      num);
  //  assign this parameter table to the air 
  defaultMaterial->SetMaterialPropertiesTable(AirMPT);


  // KamLandOil                                                
  G4MaterialPropertiesTable* KamLandOilMPT = new G4MaterialPropertiesTable();
  KamLandOilMPT->AddProperty("RINDEX", PhotonEnergy, KamLandOilRefractiveIndex, num);
  KamLandOilMPT->AddProperty("ABSLENGTH", PhotonEnergy, KamLandOilAbsorption, num);
  // assing this parameter table  to the KamLandOil
  OilMaterial->SetMaterialPropertiesTable(KamLandOilMPT);  

  // N-Lak 33a                                                
  G4MaterialPropertiesTable* Nlak33aMPT = new G4MaterialPropertiesTable();
  Nlak33aMPT->AddProperty("RINDEX", PhotonEnergyNlak33a, Nlak33aRefractiveIndex, 76);
  Nlak33aMPT->AddProperty("ABSLENGTH",PhotonEnergyNlak33a, Nlak33aAbsorption, 76);
  Nlak33aMaterial->SetMaterialPropertiesTable(Nlak33aMPT);  

  // Epotek Glue                                        
  G4MaterialPropertiesTable* EpotekMPT = new G4MaterialPropertiesTable();
  EpotekMPT->AddProperty("RINDEX", PhotonEnergy, EpotekRefractiveIndex, num);
  EpotekMPT->AddProperty("ABSLENGTH", PhotonEnergy, EpotekAbsorption, num);
  // assign this parameter table to the epotek
  epotekMaterial->SetMaterialPropertiesTable(EpotekMPT);

}

void PrtDetectorConstruction::SetVisualization(){

  G4Colour blue = G4Colour(0.0,0.0,1.0);
  G4Colour green = G4Colour(0.0,1.0,.0);
  G4Colour red = G4Colour(1.0,0.0,.0); 
  G4Colour DircColour = G4Colour(1.,1.0,0.);

  G4VisAttributes *waExpHall = new G4VisAttributes(DircColour);
  waExpHall->SetVisibility(false);
  lExpHall->SetVisAttributes(waExpHall);

  G4VisAttributes *waDirc = new G4VisAttributes(DircColour);
  waDirc->SetVisibility(false);
  lDirc->SetVisAttributes(waDirc);

  G4VisAttributes *waBar = new G4VisAttributes(G4Colour(0.,1.,0.9,0.05));
  waBar->SetVisibility(true);
  lBar->SetVisAttributes(waBar);

  G4VisAttributes *waMirror = new G4VisAttributes(G4Colour(1.,1.,0.9,0.2));
  waMirror->SetVisibility(true);
  lMirror->SetVisAttributes(waMirror);

  G4VisAttributes * vaLens = new G4VisAttributes(G4Colour(0.,1.,1.));
  vaLens->SetForceWireframe(true);
  //vaLens->SetForceAuxEdgeVisible(true);
  lLens1->SetVisAttributes(vaLens);
  lLens2->SetVisAttributes(vaLens);
  if(PrtManager::Instance()->GetLens()==3) lLens3->SetVisAttributes(vaLens);

  G4VisAttributes *waPrizm = new G4VisAttributes(G4Colour(0.,0.9,0.9,0.4));
  //waPrizm->SetForceAuxEdgeVisible(true);
  //waPrizm->SetForceSolid(true);
  lPrizm->SetVisAttributes(waPrizm);

  G4VisAttributes *waMcp = new G4VisAttributes(G4Colour(0.1,0.1,0.9,0.3));
  //  waMcp->SetForceWireframe(true);
  waMcp->SetForceSolid(true);
  lMcp->SetVisAttributes(waMcp);

  G4VisAttributes *waPixel = new G4VisAttributes(G4Colour(1.0,0.0,0.1,0.1));
  //waPixel->SetForceWireframe(true);
  waPixel->SetVisibility(false);
  lPixel->SetVisAttributes(waPixel);

}

void PrtDetectorConstruction::ConstructSDandField(){
  // Sensitive detectors
  PrtPixelSD* pixelSD = new PrtPixelSD("PixelSD", "PixelHitsCollection", 0);
  SetSensitiveDetector("lPixel",pixelSD);
  PrtPrizmSD* prizmSD = new PrtPrizmSD("PrizmSD", "PrizmHitsCollection", 0);
  SetSensitiveDetector("lPrizm",prizmSD);
  // Magnetic field
}

void PrtDetectorConstruction::SetRotation(G4double angle){
  fPrtRot->rotateY(-fRotAngle);
  fPrtRot->rotateY(angle);
  fRotAngle=angle;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}


void PrtDetectorConstruction::DrawHitBox(int id){
  Int_t dispopt = PrtManager::Instance()->GetDisplayOpt();
 
  if(dispopt!=1 && dispopt!=2 ) return;

  PrtEvent *event = new PrtEvent();
  PrtHit hit;

  TChain *tree;
  if(dispopt==1) tree = (TChain *) PrtManager::Instance()->GetTree();
  else{
    tree = new TChain("data");
    tree->Add("occuhits.root");
  }

  TBranch *branch = tree->GetBranch("PrtEvent");
  if(!branch) return;
  branch->SetAddress(&event);
  
  Int_t nevent = tree->GetEntries();
  std::cout<<"nevent  "<< nevent <<std::endl;
  
  int pix(0),prizm(0);
  const int maxpixel(50000);
  int load[16][maxpixel]={{0}};
  for (Int_t i=0;i<nevent;i++) {
    tree->GetEvent(i);    
    for(Int_t h=0; h<event->GetHitSize(); h++){
      hit = event->GetHit(h);
      pix=hit.GetPixelId();
      prizm=hit.GetPrizmId();
      if(prizm<0) continue;
      load[prizm][pix]++;
    }
  }

  int maxload[16]={0};
  for(Int_t p=0; p<16; p++){
    for(Int_t i=0; i<maxpixel; i++){
      if(maxload[p]<load[p][i]) maxload[p] =  load[p][i];
    }
  }

  int rgbcolors[256][3] = {{0, 0, 130},{0, 2, 131},{0, 4, 132},{0, 7, 134},{0, 9, 135},{0, 11, 137},{0, 14, 138},{0, 16, 140},{0, 18, 141},{0, 21, 142},{0, 23, 144},{0, 26, 145},{0, 28, 147},{0, 30, 148},{0, 33, 150},{0, 35, 151},{0, 37, 153},{0, 40, 154},{0, 42, 155},{0, 45, 157},{0, 47, 158},{0, 49, 160},{0, 52, 161},{0, 54, 163},{0, 56, 164},{0, 59, 165},{0, 61, 167},{0, 64, 168},{0, 66, 170},{0, 68, 171},{0, 71, 173},{0, 73, 174},{0, 75, 176},{0, 78, 177},{0, 80, 178},{0, 83, 180},{0, 85, 181},{0, 87, 183},{0, 90, 184},{0, 92, 186},{0, 94, 187},{0, 97, 188},{0, 99, 190},{0, 102, 191},{0, 104, 193},{0, 106, 194},{0, 109, 196},{0, 111, 197},{0, 113, 198},{0, 116, 200},{0, 118, 201},{0, 121, 203},{0, 123, 204},{0, 125, 206},{0, 128, 207},{0, 130, 209},{0, 132, 210},{0, 135, 211},{0, 137, 213},{0, 140, 214},{0, 142, 216},{0, 144, 217},{0, 147, 219},{0, 149, 220},{0, 151, 221},{0, 154, 223},{0, 156, 224},{0, 159, 226},{0, 161, 227},{0, 163, 229},{0, 166, 230},{0, 168, 232},{0, 170, 233},{0, 173, 234},{0, 175, 236},{0, 178, 237},{0, 180, 239},{0, 182, 240},{0, 185, 242},{0, 187, 243},{0, 189, 244},{0, 192, 246},{0, 194, 247},{0, 197, 249},{0, 199, 250},{0, 201, 252},{0, 204, 253},{0, 206, 255},{3, 207, 251},{6, 207, 248},{9, 208, 245},{12, 209, 241},{16, 210, 238},{19, 210, 235},{22, 211, 232},{25, 212, 228},{28, 212, 225},{32, 213, 222},{35, 214, 219},{38, 214, 215},{41, 215, 212},{45, 216, 209},{48, 217, 206},{51, 217, 202},{54, 218, 199},{57, 219, 196},{61, 219, 193},{64, 220, 189},{67, 221, 186},{70, 221, 183},{73, 222, 180},{77, 223, 176},{80, 224, 173},{83, 224, 170},{86, 225, 167},{90, 226, 163},{93, 226, 160},{96, 227, 157},{99, 228, 154},{102, 229, 150},{106, 229, 147},{109, 230, 144},{112, 231, 141},{115, 231, 137},{118, 232, 134},{122, 233, 131},{125, 233, 128},{128, 234, 124},{131, 235, 121},{135, 236, 118},{138, 236, 115},{141, 237, 111},{144, 238, 108},{147, 238, 105},{151, 239, 102},{154, 240, 98},{157, 240, 95},{160, 241, 92},{163, 242, 89},{167, 243, 85},{170, 243, 82},{173, 244, 79},{176, 245, 76},{180, 245, 72},{183, 246, 69},{186, 247, 66},{189, 247, 63},{192, 248, 59},{196, 249, 56},{199, 250, 53},{202, 250, 50},{205, 251, 46},{208, 252, 43},{212, 252, 40},{215, 253, 37},{218, 254, 33},{221, 255, 30},{222, 251, 30},{222, 248, 29},{223, 244, 29},{224, 241, 28},{224, 237, 28},{225, 234, 27},{225, 230, 26},{226, 227, 26},{226, 223, 25},{227, 220, 25},{228, 216, 24},{228, 213, 24},{229, 210, 23},{229, 206, 23},{230, 203, 22},{230, 199, 22},{231, 196, 21},{231, 192, 21},{232, 189, 20},{233, 185, 20},{233, 182, 19},{234, 178, 19},{234, 175, 18},{235, 172, 18},{235, 168, 17},{236, 165, 17},{237, 161, 16},{237, 158, 16},{238, 154, 15},{238, 151, 15},{239, 147, 14},{239, 144, 14},{240, 140, 13},{240, 137, 12},{241, 133, 12},{242, 130, 11},{242, 127, 11},{243, 123, 10},{243, 120, 10},{244, 116, 9},{244, 113, 9},{245, 109, 8},{246, 106, 8},{246, 102, 7},{247, 99, 7},{247, 95, 6},{248, 92, 6},{248, 89, 5},{249, 85, 5},{249, 82, 4},{250, 78, 4},{251, 75, 3},{251, 71, 3},{252, 68, 2},{252, 64, 2},{253, 61, 1},{253, 57, 1},{254, 54, 0},{255, 51, 0},{251, 49, 0},{248, 48, 0},{245, 47, 0},{242, 46, 0},{239, 44, 0},{236, 43, 0},{233, 42, 0},{230, 41, 0},{227, 39, 0},{224, 38, 0},{221, 37, 0},{218, 36, 0},{215, 34, 0},{212, 33, 0},{209, 32, 0},{206, 31, 0},{203, 29, 0},{200, 28, 0},{197, 27, 0},{194, 26, 0},{191, 24, 0},{187, 23, 0},{184, 22, 0},{181, 21, 0},{178, 19, 0},{175, 18, 0},{172, 17, 0},{169, 16, 0},{166, 14, 0},{163, 13, 0},{160, 12, 0},{157, 11, 0},{154, 9, 0},{151, 8, 0},{148, 7, 0},{145, 6, 0},{142, 4, 0},{139, 3, 0},{136, 2, 0},{133, 1, 0}};


  G4VisAttributes *waDircHit = new G4VisAttributes(G4Colour(1.,1.,0.9,0.2));
  waDircHit->SetVisibility(false);
	  
  G4double radius = 1000;
  G4double tphi, dphi = 22.5*deg;
  G4int BoxId = 0;
  G4LogicalVolume *lDircHit[16]; 
  for(int i=0; i<16; i++){
    // if(i!=1) continue;
    tphi = dphi*i; 
    G4double dx = (radius+(fPrizm[2]/2.-fPrizm[3]/2.)) * cos(tphi);
    G4double dy = (radius+(fPrizm[2]/2.-fPrizm[3]/2.)) * sin(tphi);
    G4ThreeVector dirc = G4ThreeVector(dx,dy,-500+fBar[2]/2.+fPrizm[1]+fMcpActive[2]/2.+fLens[2]);

    G4Box* gDircHit = new G4Box("gDircHit",400.,300.,10);
    lDircHit[i] = new G4LogicalVolume(gDircHit,defaultMaterial,Form("lDircHit%d",i),0,0,0);
    lDircHit[i]->SetVisAttributes(waDircHit);

    G4RotationMatrix *tRot = new G4RotationMatrix();
    tRot->rotateZ(-tphi);     
    new G4PVPlacement(tRot,dirc,lDircHit[i],"wDircHit",lExpHall,false,i);
    
  }



  G4LogicalVolume * lHit;
  G4VisAttributes *waHit;
  G4Box* gHit;
  for(int p=0; p<16; p++){
    int pixelId=0;
    for(int i=0; i<fNpix1; i++){
      for(int j=0; j<fNpix2; j++){
	double colorid = load[p][pixelId]/(double)maxload[p];
	if(colorid>1) colorid=1;
	double hight = colorid*5;
	int cid = colorid*255;
	if(cid>1){
	
	  gHit = new G4Box("gHit",0.5*fPrizm[2]/fNpix1,0.5*fPrizm[0]/fNpix2,hight+0.2);
	  lHit = new G4LogicalVolume(gHit,BarMaterial,"lHit",0,0,0);
	  waHit = new G4VisAttributes(G4Colour(rgbcolors[cid][0]/255.,rgbcolors[cid][1]/255.,rgbcolors[cid][2]/255.,1.0));
	  waHit->SetForceSolid(true);
	  lHit->SetVisAttributes(waHit);

	  double shiftx = i*(fPrizm[2]/fNpix1) - fPrizm[2]/2.+fPrizm[2]/(2*fNpix1);
	  double shifty = j*(fPrizm[0]/fNpix2) - fPrizm[0]/2.+fPrizm[0]/(2*fNpix2);
	  //	lDirc->GetDaughter(0)
	  //	  pDirc[2]->GetLogicalVolume()->SetVisAttributes(waHit);



	  G4ThreeVector pix = G4ThreeVector(shiftx,shifty,hight);
	  new G4PVPlacement(0,pix,lHit,"wDircHit",lDircHit[p],false,pixelId);
 
	}
	pixelId++;
      }
    }
  }
  
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
}


void PrtDetectorConstruction::SetLens(G4int id){
 
  // if(id==0){
  //   fPrismShift.setZ(fPrismShift.z()-fLens[2]);
  // }
  // std::cout<<"id  "<<id <<std::endl;
  
  // G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void PrtDetectorConstruction::SetQuantumEfficiency(G4int id){
  const G4int num = 36;
  //ideal pmt quantum efficiency
  G4double QuantumEfficiencyIdial[num]=
    {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     1.0,1.0,1.0,1.0,1.0,1.0};

  // Burle PMT's 
  G4double QuantumEfficiencyB[num] =
    {0.,0.001,0.002,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.06,
     0.07,0.09,0.1,0.13,0.15,0.17,0.2,0.24,0.26,0.28,0.282,0.284,0.286,
     0.288,0.29,0.28,0.26,0.24,0.22,0.20,0.18,0.15,0.13,0.12,0.10};
  
  //hamamatsu pmt quantum efficiency
  G4double QuantumEfficiencyPMT[num]=
    {0.001,0.002,0.004,0.007,0.011,0.015,0.020,0.026,0.033,0.040,0.045,
     0.056,0.067,0.085,0.109,0.129,0.138,0.147,0.158,0.170,
     0.181,0.188,0.196,0.203,0.206,0.212,0.218,0.219,0.225,0.230,
     0.228,0.222,0.217,0.210,0.199,0.177};
  
  if(id == 0 ) fQuantumEfficiency = QuantumEfficiencyIdial;
  if(id == 1 ) fQuantumEfficiency = QuantumEfficiencyPMT;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
