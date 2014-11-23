#include "PrtDetectorConstructionMessenger.h"
#include "PrtPrimaryGeneratorAction.h"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"

#include "PrtManager.h"

PrtDetectorConstructionMessenger::PrtDetectorConstructionMessenger
(PrtDetectorConstruction *PrtGeom): G4UImessenger(),
    fPrtGeom(PrtGeom)
{
  fGeomDir = new G4UIdirectory("/Prt/geom/");
  fGeomDir->SetGuidance("Geometry control");

  fAngleCmd = new G4UIcmdWithADoubleAndUnit("/Prt/geom/prtRotation",this);
  fAngleCmd->SetGuidance("Rotation angle of the prototype");
  fAngleCmd->SetParameterName("angle",true);
  fAngleCmd->SetRange("angle>=0. && angle<180.");
  fAngleCmd->SetDefaultValue(0.);
  fAngleCmd->SetDefaultUnit("deg");

  fLensIdCmd = new G4UIcmdWithAnInteger("/Prt/geom/lensId",this);
  fLensIdCmd->SetGuidance("Lens Id");
  fLensIdCmd->SetParameterName("lenseId",true);
  fLensIdCmd->SetRange("lenseId>=0");
  fLensIdCmd->SetDefaultValue(1);

  fDetEffType = new G4UIcmdWithAnInteger("/Prt/geom/detEffType",this);
  fDetEffType->SetGuidance("Type of the detector efficiency");
  fDetEffType->SetParameterName("detEffType",true);
  fDetEffType->SetRange("detEffType>=0");
  fDetEffType->SetDefaultValue(0);
}

PrtDetectorConstructionMessenger::~PrtDetectorConstructionMessenger()
{
  delete fAngleCmd;
  delete fGeomDir;
}

void PrtDetectorConstructionMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{ 
  if( command == fAngleCmd ) {
    G4double angle = 180*deg - fAngleCmd->GetNewDoubleValue(newValue);// + 90*deg;
    fPrtGeom->SetRotation(angle);
    PrtManager::Instance()->SetAngle(angle);
    std::cout<<"angle deg   "<<angle/deg <<std::endl;
    
  } 
  
  if( command == fLensIdCmd ) {
    G4int id = fLensIdCmd->GetNewIntValue(newValue);
    PrtManager::Instance()->SetLens(id);
    fPrtGeom->SetLens(id);
  }

  if( command == fDetEffType ) {
    G4int id = fDetEffType->GetNewIntValue(newValue);
    //PrtManager::Instance()->SetQuantumEfficiency(id);
    //PrtManager::Instance()->SetDetEffType();
  }
}
