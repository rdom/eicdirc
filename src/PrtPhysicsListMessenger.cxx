#include "PrtPhysicsListMessenger.h"

#include "PrtPhysicsList.h"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

PrtPhysicsListMessenger::
  PrtPhysicsListMessenger(PrtPhysicsList* pPhys) 
  : G4UImessenger(),
    fPhysicsList(pPhys)
{
  fPrtDir = new G4UIdirectory("/Prt/");
  fPrtDir->SetGuidance("UI commands of this example");

  fPhysDir = new G4UIdirectory("/Prt/phys/");
  fPhysDir->SetGuidance("PhysicsList control");
 
  fVerboseCmd = new G4UIcmdWithAnInteger("/Prt/phys/verbose",this);
  fVerboseCmd->SetGuidance("set verbose for physics processes");
  fVerboseCmd->SetParameterName("verbose",true);
  fVerboseCmd->SetDefaultValue(1);
  fVerboseCmd->SetRange("verbose>=0");
  fVerboseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  fCerenkovCmd =
           new G4UIcmdWithAnInteger("/Prt/phys/cerenkovMaxPhotons",this);
  fCerenkovCmd->SetGuidance("set max nb of photons per step");
  fCerenkovCmd->SetParameterName("MaxNumber",false);
  fCerenkovCmd->SetRange("MaxNumber>=0");
  fCerenkovCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrtPhysicsListMessenger::~PrtPhysicsListMessenger()
{
  delete fVerboseCmd;
  delete fCerenkovCmd;
  delete fPhysDir;
  delete fPrtDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrtPhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{
  if( command == fVerboseCmd )
   {fPhysicsList->SetVerbose(fVerboseCmd->GetNewIntValue(newValue));}

  if( command == fCerenkovCmd )
   {fPhysicsList->
              SetNbOfPhotonsCerenkov(fCerenkovCmd->GetNewIntValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
