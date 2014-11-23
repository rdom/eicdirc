#ifndef PrtPhysicsListMessenger_h
#define PrtPhysicsListMessenger_h

#include "globals.hh"
#include "G4UImessenger.hh"

class PrtPhysicsList;
class G4UIdirectory;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrtPhysicsListMessenger: public G4UImessenger
{
  public:
    PrtPhysicsListMessenger(PrtPhysicsList* );
    virtual ~PrtPhysicsListMessenger();
 
    virtual void SetNewValue(G4UIcommand*, G4String);
 
  private:
    PrtPhysicsList*  fPhysicsList;
 
    G4UIdirectory*        fPrtDir;
    G4UIdirectory*        fPhysDir;
    G4UIcmdWithAnInteger* fVerboseCmd;
    G4UIcmdWithAnInteger* fCerenkovCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
