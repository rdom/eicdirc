#ifndef PrtPrimaryGeneratorMessenger_h
#define PrtPrimaryGeneratorMessenger_h

#include "G4UImessenger.hh"
#include "globals.hh"

class PrtPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrtPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    PrtPrimaryGeneratorMessenger(PrtPrimaryGeneratorAction* );
    virtual ~PrtPrimaryGeneratorMessenger();
 
    virtual void SetNewValue(G4UIcommand*, G4String);
 
  private:
    PrtPrimaryGeneratorAction* fPrtAction;
    G4UIdirectory*                  fGunDir;
    G4UIcmdWithADoubleAndUnit*      fPolarCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
