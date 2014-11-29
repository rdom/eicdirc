#ifndef PrtDetectorConstructionMessenger_h
#define PrtDetectorConstructionMessenger_h

#include "G4UImessenger.hh"
#include "globals.hh"

#include "PrtDetectorConstruction.h"

class PrtDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;


class PrtDetectorConstructionMessenger: public G4UImessenger
{ 
public:
  PrtDetectorConstructionMessenger(PrtDetectorConstruction*);
  virtual ~PrtDetectorConstructionMessenger();
 
  virtual void SetNewValue(G4UIcommand*, G4String);
 
private:
  PrtDetectorConstruction*        fPrtGeom;
  G4UIdirectory*                  fGeomDir;
  G4UIcmdWithADoubleAndUnit*      fAngleCmd;
  G4UIcmdWithAnInteger*           fLensIdCmd;  
  G4UIcmdWithAnInteger*           fDetEffType;
  G4UIcmdWithAnInteger*           fDrawHits;
};

#endif
