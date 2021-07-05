#ifndef PrtActionInitialization_h
#define PrtActionInitialization_h 1

#include "G4String.hh"
#include "G4VUserActionInitialization.hh"

class PrtActionInitialization : public G4VUserActionInitialization
{
public:
  PrtActionInitialization();
  virtual ~PrtActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;

  virtual G4VSteppingVerbose* InitializeSteppingVerbose() const;

private:
  G4String fOutFile;
};

#endif
