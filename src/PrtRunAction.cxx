// Make this appear first!
#include "G4Timer.hh"
#include "PrtRunAction.h"
#include "G4Run.hh"
#include "PrtManager.h"
#include "G4UImanager.hh"

PrtRunAction::PrtRunAction(G4String outfile)
 : G4UserRunAction(),
   fTimer(0)
{
  fTimer = new G4Timer;
  fOutFile = outfile;
}

PrtRunAction::~PrtRunAction()
{
  delete fTimer;
}

void PrtRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fTimer->Start();


  PrtManager::Instance(fOutFile);
  
}

void PrtRunAction::EndOfRunAction(const G4Run* aRun)
{
  
  G4UImanager::GetUIpointer()->ApplyCommand("/Prt/geom/drawHits");
  PrtManager::Instance()->Save();
  
  fTimer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent()
         << " " << *fTimer << G4endl;
}

