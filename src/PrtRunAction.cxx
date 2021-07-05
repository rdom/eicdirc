// Make this appear first!
#include "G4Timer.hh"
#include "PrtRunAction.h"
#include "G4Run.hh"
#include "PrtManager.h"


PrtRunAction::PrtRunAction()
 : G4UserRunAction(),
   fTimer(0)
{
  fTimer = new G4Timer;
}

PrtRunAction::~PrtRunAction()
{
  delete fTimer;
}

void PrtRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fTimer->Start();  
}

void PrtRunAction::EndOfRunAction(const G4Run* aRun)
{
  if(IsMaster())
    {
      PrtManager::Instance()->fillLut();
      PrtManager::Instance()->save();
  
      fTimer->Stop();
      G4cout << "number of event = " << aRun->GetNumberOfEvent()
	     << " " << *fTimer << G4endl;
    }
}

