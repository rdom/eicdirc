#include "globals.hh"
#include "PrtPhysicsList.h"
#include "PrtPhysicsListMessenger.h"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ProcessManager.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "PrtOpBoundaryProcess.h"
#include "PrtCherenkovProcess.h"
#include "PrtManager.h"

PrtPhysicsList::PrtPhysicsList() 
 : G4VUserPhysicsList(),
   fCerenkovProcess(NULL),
   fScintillationProcess(NULL),
   fAbsorptionProcess(NULL),
   fRayleighScatteringProcess(NULL),
   fMieHGScatteringProcess(NULL),
   fBoundaryProcess(NULL),
   fMessenger(0) 
{
  fMessenger = new PrtPhysicsListMessenger(this);
  SetVerboseLevel(0);
  fPhysList =  PrtManager::Instance()->GetPhysList();
}

PrtPhysicsList::~PrtPhysicsList() { delete fMessenger; }


void PrtPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  G4BosonConstructor bConstructor;
  bConstructor.ConstructParticle();

  G4LeptonConstructor lConstructor;
  lConstructor.ConstructParticle();

  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

  G4BaryonConstructor rConstructor;
  rConstructor.ConstructParticle();

  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle(); 
}

void PrtPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructDecay();
  ConstructEM();
  ConstructOp();
}

#include "G4Decay.hh"

void PrtPhysicsList::ConstructDecay()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) {
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

void PrtPhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      if(fPhysList !=1 && fPhysList != 11) pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation(),       -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung(),   -1, 3, 3);

    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
      if(fPhysList !=1 && fPhysList != 11) pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation(),       -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung(),   -1, 3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation(),  0,-1, 4);

    } else if( particleName == "mu+" || particleName == "mu-"    ) {
      //muon
      // Construct processes for muon
     if(fPhysList !=1 && fPhysList != 11) pmanager->AddProcess(new G4MuMultipleScattering(),-1, 1, 1);
     pmanager->AddProcess(new G4MuIonisation(),      -1, 2, 2);
     pmanager->AddProcess(new G4MuBremsstrahlung(),  -1, 3, 3);
     pmanager->AddProcess(new G4MuPairProduction(),  -1, 4, 4);

    } else {
      if ((particle->GetPDGCharge() != 0.0) &&
          (particle->GetParticleName() != "chargedgeantino") &&
          !particle->IsShortLived()) {
       // all others charged particles except geantino
	if(fPhysList !=1 && fPhysList != 11) pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
	pmanager->AddProcess(new G4hIonisation(),       -1,2,2);
     }
    }
  }
}

void PrtPhysicsList::ConstructOp() 
{
  fCerenkovProcess           = new G4Cerenkov("Cerenkov");
  fCerenkovProcess0          = new PrtCherenkovProcess("Cerenkov");
  fScintillationProcess      = new G4Scintillation("Scintillation");
  fAbsorptionProcess         = new G4OpAbsorption();
  fRayleighScatteringProcess = new G4OpRayleigh();
  fMieHGScatteringProcess    = new G4OpMieHG();
  fBoundaryProcess           = new PrtOpBoundaryProcess();

//  fCerenkovProcess->DumpPhysicsTable();
//  fScintillationProcess->DumpPhysicsTable();
//  fRayleighScatteringProcess->DumpPhysicsTable();

  SetVerbose(0);
  
  fCerenkovProcess->SetMaxNumPhotonsPerStep(20);
  fCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  fCerenkovProcess->SetTrackSecondariesFirst(true);

  fCerenkovProcess0->SetMaxNumPhotonsPerStep(20);
  fCerenkovProcess0->SetMaxBetaChangePerStep(10.0);
  fCerenkovProcess0->SetTrackSecondariesFirst(true);
  
  fScintillationProcess->SetScintillationYieldFactor(1.);
  fScintillationProcess->SetTrackSecondariesFirst(true);

  // Use Birks Correction in the Scintillation process

  G4EmSaturation* emSaturation =
                               G4LossTableManager::Instance()->EmSaturation();
  fScintillationProcess->AddSaturation(emSaturation);

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (fCerenkovProcess->IsApplicable(*particle)) {
      if(fPhysList < 10){
	pmanager->AddProcess(fCerenkovProcess);
	pmanager->SetProcessOrdering(fCerenkovProcess,idxPostStep);
      }else{
	pmanager->AddProcess(fCerenkovProcess0);
	pmanager->SetProcessOrdering(fCerenkovProcess0,idxPostStep);
      }
    }
  
    if (fScintillationProcess->IsApplicable(*particle)) {
      //pmanager->AddProcess(fScintillationProcess);
      //pmanager->SetProcessOrderingToLast(fScintillationProcess, idxAtRest);
      //pmanager->SetProcessOrderingToLast(fScintillationProcess, idxPostStep);
    }
    if (particleName == "opticalphoton") {
      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
      pmanager->AddDiscreteProcess(fAbsorptionProcess);
      pmanager->AddDiscreteProcess(fRayleighScatteringProcess);
      pmanager->AddDiscreteProcess(fMieHGScatteringProcess);
      pmanager->AddDiscreteProcess(fBoundaryProcess);
    }
  }
}

void PrtPhysicsList::SetVerbose(G4int verbose)
{
  fCerenkovProcess->SetVerboseLevel(verbose);
  fCerenkovProcess0->SetVerboseLevel(verbose);
  fScintillationProcess->SetVerboseLevel(verbose);
  fAbsorptionProcess->SetVerboseLevel(verbose);
  fRayleighScatteringProcess->SetVerboseLevel(verbose);
  fMieHGScatteringProcess->SetVerboseLevel(verbose);
  fBoundaryProcess->SetVerboseLevel(verbose);
}

void PrtPhysicsList::SetNbOfPhotonsCerenkov(G4int MaxNumber)
{
  fCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumber);
  fCerenkovProcess0->SetMaxNumPhotonsPerStep(MaxNumber);
}

void PrtPhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  //
  SetCutsWithDefault();

  if (verboseLevel>0) DumpCutValuesTable();
}

