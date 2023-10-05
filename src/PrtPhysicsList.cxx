#include "globals.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4FastSimulationManagerProcess.hh"


#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "PrtOpBoundaryProcess.h"
#include "PrtCherenkovProcess.h"
#include "PrtManager.h"

PrtPhysicsList::PrtPhysicsList()
  : G4VUserPhysicsList(), fCerenkovProcess(NULL), fScintillationProcess(NULL),
    fAbsorptionProcess(NULL), fRayleighScatteringProcess(NULL), fMieHGScatteringProcess(NULL),
    fBoundaryProcess(NULL), fMessenger(0) {
  fMessenger = new PrtPhysicsListMessenger(this);
  SetVerboseLevel(0);
  fPhysList = PrtManager::Instance()->getRun()->getPhysList();
  fRunType = PrtManager::Instance()->getRun()->getRunType();
}

PrtPhysicsList::~PrtPhysicsList() {
  delete fMessenger;
}

void PrtPhysicsList::ConstructParticle() {
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

  G4ShortLivedConstructor slConstructor;
  slConstructor.ConstructParticle();
}

void PrtPhysicsList::ConstructProcess() {
  AddTransportation();
  ConstructDecay();
  ConstructEM();
  if(fPhysList == 2) ConstructHad();
  ConstructOp();
  AddParameterisation();
}

#include "G4Decay.hh"
void PrtPhysicsList::ConstructDecay() {
  // Add Decay Process
  G4Decay *theDecayProcess = new G4Decay();
  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()) {
    G4ParticleDefinition *particle = particleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) {
      pmanager->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
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

void PrtPhysicsList::ConstructEM() {
  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()) {
    G4ParticleDefinition *particle = particleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
      // electron
      //  Construct processes for electron
      if (fPhysList != 1 && fPhysList != 11) {
        pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
        pmanager->AddProcess(new G4eBremsstrahlung(), -1, 3, 3);
      }
      pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);
    } else if (particleName == "e+") {
      // positron
      //  Construct processes for positron
      if (fPhysList != 1 && fPhysList != 11) {
        pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
        pmanager->AddProcess(new G4eBremsstrahlung(), -1, 3, 3);
        pmanager->AddProcess(new G4eplusAnnihilation(), 0, -1, 4);
      }
      pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);
    } else if (particleName == "mu+" || particleName == "mu-") {
      // muon
      //  Construct processes for muon
      if (fPhysList != 1 && fPhysList != 11) {
        pmanager->AddProcess(new G4MuMultipleScattering(), -1, 1, 1);
        pmanager->AddProcess(new G4MuBremsstrahlung(), -1, 3, 3);
        pmanager->AddProcess(new G4MuPairProduction(), -1, 4, 4);
      }
      pmanager->AddProcess(new G4MuIonisation(), -1, 2, 2);
    } else {
      if ((particle->GetPDGCharge() != 0.0) && (particle->GetParticleName() != "chargedgeantino") &&
          !particle->IsShortLived()) {
        // all others charged particles except geantino
        if (fPhysList != 1 && fPhysList != 11) {
          pmanager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
        }
        pmanager->AddProcess(new G4hIonisation(), -1, 2, 2);
      }
    }
  }
}

// Hadronic processes ////////////////////////////////////////////////////////
// Elastic processes:
#include "G4HadronElasticProcess.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"

// Inelastic processes:
// #include "G4PionPlusInelasticProcess.hh"
// #include "G4PionMinusInelasticProcess.hh"
// #include "G4KaonPlusInelasticProcess.hh"
// #include "G4KaonZeroSInelasticProcess.hh"
// #include "G4KaonZeroLInelasticProcess.hh"
// #include "G4KaonMinusInelasticProcess.hh"
// #include "G4ProtonInelasticProcess.hh"
// #include "G4AntiProtonInelasticProcess.hh"
// #include "G4NeutronInelasticProcess.hh"
// #include "G4AntiNeutronInelasticProcess.hh"
// #include "G4DeuteronInelasticProcess.hh"
// #include "G4TritonInelasticProcess.hh"
// #include "G4AlphaInelasticProcess.hh"

// // High energy FTFP model and Bertini cascade
// #include "G4FTFModel.hh"
// #include "G4LundStringFragmentation.hh"
// #include "G4ExcitedStringDecay.hh"
// #include "G4PreCompoundModel.hh"
// #include "G4GeneratorPrecompoundInterface.hh"
// #include "G4TheoFSGenerator.hh"
// #include "G4CascadeInterface.hh"

// #include "G4HadronicParameters.hh"

// // Cross sections
// #include "G4VCrossSectionDataSet.hh"
// #include "G4CrossSectionDataSetRegistry.hh"

// #include "G4CrossSectionElastic.hh"
// #include "G4CrossSectionInelastic.hh"
// #include "G4BGGPionElasticXS.hh"
// #include "G4BGGPionInelasticXS.hh"
// #include "G4AntiNuclElastic.hh"

// #include "G4CrossSectionInelastic.hh"
// #include "G4CrossSectionPairGG.hh"
// #include "G4BGGNucleonInelasticXS.hh"
// #include "G4BGGNucleonElasticXS.hh"
// #include "G4NeutronInelasticXS.hh"
// #include "G4NeutronElasticXS.hh"
// #include "G4ComponentAntiNuclNuclearXS.hh"
// #include "G4ComponentGGNuclNuclXsc.hh"
// #include "G4ComponentGGHadronNucleusXsc.hh"

// #include "G4HadronElastic.hh"
// #include "G4HadronCaptureProcess.hh"

// // Stopping processes
// #include "G4PiMinusAbsorptionBertini.hh"
// #include "G4KaonMinusAbsorptionBertini.hh"
// #include "G4AntiProtonAbsorptionFritiof.hh"
// #include "G4HadronicParameters.hh"

// #include "G4Decay.hh"
// #include "G4RadioactiveDecayBase.hh"
// #include "G4PhysicsListHelper.hh"
// #include "G4NuclideTable.hh"
// #include "G4NuclearLevelData.hh"

// // Neutron high-precision models: <20 MeV
// #include "G4ParticleHPElastic.hh"
// #include "G4ParticleHPElasticData.hh"
// #include "G4ParticleHPCapture.hh"
// #include "G4ParticleHPCaptureData.hh"
// #include "G4ParticleHPInelastic.hh"
// #include "G4ParticleHPInelasticData.hh"


void PrtPhysicsList::ConstructHad() 
{
  // //Elastic models
  // G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  // G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  // G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE(); 
  
  // // Inelastic scattering
  // const G4double theFTFMin0 =    0.0*GeV;
  // const G4double theFTFMin1 =    3.0*GeV;
  // const G4double theFTFMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  // const G4double theBERTMin0 =   0.0*GeV;
  // const G4double theBERTMin1 =  19.0*MeV;
  // const G4double theBERTMax =    6.0*GeV;
  
  // G4FTFModel * theStringModel = new G4FTFModel;
  // G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
  // theStringModel->SetFragmentationModel( theStringDecay );
  // G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  // G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  // G4TheoFSGenerator * theFTFModel0 = new G4TheoFSGenerator( "FTFP" );
  // theFTFModel0->SetHighEnergyGenerator( theStringModel );
  // theFTFModel0->SetTransport( theCascade );
  // theFTFModel0->SetMinEnergy( theFTFMin0 );
  // theFTFModel0->SetMaxEnergy( theFTFMax );

  // G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
  // theFTFModel1->SetHighEnergyGenerator( theStringModel );
  // theFTFModel1->SetTransport( theCascade );
  // theFTFModel1->SetMinEnergy( theFTFMin1 );
  // theFTFModel1->SetMaxEnergy( theFTFMax );

  // G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
  // theBERTModel0->SetMinEnergy( theBERTMin0 );
  // theBERTModel0->SetMaxEnergy( theBERTMax );

  // G4CascadeInterface * theBERTModel1 = new G4CascadeInterface;
  // theBERTModel1->SetMinEnergy( theBERTMin1 );
  // theBERTModel1->SetMaxEnergy( theBERTMax );

  // G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  // G4ComponentGGHadronNucleusXsc * ggHNXsec = new G4ComponentGGHadronNucleusXsc();
  // G4VCrossSectionDataSet * theGGHNEl = new G4CrossSectionElastic(ggHNXsec);
  // G4VCrossSectionDataSet * theGGHNInel = new G4CrossSectionInelastic(ggHNXsec);

  // auto particleIterator=GetParticleIterator();
  // particleIterator->reset();
  // while ((*particleIterator)())
  //   {
  //     G4ParticleDefinition* particle = particleIterator->value();
  //     G4ProcessManager* pmanager = particle->GetProcessManager();
  //     G4String particleName = particle->GetParticleName();

  //     if (particleName == "pi+") 
  // 	{
  // 	  // Elastic scattering
  //         G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  //         theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
  //         theElasticProcess->RegisterMe( elastic_he );
  // 	  pmanager->AddDiscreteProcess( theElasticProcess );
  // 	  // //Inelastic scattering
  // 	  // G4PionPlusInelasticProcess* theInelasticProcess = 
  // 	  //   new G4PionPlusInelasticProcess("inelastic");
  // 	  // theInelasticProcess->AddDataSet( new G4BGGPionInelasticXS( particle ) );
  // 	  // theInelasticProcess->RegisterMe( theFTFModel1 );
  //         // theInelasticProcess->RegisterMe( theBERTModel0 );
  // 	  // pmanager->AddDiscreteProcess( theInelasticProcess );
  // 	} 

  //     else if (particleName == "pi-") 
  // 	{
  // 	  // Elastic scattering
  //         G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  //         theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
  //         theElasticProcess->RegisterMe( elastic_he );
  // 	  pmanager->AddDiscreteProcess( theElasticProcess );
  // 	  // //Inelastic scattering
  // 	  // G4PionMinusInelasticProcess* theInelasticProcess = 
  // 	  //   new G4PionMinusInelasticProcess("inelastic");
  // 	  // theInelasticProcess->AddDataSet( new G4BGGPionInelasticXS( particle ) );
  // 	  // theInelasticProcess->RegisterMe( theFTFModel1 );
  //         // theInelasticProcess->RegisterMe( theBERTModel0 );
  // 	  // pmanager->AddDiscreteProcess( theInelasticProcess );	  
  // 	  //Absorption
  // 	  pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);
  // 	}
  //     else if (particleName == "kaon+") 
  // 	{
  // 	  // Elastic scattering
  //         G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  // 	  theElasticProcess->AddDataSet( theGGHNEl );
  //         theElasticProcess->RegisterMe( elastic_lhep0 );
  // 	  pmanager->AddDiscreteProcess( theElasticProcess );
  //         // // Inelastic scattering	
  // 	  // G4KaonPlusInelasticProcess* theInelasticProcess = 
  // 	  //   new G4KaonPlusInelasticProcess("inelastic");
  // 	  // theInelasticProcess->AddDataSet( theGGHNInel );
  // 	  // theInelasticProcess->RegisterMe( theFTFModel1 );
  //         // theInelasticProcess->RegisterMe( theBERTModel0 );
  // 	  // pmanager->AddDiscreteProcess( theInelasticProcess );
  // 	}      
  //     else if (particleName == "kaon0S") 
  // 	{
  // 	  // Elastic scattering
  //         G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  // 	  theElasticProcess->AddDataSet( theGGHNEl );
  //         theElasticProcess->RegisterMe( elastic_lhep0 );
  // 	  pmanager->AddDiscreteProcess( theElasticProcess );
  //         // Inelastic scattering	 
  // 	  G4KaonZeroSInelasticProcess* theInelasticProcess = 
  // 	    new G4KaonZeroSInelasticProcess("inelastic");
  // 	  theInelasticProcess->AddDataSet( theGGHNInel );
  // 	  theInelasticProcess->RegisterMe( theFTFModel1 );
  //         theInelasticProcess->RegisterMe( theBERTModel0 );
  // 	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
  // 	}

  //     else if (particleName == "kaon0L") 
  // 	{
  // 	  // Elastic scattering
  //         G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  // 	  theElasticProcess->AddDataSet( theGGHNEl );
  //         theElasticProcess->RegisterMe( elastic_lhep0 );
  // 	  pmanager->AddDiscreteProcess( theElasticProcess );
  // 	  // Inelastic scattering
  // 	  G4KaonZeroLInelasticProcess* theInelasticProcess = 
  // 	    new G4KaonZeroLInelasticProcess("inelastic");
  // 	  theInelasticProcess->AddDataSet( theGGHNInel );
  // 	  theInelasticProcess->RegisterMe( theFTFModel1 );
  //         theInelasticProcess->RegisterMe( theBERTModel0 ); 
  // 	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
  // 	}

  //     else if (particleName == "kaon-") 
  // 	{
  // 	  // Elastic scattering
  //         G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  // 	  theElasticProcess->AddDataSet( theGGHNEl );
  //         theElasticProcess->RegisterMe( elastic_lhep0 );
  // 	  pmanager->AddDiscreteProcess( theElasticProcess );
  //         // Inelastic scattering
  // 	  G4KaonMinusInelasticProcess* theInelasticProcess = 
  // 	    new G4KaonMinusInelasticProcess("inelastic");	
  //         theInelasticProcess->AddDataSet( theGGHNInel );
  // 	  theInelasticProcess->RegisterMe( theFTFModel1 );
  //         theInelasticProcess->RegisterMe( theBERTModel0 );
  // 	  pmanager->AddDiscreteProcess( theInelasticProcess );
  // 	  pmanager->AddRestProcess(new G4KaonMinusAbsorptionBertini, ordDefault);
  // 	}

  //     else if (particleName == "proton") 
  // 	{
  // 	  // Elastic scattering
  //         G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  //         theElasticProcess->AddDataSet( new G4BGGNucleonElasticXS( G4Proton::Proton() ) );
  //         theElasticProcess->RegisterMe( elastic_chip );
  // 	  pmanager->AddDiscreteProcess( theElasticProcess );
  // 	  // Inelastic scattering
  // 	  G4ProtonInelasticProcess* theInelasticProcess = 
  // 	    new G4ProtonInelasticProcess("inelastic");
  // 	  theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
  // 	  theInelasticProcess->RegisterMe( theFTFModel1 );
  //         theInelasticProcess->RegisterMe( theBERTModel0 );
  // 	  pmanager->AddDiscreteProcess( theInelasticProcess );
  // 	}
  //     else if (particleName == "anti_proton") 
  // 	{
  // 	  // Elastic scattering
  //         const G4double elastic_elimitAntiNuc = 100.0*MeV;
  //         G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
  //         elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
  //         G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic( elastic_anuc->GetComponentCrossSection() );
  //         G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
  //         elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );
  //         G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  //         theElasticProcess->AddDataSet( elastic_anucxs );
  //         theElasticProcess->RegisterMe( elastic_lhep2 );
  //         theElasticProcess->RegisterMe( elastic_anuc );
  // 	  pmanager->AddDiscreteProcess( theElasticProcess );
  // 	  // Inelastic scattering
  // 	  G4AntiProtonInelasticProcess* theInelasticProcess = 
  // 	    new G4AntiProtonInelasticProcess("inelastic");
  // 	  theInelasticProcess->AddDataSet( theAntiNucleonData );
  // 	  theInelasticProcess->RegisterMe( theFTFModel0 );
  // 	  pmanager->AddDiscreteProcess( theInelasticProcess );
  // 	  // Absorption
  // 	  pmanager->AddRestProcess(new G4AntiProtonAbsorptionFritiof, ordDefault);
  // 	}
  //   }
}

void PrtPhysicsList::ConstructOp() {
  fCerenkovProcess = new G4Cerenkov("Cerenkov");
  fCerenkovProcess0 = new PrtCherenkovProcess("Cerenkov");
  fScintillationProcess = new G4Scintillation("Scintillation");
  fAbsorptionProcess = new G4OpAbsorption();
  fRayleighScatteringProcess = new G4OpRayleigh();
  fMieHGScatteringProcess = new G4OpMieHG();
  fBoundaryProcess = new PrtOpBoundaryProcess();

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

  // fScintillationPrSetScintillationYieldFactor(1.);
  fScintillationProcess->SetTrackSecondariesFirst(true);

  // Use Birks Correction in the Scintillation process

  G4EmSaturation *emSaturation = G4LossTableManager::Instance()->EmSaturation();
  fScintillationProcess->AddSaturation(emSaturation);

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()) {
    G4ParticleDefinition *particle = particleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (fCerenkovProcess->IsApplicable(*particle)) {
      if (fPhysList < 10) {
        pmanager->AddProcess(fCerenkovProcess);
        pmanager->SetProcessOrdering(fCerenkovProcess, idxPostStep);
      } else {
        pmanager->AddProcess(fCerenkovProcess0);
        pmanager->SetProcessOrdering(fCerenkovProcess0, idxPostStep);
      }
    }

    if (fScintillationProcess->IsApplicable(*particle)) {
      // pmanager->AddProcess(fScintillationProcess);
      // pmanager->SetProcessOrderingToLast(fScintillationProcess, idxAtRest);
      // pmanager->SetProcessOrderingToLast(fScintillationProcess, idxPostStep);
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

void PrtPhysicsList::SetVerbose(G4int verbose) {
  fCerenkovProcess->SetVerboseLevel(verbose);
  fCerenkovProcess0->SetVerboseLevel(verbose);
  fScintillationProcess->SetVerboseLevel(verbose);
  fAbsorptionProcess->SetVerboseLevel(verbose);
  fRayleighScatteringProcess->SetVerboseLevel(verbose);
  fMieHGScatteringProcess->SetVerboseLevel(verbose);
  fBoundaryProcess->SetVerboseLevel(verbose);
}

void PrtPhysicsList::SetNbOfPhotonsCerenkov(G4int MaxNumber) {
  fCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumber);
  fCerenkovProcess0->SetMaxNumPhotonsPerStep(MaxNumber);
}

void PrtPhysicsList::SetCuts() {
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  //
  SetCutsWithDefault();

  if (verboseLevel > 0) DumpCutValuesTable();
}

void PrtPhysicsList::AddParameterisation() {

  auto fastsim = new G4FastSimulationManagerProcess();
  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()) {
    G4ParticleDefinition *particle = particleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    pmanager->AddProcess(fastsim, -1, 0, 0);
  }
}
