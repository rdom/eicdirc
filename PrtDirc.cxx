// simulation software for the Eic DIRC prototype
// author: Roman Dzhygadlo
// contact: r.dzhygadlo at gsi.de

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#include "PrtPhysicsList.h"
#include "PrtDetectorConstruction.h"

#include "PrtActionInitialization.h"
#include "time.h"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4SystemOfUnits.hh"
#include "PrtManager.h"
#include "TVector3.h"
#include "TString.h"

namespace {
  void PrintUsage() {
    G4cerr<<" Usage: "<<G4endl;
    G4cerr<<" Prt [-m macro ] [-u UIsession] [-t nThreads] [-r seed] "<<G4endl;
    G4cerr<<"   note: -t option is available only for multi-threaded mode."<<G4endl;
  }
}

int main(int argc,char** argv)
{
  // Evaluate arguments
  if ( argc > 50 ) {
    PrintUsage();
    return 1;
  }

  G4String macro, outfile, events, geometry, physlist;
  G4String session,geomAng,batchmode,lensId,particle,momentum;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif

  G4long myseed = 3453541;
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro     = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session   = argv[i+1];
    else if ( G4String(argv[i]) == "-r" ) myseed    = atol(argv[i+1]);
    else if ( G4String(argv[i]) == "-o" ) outfile   = argv[i+1];
    else if ( G4String(argv[i]) == "-g" ) geometry  = argv[i+1];
    else if ( G4String(argv[i]) == "-a" ) geomAng   = argv[i+1];
    else if ( G4String(argv[i]) == "-b" ) batchmode = argv[i+1];
    else if ( G4String(argv[i]) == "-e" ) events    = argv[i+1];
    else if ( G4String(argv[i]) == "-l" ) lensId    = argv[i+1];
    else if ( G4String(argv[i]) == "-x" ) particle  = argv[i+1];
    else if ( G4String(argv[i]) == "-p" ) momentum  = argv[i+1];
    else if ( G4String(argv[i]) == "-w" ) physlist  = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  // Seed the random number generator manually
  // if(myseed==345354) myseed = time(NULL);
  
  std::cout<<"SEED "<<myseed <<std::endl;
  G4Random::setTheSeed(myseed);

  if(outfile=="") outfile = "hits.root";
  PrtManager::Instance(outfile);

  if(physlist.size()) PrtManager::Instance()->SetPhysList(atoi(physlist));
  if(geometry.size()) PrtManager::Instance()->SetGeometry(atoi(geometry));
  if(lensId.size())   PrtManager::Instance()->SetLens(atoi(lensId));
 
  // Detector construction
  runManager-> SetUserInitialization(new PrtDetectorConstruction());
  // Physics list
  runManager-> SetUserInitialization(new PrtPhysicsList());
  // User action initialization
  runManager->SetUserInitialization(new PrtActionInitialization(outfile));
  // Initialize G4 kernel
  runManager->Initialize();

#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer(); 
   
  if ( macro.size() ) {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  } else { 
    UImanager->ApplyCommand("/control/execute PrtDirc.mac");
  }
  
  if ( geomAng.size() ) {
    G4String command = "/Prt/geom/prtRotation ";
    UImanager->ApplyCommand(command+geomAng);
    G4String command1 = "/gun/direction ";
    TVector3 rot(1,0,0);
    rot.RotateY((-atoi(geomAng)+90)*deg);
    UImanager->ApplyCommand(command1+Form("%f %f %f",rot.X(),rot.Y(),rot.Z()));

    std::cout<<"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF    "<<atoi(geomAng)*deg<<"  "<<rot.X()<<"  "<<rot.Y()<<"  "<<rot.Z() <<std::endl;
    
  }

  if ( lensId.size() ) {
    G4String command = "/Prt/geom/lensId ";
    UImanager->ApplyCommand(command+lensId);
  }
 
  if ( particle.size() ) {
    G4String command = "/gun/particle ";
    UImanager->ApplyCommand(command+particle);
    int pdgid = 0;
    if(particle=="proton") pdgid = 2212;
    if(particle=="pi+") pdgid = 211;
    if(particle=="pi0") pdgid = 111;
    if(particle=="kaon+") pdgid = 321;
    if(particle=="mu-") pdgid = 13;
    if(particle=="e-") pdgid = 11;

    PrtManager::Instance()->SetParticle(pdgid);
  }

  if ( momentum.size() ) {
    G4String command = "/gun/energy ";
    UImanager->ApplyCommand(command+momentum);
    PrtManager::Instance()->SetMomentum(atof(momentum));
  }

  if ( batchmode.size() ) { // batch mode
    if ( events.size() ) {
      G4String command = "/run/beamOn ";
      UImanager->ApplyCommand(command+events);
    }else{
      UImanager->ApplyCommand("/run/beamOn 1");
    }
  } else {  // UI session for interactive mode
#ifdef G4UI_USE
    G4UIExecutive * ui = new G4UIExecutive(argc,argv,session);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute ../vis.mac");
#endif
    if (ui->IsGUI()) UImanager->ApplyCommand("/control/execute gui.mac");
    if ( events.size() ) {
      G4String command = "/run/beamOn ";
      UImanager->ApplyCommand(command+events);
    }else{
      UImanager->ApplyCommand("/run/beamOn 1");
    }
    ui->SessionStart();
    delete ui;
#endif
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

