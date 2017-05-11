// simulation software for the Eic DIRC prototype
// author: Roman Dzhygadlo
// contact: r.dzhygadlo at gsi.de

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "TROOT.h"

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

#include "TApplication.h"

#include "PrtManager.h"
#include "PrtLutReco.h"



namespace {
  void PrintUsage() {
    G4cerr<<" read README.md "<<G4endl;
  }
}

int main(int argc,char** argv)
{
  for ( G4int i=1; i<argc; i=i+2 ) std::cout<< argv[i] << "  "<<argv[i+1] <<std::endl;

  // Evaluate arguments
  if ( argc > 50 ) {
    PrintUsage();
    return 1;
  }

  TApplication theApp("App", 0, 0);
  
  G4String macro, events, geometry,evType, radiator, physlist, outfile, 
    session,geomAng,batchmode,lensId,particle("pi+"),momentum("3.5 GeV"),timeRes,displayOpt,
    beamDimension, mcpLayout, infile = "hits.root", lutfile = "../data/lut_avr.root";
  G4int firstevent(0), runtype(0), verbose(0);

  G4long myseed = 345354;
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro     = argv[i+1];
    //    else if ( G4String(argv[i]) == "-u" ) session   = argv[i+1];
    else if ( G4String(argv[i]) == "-r" ) myseed    = atol(argv[i+1]);
    else if ( G4String(argv[i]) == "-o" ) outfile   = argv[i+1];
    else if ( G4String(argv[i]) == "-i" ) infile    = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) lutfile   = argv[i+1];
    else if ( G4String(argv[i]) == "-g" ) geometry  = argv[i+1];
    else if ( G4String(argv[i]) == "-h" ) radiator  = argv[i+1];
    else if ( G4String(argv[i]) == "-a" ) geomAng   = argv[i+1];
    else if ( G4String(argv[i]) == "-b" ) batchmode = argv[i+1];
    else if ( G4String(argv[i]) == "-f" ) firstevent= atoi(argv[i+1]);
    else if ( G4String(argv[i]) == "-e" ) events    = argv[i+1];
    else if ( G4String(argv[i]) == "-l" ) lensId    = argv[i+1];
    else if ( G4String(argv[i]) == "-x" ) particle  = argv[i+1];
    else if ( G4String(argv[i]) == "-p" ) momentum  = argv[i+1];
    else if ( G4String(argv[i]) == "-w" ) physlist  = argv[i+1];
    else if ( G4String(argv[i]) == "-s" ) runtype   = atoi(argv[i+1]);
    else if ( G4String(argv[i]) == "-z" ) beamDimension  = argv[i+1];
    else if ( G4String(argv[i]) == "-c" ) mcpLayout = argv[i+1];
    else if ( G4String(argv[i]) == "-d" ) displayOpt= argv[i+1];
    else if ( G4String(argv[i]) == "-tr" ) timeRes   = argv[i+1];
    else if ( G4String(argv[i]) == "-v" ) verbose   = atoi(argv[i+1]);
    else {
      PrintUsage();
      return 1;
    }
  }

  if(outfile=="" && (runtype == 0 || runtype == 10)) outfile = "hits.root"; // simulation
  if(outfile=="" && (runtype == 1 || runtype == 5)) outfile = "../data/lut.root"; // lookup table generation
  if(outfile=="" && (runtype == 2 || runtype==3 || runtype==4)) outfile = "reco.root"; // reconstruction

  if(batchmode.size()) gROOT->SetBatch(kTRUE);
  if(!events.size()) events = "0";
  PrtManager::Instance(outfile,runtype);

  if(physlist.size()) PrtManager::Instance()->SetPhysList(atoi(physlist));
  if(geometry.size()) PrtManager::Instance()->SetGeometry(atoi(geometry));
  if(evType.size())   PrtManager::Instance()->SetEvType(atoi(evType));
  if(radiator.size()) PrtManager::Instance()->SetRadiator(atoi(radiator));
  if(lensId.size())   PrtManager::Instance()->SetLens(atoi(lensId));
  if(mcpLayout.size())PrtManager::Instance()->SetMcpLayout(atoi(mcpLayout));
  if(beamDimension.size()) PrtManager::Instance()->SetBeamDimension(atoi(beamDimension));
  if(displayOpt.size())   PrtManager::Instance()->SetDisplayOpt(atoi(displayOpt));
  if(timeRes.size())   PrtManager::Instance()->SetTimeRes(atof(timeRes));
  if(geomAng.size())   PrtManager::Instance()->SetAngle(atof(geomAng));

  if(runtype == 2 || runtype==3 || runtype==4){
    PrtLutReco * reco = new PrtLutReco(infile.c_str(),lutfile.c_str(),verbose);
    reco->Run(firstevent, atoi(events));
    return 0;
  }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  std::cout<<"SEED "<<myseed <<std::endl;
  G4Random::setTheSeed(myseed);
  G4RunManager * runManager = new G4RunManager;


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
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer(); 
   
  if ( macro.size() ) {
    UImanager->ApplyCommand("/control/execute "+macro);
  } else { 
    //UImanager->ApplyCommand("/control/execute ../prt.mac");
  }
  
  if ( geomAng.size() ) {
    UImanager->ApplyCommand("/Prt/geom/prtRotation "+geomAng);
  }

  if ( lensId.size() ) {
    UImanager->ApplyCommand("/Prt/geom/lensId "+lensId);
  }
 
  int pdgid = 2212;
  if(particle=="proton") pdgid = 2212;
  if(particle=="pi+") pdgid = 211;
  if(particle=="pi0") pdgid = 111;
  if(particle=="kaon+") pdgid = 321;
  if(particle=="mu-") pdgid = 13;
  if(particle=="e-") pdgid = 11;
  if(particle=="mix"){
    PrtManager::Instance()->SetMixPiK(true);
    particle="pi+";
  }
  
  PrtManager::Instance()->SetParticle(pdgid);
  UImanager->ApplyCommand("/gun/particle "+particle);
  UImanager->ApplyCommand("/gun/momentumAmp "+momentum);

  if(batchmode.size()){ // batch mode
    if(runtype==10){
      PrtManager::Instance()->SetParticle(211);
      UImanager->ApplyCommand("/gun/particle pi+");
      UImanager->ApplyCommand("/gun/momentumAmp "+momentum);
      UImanager->ApplyCommand("/run/beamOn "+events);

      PrtManager::Instance()->SetParticle(321);
      UImanager->ApplyCommand("/gun/particle kaon+");
      UImanager->ApplyCommand("/gun/momentumAmp "+momentum);
    }
    
    UImanager->ApplyCommand("/run/beamOn "+events);
  }else{  // UI session for interactive mode

#ifdef G4UI_USE
    G4UIExecutive * ui = new G4UIExecutive(argc,argv,session);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute ../vis.mac");
#endif
    if (ui->IsGUI()) UImanager->ApplyCommand("/control/execute gui.mac");
    UImanager->ApplyCommand("/run/beamOn "+events);
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

