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

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "TApplication.h"

#include "PrtManager.h"
#include "PrtLutReco.h"
#include "../../prttools/PrtTools.h"

#include "G4PhysListFactory.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"
#include "FTFP_BERT.hh"



namespace {
void PrintUsage() {
  G4cerr << " read README.md " << G4endl;
}
} // namespace

int main(int argc, char **argv) {
  for (G4int i = 1; i < argc; i = i + 2) std::cout << argv[i] << "  " << argv[i + 1] << std::endl;

  // Evaluate arguments
  if (argc > 50) {
    PrintUsage();
    return 1;
  }

#ifdef G4MULTITHREADED
  G4int nThreads = 1;
#endif
  TApplication theApp("EicDirc", 0, 0);

  G4String macro, events, field, ev, radiator, physlist, session, geomTheta, geomPhi, displayOpt,
    lensId, particle = "", testVal1, testVal2, testVal3, prismStepX, prismStepY,
            beamX, timeSigma, timeCut, mcpLayout;
  TString infile = "", lutfile = "", pdffile = "", nnfile = "", outfile = "";
  int geometry(-1), firstevent(0), runtype(0), study(0), fid(0), verbose(0), batchmode(0), physlistid(0);
  double momentum(-1), beamZ(20000), trackingres(-1);

  G4long myseed = 0;
  for (G4int i = 1; i < argc; i = i + 2) {
    if (G4String(argv[i]) == "-m") macro = argv[i + 1];
    else if (G4String(argv[i]) == "-seed") myseed = atol(argv[i + 1]);
    else if (G4String(argv[i]) == "-o") outfile = argv[i + 1];
    else if (G4String(argv[i]) == "-i") infile = argv[i + 1];
    else if (G4String(argv[i]) == "-u") lutfile = argv[i + 1];
    else if (G4String(argv[i]) == "-pdf") pdffile = argv[i + 1];
    else if (G4String(argv[i]) == "-nn") nnfile = argv[i + 1];
    else if (G4String(argv[i]) == "-field") field = argv[i + 1];
    else if (G4String(argv[i]) == "-g") geometry = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-ev") ev = argv[i + 1];
    else if (G4String(argv[i]) == "-h") radiator = argv[i + 1];
    else if (G4String(argv[i]) == "-theta") geomTheta = argv[i + 1];
    else if (G4String(argv[i]) == "-phi") geomPhi = argv[i + 1];
    else if (G4String(argv[i]) == "-b") batchmode = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-f") firstevent = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-e") events = argv[i + 1];
    else if (G4String(argv[i]) == "-l") lensId = argv[i + 1];
    else if (G4String(argv[i]) == "-x") particle = argv[i + 1];
    else if (G4String(argv[i]) == "-p") momentum =  atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-w") physlist = argv[i + 1];
    else if (G4String(argv[i]) == "-r") runtype = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-study") study = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-fid") fid = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-c") mcpLayout = argv[i + 1];
    else if (G4String(argv[i]) == "-t1") testVal1 = argv[i + 1];
    else if (G4String(argv[i]) == "-t2") testVal2 = argv[i + 1];
    else if (G4String(argv[i]) == "-t3") testVal3 = argv[i + 1];
    else if (G4String(argv[i]) == "-gsx") prismStepX = argv[i + 1];
    else if (G4String(argv[i]) == "-gsy") prismStepY = argv[i + 1];
    else if (G4String(argv[i]) == "-gz") beamZ = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-gx") beamX = argv[i + 1];
    else if (G4String(argv[i]) == "-timeres") timeSigma = argv[i + 1];
    else if (G4String(argv[i]) == "-timecut") timeCut = argv[i + 1];
    else if (G4String(argv[i]) == "-trackingres") trackingres = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-v") verbose = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-d") displayOpt = argv[i + 1];

#ifdef G4MULTITHREADED
    else if (G4String(argv[i]) == "-t") {
      nThreads = G4UIcommand::ConvertToInt(argv[i + 1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }

  // default values
  if (runtype == 0 || runtype == 1) {
    if (geometry == 0 || geometry == 10) {
      beamZ = -600;
    }
    if (geometry == 1 || geometry == 11) {
      beamZ = -0;
    }
    if (geometry == 2 || geometry == 12) {
      beamZ = 10;
    }
    if (momentum < 0) {
      if (particle == "mix_pie") momentum = 1.2;
      else if (particle == "mix_kp") momentum = 12;
      else momentum = 6;
      if (trackingres < 0) trackingres = 0.0005;
    }
    if (runtype == 1) {
      particle = "opticalphoton";
      momentum = 3.18e-09;
      geomTheta = "180";
      geomPhi = "0";
    }
  }

  // if (runtype == 2 && trackingres < 0) {      
  //   if (particle == "mix_pie") trackingres = 0.0022;
  //   else trackingres = 0.0005;
  // }

  if (runtype == 2 && trackingres > 50) {
    trackingres = -1;
  }

  if (outfile == "" && (runtype == 0 || runtype == 10)) outfile = "hits.root"; // simulation
  if (outfile == "" && (runtype == 1 || runtype == 5))
    outfile = "../data/lut.root"; // lookup table generation
  if (outfile == "" && (runtype == 2 || runtype == 3 || runtype == 4))
    outfile = "reco.root"; // reconstruction

  if (batchmode == 1) gROOT->SetBatch(kTRUE);
  if (!events.size()) events = "0";

  PrtTools t;
  PrtRun *run = t.find_run(study, fid);

  if (runtype == 2 || runtype == 3 || runtype == 4) {
    if (infile == "") {
      infile = t.get_inpath();
      std::cout << "--- infile  " << infile << std::endl;
    }
    run = t.get_run(infile);
    if (outfile == "") {
      outfile = t.get_outpath();
      if (run->getStudy() == 0) outfile = "reco.root";
      std::cout << "--- outfile  " << outfile << std::endl;
    }
  }

  run->setRunType(runtype);

  if (momentum > -1) run->setMomentum(momentum);
  if (physlist.size()) {
    physlistid = atoi(physlist);
    run->setPhysList(physlistid);
  }
  if (field.size()) run->setField(atoi(field));
  if (geometry > -1) run->setGeometry(geometry);
  if (ev.size()) run->setEv(atoi(ev));
  if (radiator.size()) run->setRadiator(atoi(radiator));
  if (lensId.size()) run->setLens(atoi(lensId));
  if (mcpLayout.size()) run->setPmtLayout(atoi(mcpLayout));
  if (trackingres > -1) {
    run->setTrackingResTheta(trackingres);
    run->setTrackingResPhi(trackingres);
  }
  if (testVal1.size()) run->setTest1(atof(testVal1));
  if (testVal2.size()) run->setTest2(atof(testVal2));
  if (testVal3.size()) run->setTest3(atof(testVal3));
  if (geomTheta.size()) run->setTheta(atof(geomTheta));
  if (geomPhi.size()) run->setPhi(atof(geomPhi));
  if (prismStepX.size()) run->setPrismStepX(atof(prismStepX));
  if (prismStepY.size()) run->setPrismStepY(atof(prismStepY));
  if (beamX.size()) run->setBeamX(atof(beamX));
  if (beamZ < 10000) run->setBeamZ(beamZ);
  if (timeSigma.size()) run->setTimeSigma(atof(timeSigma));
  if (timeCut.size()) run->setTimeCut(atof(timeCut));
  if (displayOpt.size()) PrtManager::Instance()->setDisplayOpt(atoi(displayOpt));

  PrtManager::Instance(outfile, run);

  if (particle.size()) {
    int pdgid = 0;
    if (particle == "proton") pdgid = 2212;
    if (particle == "pi+") pdgid = 211;
    if (particle == "pi0") pdgid = 111;
    if (particle == "kaon+") pdgid = 321;
    if (particle == "kaon-") pdgid = -321;
    if (particle == "mu-") pdgid = -13;
    if (particle == "e-") pdgid = -11;
    if (particle == "mu+") pdgid = 13;
    if (particle == "e+") pdgid = 11;
    if (particle == "opticalphoton") pdgid = 0;
    if (particle == "mix_pie") pdgid = 10000;
    if (particle == "mix_pimu") pdgid = 10001;
    if (particle == "mix_pik") pdgid = 10003;
    if (particle == "mix_pip") pdgid = 10004;
    if (particle == "mix_kp") pdgid = 10005;
    PrtManager::Instance()->getRun()->setPid(pdgid);
  }

  std::cout << "=== Run info:  " << std::endl << run->getInfo() << std::endl;

  if (runtype == 2 || runtype == 3 || runtype == 4) {
    PrtLutReco *reco = new PrtLutReco(infile, lutfile, pdffile, nnfile, verbose);
    reco->Run(firstevent, atoi(events));
    return 0;
  }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
// #ifdef G4MULTITHREADED
//   G4MTRunManager *runManager = new G4MTRunManager;
//   if (nThreads > 0) runManager->SetNumberOfThreads(nThreads);
// #else
  G4RunManager *runManager = new G4RunManager;
// #endif
 
  G4Random::setTheSeed(myseed);

  // G4PhysListFactory *physListFactory = new G4PhysListFactory();
  // G4VUserPhysicsList *physicsList = physListFactory->GetReferencePhysList("QGSP_BERT");
  // for(auto v : physListFactory->AvailablePhysLists()) std::cout << "v " << v << std::endl;

  if (physlistid == 3) {
    G4VModularPhysicsList *physicsList = new FTFP_BERT;
    G4OpticalPhysics *opticalPhysics = new G4OpticalPhysics();
    auto opticalParams = G4OpticalParameters::Instance();
    opticalParams->SetCerenkovMaxPhotonsPerStep(20);
    opticalParams->SetCerenkovMaxBetaChange(10.0);
    opticalParams->SetCerenkovTrackSecondariesFirst(true);
    physicsList->RegisterPhysics(opticalPhysics);
    runManager->SetUserInitialization(physicsList);
  }else{
    runManager->SetUserInitialization(new PrtPhysicsList());
  }
  
  runManager->SetUserInitialization(new PrtDetectorConstruction());
  runManager->SetUserInitialization(new PrtActionInitialization());
  runManager->Initialize();

 
  // Initialize visualization
  G4VisManager *visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager *UImanager = G4UImanager::GetUIpointer();

  if (macro.size()) {
    UImanager->ApplyCommand("/control/execute " + macro);
  } else {
    // UImanager->ApplyCommand("/control/execute ../prt.mac");
  }

  if (batchmode == 1) { // batch mode
    UImanager->ApplyCommand("/run/beamOn " + events);    
  } else { // UI session for interactive mode

    G4UIExecutive *ui = new G4UIExecutive(argc, argv, session);
    UImanager->ApplyCommand("/control/execute ../vis.mac");
    if (ui->IsGUI()) UImanager->ApplyCommand("/control/execute gui.mac");
    UImanager->ApplyCommand("/run/beamOn " + events);
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;
}
