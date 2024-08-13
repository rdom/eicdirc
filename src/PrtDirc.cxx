// simulation software for the Eic DIRC prototype
// author: Roman Dzhygadlo
// contact: r.dzhygadlo at gsi.de

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4PhysListFactory.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"
#include "FTFP_BERT.hh"
#include "G4FastSimulationManagerProcess.hh"

#include "TROOT.h"
#include "TApplication.h"

#include "PrtManager.h"
#include "PrtLutReco.h"
#include "PrtPhysicsList.h"
#include "PrtDetectorConstruction.h"
#include "PrtActionInitialization.h"

#include "time.h"

int main(int argc, char **argv) {
  std::cout << "----------------------------------------------- used args ---\n";
  for (int i = 1; i < argc; i = i + 2) std::cout << argv[i] << "  " << argv[i + 1] << std::endl;

  TApplication theApp("EicDirc", 0, 0);

  TString macro, particle = "pi+", infile = "", lutfile = "", pdffile = "", nnfile = "",
                 outfile = "";
  int events(0), pdgid(0), geometry(1), firstevent(0), runtype(0), study(0), fid(0), verbose(0),
    batchmode(0), physlist(0), pmtLayout(2031), correction(2), field(0), ev(0), radiator(0),
    lensId(3), displayOpt(0);
  double momentum(0), theta(25), phi(0), beamZ(0), trackingres(0.0005), dark_noise(0),
    prismStepX(0), prismStepY, beamX(0), timeSigma(0.1), timeCut(0.5), testVal1(0), testVal2(0),
    testVal3(0);
  long seed = 0;

  for (int i = 1; i < argc; i = i + 2) {
    if (G4String(argv[i]) == "-m") macro = argv[i + 1];
    else if (G4String(argv[i]) == "-seed") seed = atol(argv[i + 1]);
    else if (G4String(argv[i]) == "-o") outfile = argv[i + 1];
    else if (G4String(argv[i]) == "-i") infile = argv[i + 1];
    else if (G4String(argv[i]) == "-u") lutfile = argv[i + 1];
    else if (G4String(argv[i]) == "-pdf") pdffile = argv[i + 1];
    else if (G4String(argv[i]) == "-nn") nnfile = argv[i + 1];
    else if (G4String(argv[i]) == "-field") field = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-g") geometry = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-ev") ev = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-h") radiator = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-theta") theta = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-phi") phi = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-b") batchmode = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-f") firstevent = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-e") events = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-l") lensId = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-x") particle = argv[i + 1];
    else if (G4String(argv[i]) == "-p") momentum = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-w") physlist = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-r") runtype = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-study") study = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-fid") fid = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-c") pmtLayout = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-t1") testVal1 = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-t2") testVal2 = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-t3") testVal3 = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-gsx") prismStepX = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-gsy") prismStepY = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-gz") beamZ = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-gx") beamX = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-timeres") timeSigma = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-timecut") timeCut = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-trackingres") trackingres = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-v") verbose = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-d") displayOpt = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-dn") dark_noise = atof(argv[i + 1]);
    else if (G4String(argv[i]) == "-cor") correction = atoi(argv[i + 1]);
    else {
      G4cerr << "read README.md" << G4endl;
      return 1;
    }
  }

  // default values
  if (runtype < 2) {
    if (geometry == 0 || geometry == 10) {
      beamZ = -600;
    }
    if (geometry == 1 || geometry == 11) {
      beamZ = 0;
    }
    if (geometry == 2 || geometry == 12) {
      beamZ = 10;
    }
    if (momentum < 1e-20) {
      if (particle == "mix_pie") momentum = 1.2;
      else if (particle == "mix_kp") momentum = 12;
      else momentum = 6;
    }
    if (runtype == 1) {
      particle = "opticalphoton";
      momentum = 3.18e-09;
      theta = 180;
    }
  }

  if (outfile == "") {
    outfile = "reco.root";                                          // reconstruction
    if (runtype == 0 || runtype == 10) outfile = "hits.root";       // simulation
    if (runtype == 1 || runtype == 5) outfile = "../data/lut.root"; // lookup table generation
  }

  if (particle == "proton") pdgid = 2212;
  if (particle == "pi+") pdgid = 211;
  if (particle == "pi-") pdgid = -211;
  if (particle == "pi0") pdgid = 111;
  if (particle == "kaon+") pdgid = 321;
  if (particle == "kaon-") pdgid = -321;
  if (particle == "e-") pdgid = 11;
  if (particle == "e+") pdgid = -11;
  if (particle == "mu-") pdgid = 13;
  if (particle == "mu+") pdgid = -13;
  if (particle == "opticalphoton") pdgid = 0;
  if (particle == "mix_pie") pdgid = 10000;
  if (particle == "mix_pimu") pdgid = 10001;
  if (particle == "mix_pik") pdgid = 10003;
  if (particle == "mix_pip") pdgid = 10004;
  if (particle == "mix_kp") pdgid = 10005;

  if (batchmode == 1) gROOT->SetBatch(kTRUE);

  PrtTools t;
  PrtRun *run = t.set_run();

  if (infile != "") run = t.get_run(infile);

  run->setRunType(runtype);
  if (runtype < 2) {
    run->setMomentum(momentum);
    run->setPhysList(physlist);
    run->setField(field);
    run->setGeometry(geometry);
    run->setEv(ev);
    run->setRadiator(radiator);
    run->setLens(lensId);
    run->setPmtLayout(pmtLayout);
    if (pmtLayout == 4) {
      run->setNpmt(6);
      run->setNpix(4096);
    }
    if (pmtLayout == 3) {
      run->setNpmt(1);
      run->setNpix(100000); // max number of pixels
    }
    run->setTrackingResTheta(trackingres);
    run->setTrackingResPhi(trackingres);
    run->setDarkNoise(dark_noise);
    run->setTheta(theta);
    run->setPhi(phi);
    run->setPrismStepX(prismStepX);
    run->setPrismStepY(prismStepY);
    run->setBeamX(beamX);
    run->setBeamZ(beamZ);
    run->setPid(pdgid);
  }

  run->setCorrection(correction);
  run->setTimeSigma(timeSigma);
  run->setTimeCut(timeCut);
  run->setTest1(testVal1);
  run->setTest2(testVal2);
  run->setTest3(testVal3);

  PrtManager::Instance(outfile, run);
  PrtManager::Instance()->setDisplayOpt(displayOpt);  
  std::cout << run->getInfo() << std::endl;

  if (runtype == 2 || runtype == 3 || runtype == 4 || runtype > 19) {
    PrtLutReco *reco = new PrtLutReco(infile, lutfile, pdffile, nnfile, verbose);
    reco->Run(firstevent, events);
    return 0;
  }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4RunManager *runManager = new G4RunManager;

  G4Random::setTheSeed(seed);

  // G4PhysListFactory *physListFactory = new G4PhysListFactory();
  // G4VUserPhysicsList *physicsList = physListFactory->GetReferencePhysList("QGSP_BERT");
  // for(auto v : physListFactory->AvailablePhysLists()) std::cout << "v " << v << std::endl;

  if (physlist == 3) {
    G4VModularPhysicsList *physicsList = new FTFP_BERT;
    G4OpticalPhysics *opticalPhysics = new G4OpticalPhysics();
    auto opticalParams = G4OpticalParameters::Instance();
    opticalParams->SetCerenkovMaxPhotonsPerStep(20);
    opticalParams->SetCerenkovMaxBetaChange(10.0);
    opticalParams->SetCerenkovTrackSecondariesFirst(true);
    physicsList->RegisterPhysics(opticalPhysics);

    // auto fastsim = new G4FastSimulationManagerProcess();
    // // const G4PTblDicIterator* particleIterator = physicsList->GetParticleIterator();
    // auto* particleIterator = G4ParticleTable::GetParticleTable()->GetIterator();
    // particleIterator->reset();
    // while ((*particleIterator)()) {
    //   G4ParticleDefinition *particle = particleIterator->value();
    //   G4ProcessManager *pmanager = particle->GetProcessManager();
    //   pmanager->AddProcess(fastsim, -1, 0, 0);
    // }

    runManager->SetUserInitialization(physicsList);
  } else {
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

  if (macro.Sizeof()) {
    UImanager->ApplyCommand("/control/execute " + macro);
  } else {
    // UImanager->ApplyCommand("/control/execute ../prt.mac");
  }

  if (physlist == 1) {
    UImanager->ApplyCommand("/process/inactivate Decay all");
    UImanager->ApplyCommand("/process/inactivate compt all");
    UImanager->ApplyCommand("/process/inactivate hIoni all");
    UImanager->ApplyCommand("/process/inactivate eIoni all");
    UImanager->ApplyCommand("/process/inactivate muIoni all");
    UImanager->ApplyCommand("/process/inactivate conv all");
  }

  if (batchmode == 1) { // batch mode
    UImanager->ApplyCommand(Form("/run/beamOn %d", events));
  } else { // UI session for interactive mode

    G4UIExecutive *ui = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute ../vis.mac");
    if (ui->IsGUI()) UImanager->ApplyCommand("/control/execute gui.mac");
    UImanager->ApplyCommand(Form("/run/beamOn %d", events));
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;
}
