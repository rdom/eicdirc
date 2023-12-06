#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
// #include "g4root.hh"

#include "Randomize.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4PathFinder.hh"
#include "G4FieldTrack.hh"
#include "G4FieldTrackUpdator.hh"
#include "G4SystemOfUnits.hh"

#include "PrtFastSimModelTracker.h"
#include "PrtManager.h"

PrtFastSimModelTracker::PrtFastSimModelTracker(G4String aModelName, G4Region *aEnvelope)
  : G4VFastSimulationModel(aModelName, aEnvelope) {

  fRun = PrtManager::Instance()->getRun();
  fTrackingRes = fRun->getTrackingResTheta();

  // set tracking resolution
  double mombins[] = {0.75, 1.25, 1.75, 2.5,  3.5,  4.5,  5.5,  6.5,  7.5,  8.5, 9.5,
                      10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5};

  double rtheta_0[] = {5.06794,  3.86015,  2.15975,  1.70208,  1.13201,  0.949914, 0.722017,
                       0.627616, 0.601248, 0.482372, 0.456024, 0.424919, 0.38758,  0.330883,
                       0.314847, 0.298489, 0.276163, 0.282891, 0.287801, 0.251317, 0.232265};
  double rtheta_1[] = {4.30083,  2.28747,  1.74149,  1.2333,   0.884043, 0.679373, 0.583677,
                       0.484924, 0.424464, 0.376458, 0.326366, 0.30053,  0.280166, 0.267392,
                       0.24391,  0.222152, 0.20277,  0.199871, 0.193495, 0.176591, 0.17664};
  double rtheta_2[] = {5.45008,  4.13439,  2.78373,  1.94307,  1.65759,  1.24321,  1.03413,
                       0.844061, 0.764819, 0.697782, 0.591453, 0.537856, 0.521389, 0.434031,
                       0.446836, 0.380808, 0.370561, 0.358088, 0.343058, 0.297783, 0.315344};

  double rphi_0[] = {18.6745,  12.426,   3.19723,  2.67958,  1.97109,  1.82281,  1.44746,
                     1.11069,  1.08232,  1.01937,  0.895427, 0.787919, 0.70771,  0.707664,
                     0.614802, 0.594137, 0.525847, 0.560032, 0.532071, 0.466807, 0.47712};
  double rphi_1[] = {6.54607,  3.47465,  2.92379,  2.04704,  1.48889,  1.17779,  0.912445,
                     0.812745, 0.697406, 0.629205, 0.544928, 0.495619, 0.460813, 0.444274,
                     0.402091, 0.362215, 0.361622, 0.334873, 0.328851, 0.301573, 0.29414};
  double rphi_2[] = {16,       11.1076,  5.4352,   3.55541, 2.73494,  2.07148,  1.98559,
                     1.49008,  1.38194,  1.26842,  1.14616, 1.02384,  1.01611,  0.943076,
                     0.871505, 0.746419, 0.735032, 0.68612, 0.670443, 0.631896, 0.575671};

  grtheta[0] = new TGraph(21, mombins, rtheta_0);
  grtheta[1] = new TGraph(21, mombins, rtheta_1);
  grtheta[2] = new TGraph(21, mombins, rtheta_2);

  grphi[0] = new TGraph(21, mombins, rphi_0);
  grphi[1] = new TGraph(21, mombins, rphi_1);
  grphi[2] = new TGraph(21, mombins, rphi_2);

  // TFile *file = TFile::Open("../data/tracking_resolution_map_220723.root");
  TFile *file = TFile::Open("../data/tracking_resolution_map_271123.root");
  fTrMapTheta = new TH2F();
  fTrMapPhi = new TH2F();
  file->GetObject("tr_theta", fTrMapTheta);
  file->GetObject("tr_phi", fTrMapPhi);
}

PrtFastSimModelTracker::~PrtFastSimModelTracker() {}

G4bool PrtFastSimModelTracker::IsApplicable(const G4ParticleDefinition &aParticleType) {
  return aParticleType.GetPDGCharge() != 0; // Applicable for all charged particles
}

G4bool PrtFastSimModelTracker::ModelTrigger(const G4FastTrack & /*aFastTrack*/) {
  return true; // No kinematical restrictions to apply the parametrisation
}

void PrtFastSimModelTracker::DoIt(const G4FastTrack &aFastTrack, G4FastStep &aFastStep) {

  // Calculate the final position (at the outer boundary of the tracking detector)
  // of the particle with the momentum at the entrance of the tracking detector.

  G4Track track = *aFastTrack.GetPrimaryTrack();

  G4FieldTrack aFieldTrack('0');
  G4FieldTrackUpdator::Update(&aFieldTrack, &track);

  G4double retSafety = -1.0;
  ELimited retStepLimited;
  G4FieldTrack endTrack('a');
  G4double currentMinimumStep = 10.0 * m; // Temporary: change that to sth connected
                                          // to particle momentum.
  G4PathFinder *fPathFinder = G4PathFinder::GetInstance();
  fPathFinder->ComputeStep(aFieldTrack, currentMinimumStep, 0,
                           aFastTrack.GetPrimaryTrack()->GetCurrentStepNumber(), retSafety,
                           retStepLimited, endTrack, aFastTrack.GetPrimaryTrack()->GetVolume());

  // Place the particle at the tracking detector exit
  // (at the place it would reach without the change of its momentum).
  aFastStep.ProposePrimaryTrackFinalPosition(endTrack.GetPosition(), false);

  auto dir = track.GetMomentumDirection();
  double mom = track.GetMomentum().mag();

  // // use 0.5 mrad smearing for PDF 
  // if(G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() > 5000 && fTrackingRes < 0.0005) fTrackingRes = 0.0005;

  double theta = dir.getTheta();
  double dtheta = 0, dphi = 0;
  if (fTrackingRes < 1) {
    dtheta = fTrackingRes;
    dphi = fTrackingRes;
  } else {
    // dtheta = 0.001 * get_res(grtheta, theta, 0.001 * mom);
    // dphi = 0.001 * get_res(grphi, theta, 0.001 * mom);

    theta = 180 - theta / deg;
    dtheta = 0.001 * fTrMapTheta->GetBinContent(fTrMapTheta->FindBin(0.001 * mom, theta));
    dphi = 0.001 * fTrMapPhi->GetBinContent(fTrMapPhi->FindBin(0.001 * mom, theta));    
  }
  
  dir.setTheta(G4RandGauss::shoot(dir.getTheta(), dtheta));
  dir.setPhi(G4RandGauss::shoot(dir.getPhi(), dphi));

  aFastStep.ProposePrimaryTrackFinalMomentumDirection(dir, false);

  fRun->setTrackingResTheta(dtheta);
  fRun->setTrackingResPhi(dphi);

  std::cout << "Smearing at the tracker:  dtheta = " << 1000 * dtheta
            << " mrad, dphi = " << 1000 * dphi << " mrad" << std::endl;
}

double PrtFastSimModelTracker::get_res(TGraph *gg[3], double theta, double mom){

  double x[3], y[3];
  x[0] = 2.0 * atan(exp(-1 * 1.5)) * TMath::RadToDeg();  // 25
  x[1] = 2.0 * atan(exp(0)) * TMath::RadToDeg();         // 90
  x[2] = 2.0 * atan(exp(-1 * -1.5)) * TMath::RadToDeg(); // 155

  y[0] = gg[0]->Eval(mom);
  y[1] = gg[1]->Eval(mom);
  y[2] = gg[2]->Eval(mom);

  TGraph * gr = new TGraph(3,x,y);    
  return gr->Eval(theta);
}

