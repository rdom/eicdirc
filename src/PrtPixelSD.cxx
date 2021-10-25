#include "PrtPixelSD.h"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include <TVector3.h>

#include "PrtEvent.h"
#include "PrtPrizmHit.h"
#include "PrtBarHit.h"

#include "PrtRunAction.h"
#include "PrtManager.h"

PrtPixelSD::PrtPixelSD(const G4String &name, const G4String &hitsCollectionName)
  : G4VSensitiveDetector(name) {
  collectionName.insert(hitsCollectionName);

  fRunType = PrtManager::Instance()->getRun()->getRunType();
  fEvType = PrtManager::Instance()->getRun()->getEv();
  fMcpLayout = PrtManager::Instance()->getRun()->getPmtLayout();
  int npmt = PrtManager::Instance()->getRun()->getNpmt();
  int npix = PrtManager::Instance()->getRun()->getNpix();
  fRadiatorL = PrtManager::Instance()->getRun()->getRadiatorL();
  fRadiatorW = PrtManager::Instance()->getRun()->getRadiatorW();
  fRadiatorH = PrtManager::Instance()->getRun()->getRadiatorH();
  // create MPC map

  std::cout << "npmt * npix " << npmt * npix << std::endl;

  for (int ch = 0; ch < npmt * npix; ch++) {
    int mcp = ch / npix;
    int pix = ch % npix;
    fMap_Mpc[mcp][pix] = ch;
  }
}

PrtPixelSD::~PrtPixelSD() {}

void PrtPixelSD::Initialize(G4HCofThisEvent *hce) {
}

G4bool PrtPixelSD::ProcessHits(G4Step *step, G4TouchableHistory *hist) {
  // // energy deposit
  // G4double edep = step->GetTotalEnergyDeposit();

  // // step length
  // G4double stepLength = 0.;
  // if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
  //   stepLength = step->GetStepLength();
  // }

  // if ( edep==0. && stepLength == 0. ) return false;

  if (step == 0) return false;

  // G4ThreeVector translation = hist->GetTranslation();
  // G4ThreeVector localpos = step->GetPreStepPoint()->GetPhysicalVolume()->GetObjectTranslation();
  G4TouchableHistory *touchable = (G4TouchableHistory *)(step->GetPostStepPoint()->GetTouchable());

  // Get cell id
  G4Track *track = step->GetTrack();
  const G4DynamicParticle *dynParticle = track->GetDynamicParticle();
  G4ParticleDefinition *particle = dynParticle->GetDefinition();
  G4String ParticleName = particle->GetParticleName();

  G4ThreeVector globalpos = step->GetPostStepPoint()->GetPosition();
  G4ThreeVector localpos = touchable->GetHistory()->GetTopTransform().TransformPoint(globalpos);
  G4ThreeVector translation =
    touchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0, 0, 0));
  G4ThreeVector inPrismpos = touchable->GetHistory()->GetTransform(1).TransformPoint(globalpos);
  G4ThreeVector g4mom = track->GetVertexMomentumDirection(); // GetMomentum();
  G4ThreeVector g4pos = track->GetVertexPosition();

  G4ThreeVector localvec = touchable->GetHistory()->GetTopTransform().TransformAxis(g4mom);

  TVector3 globalPos(inPrismpos.x(), inPrismpos.y(), inPrismpos.z());
  // TVector3 globalPos(globalpos.x(),globalpos.y(),globalpos.z());
  TVector3 localPos(localpos.x(), localpos.y(), localpos.z());
  TVector3 digiPos(translation.x(), translation.y(), translation.z());
  TVector3 momentum(g4mom.x(), g4mom.y(), g4mom.z());
  TVector3 position(g4pos.x(), g4pos.y(), g4pos.z());

  // information from prizm
  G4SDManager *fSDM = G4SDManager::GetSDMpointer();
  G4RunManager *fRM = G4RunManager::GetRunManager();
  G4int collectionID = fSDM->GetCollectionID("PrizmHitsCollection");
  const G4Event *currentEvent = fRM->GetCurrentEvent();
  G4HCofThisEvent *HCofEvent = currentEvent->GetHCofThisEvent();
  PrtPrizmHitsCollection *prizmCol = (PrtPrizmHitsCollection *)(HCofEvent->GetHC(collectionID));
  double time = step->GetPreStepPoint()->GetLocalTime();

  Long_t pathId = 0;
  int refl = 0;
  int prizmId = -1;

  double prismtime = 0;
  for (size_t i = 0; i < prizmCol->entries(); i++) {
    PrtPrizmHit *phit = (*prizmCol)[i];
    if (phit->GetTrackID() == track->GetTrackID()) {
      if (fRunType == 5 && phit->GetNormalId() == -5) {
        momentum.SetXYZ(phit->GetPos().x(), phit->GetPos().y(), phit->GetPos().z());
        prismtime = phit->GetEdep();
      }
      if (phit->GetNormalId() > 0) {
        ++refl;
        if (fEvType == 0 && refl == 1) continue;
        pathId = pathId * 10 + phit->GetNormalId();
      }
    }
  }

  // store time in prism
  if (fRunType == 5) time -= prismtime;
  
  // // information from bar
  // G4int collectionID_bar = fSDM->GetCollectionID("BarHitsCollection");
  // PrtBarHitsCollection* barCol = (PrtBarHitsCollection*)(HCofEvent->GetHC(collectionID_bar));
  // std::cout<<"==================== barCol->entries() "<<barCol->entries()<<std::endl;
  // G4ThreeVector bmom1;
  // for (G4int i=0;i<barCol->entries();i++){
  //   PrtBarHit* phit = (*barCol)[i];
  //   G4ThreeVector bmom = phit->GetMom();
  //   if(i==0)  bmom1 = phit->GetMom();
  //   std::cout<<i<<" bmom "<<bmom<< " "<<bmom.angle(bmom1) <<std::endl;
  //   // if(bmom.mag()>1000){
  //   //   PrtHit hit;
  //   //   hit.SetLeadTime(bmom.angle(bmom1)*1000);
  //   //   PrtManager::Instance()->AddHit(hit);
  //   // }
  // }

  PrtHit hit;
  // hit.SetPrizmId(prizmId);
  int mcp = touchable->GetReplicaNumber(1);
  int pix = touchable->GetReplicaNumber(0);
  int ch = fMap_Mpc[mcp][pix];
  hit.setPmt(mcp);
  hit.setPixel(pix);
  hit.setChannel(ch);
  hit.setPosition(position);
  hit.setMomentum(momentum);
  hit.setPathInPrizm(pathId);
  hit.setLeadTime(time); // time since track created 
  double wavelength = 1.2398 / (track->GetMomentum().mag() * 1E6) * 1000;
  hit.setTotTime(wavelength); // set photon wavelength

  // time since event created
  // hit.SetTrailTime(0,step->GetPreStepPoint()->GetGlobalTime()*1000);

  bool transport_efficiency(true);
  if (transport_efficiency) {
    double pi(4 * atan(1));
    double roughness(0.5); // nm
    double angleX = localvec.angle(G4ThreeVector(1, 0, 0));
    double angleY = localvec.angle(G4ThreeVector(0, 1, 0));
    if (angleX > 0.5 * pi) angleX = pi - angleX;
    if (angleY > 0.5 * pi) angleY = pi - angleY;
    double length = track->GetTrackLength() - 400; // 400 - average path in EV
    double lengthx = fabs(length * localvec.x());  // along the bar
    double lengthy = fabs(length * localvec.y());

    int nBouncesX = (int)(lengthx) / fRadiatorH;
    int nBouncesY = (int)(lengthy) / fRadiatorW;

    double ll = wavelength * wavelength;
    double n_quartz = sqrt(1. + (0.696 * ll / (ll - pow(0.068, 2))) +
                           (0.407 * ll / (ll - pow(0.116, 2))) + 0.897 * ll / (ll - pow(9.896, 2)));
    double bounce_probX = 1 - pow(4 * pi * cos(angleX) * roughness * n_quartz / wavelength, 2);
    double bounce_probY = 1 - pow(4 * pi * cos(angleY) * roughness * n_quartz / wavelength, 2);

    double totalProb = pow(bounce_probX, nBouncesX) * pow(bounce_probY, nBouncesY);

    if (G4UniformRand() > totalProb) {
      // std::cout << "photon lost in the radiator. n_bounces = [" << nBouncesX << " " << nBouncesY
      //           << "] with prob = "<<totalProb<<std::endl;
      return true;
    }
  }

  PrtManager::Instance()->addHit(hit, localPos);

  return true;
}

void PrtPixelSD::EndOfEvent(G4HCofThisEvent *) {
  int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  if (eventNumber % 1 == 0 && (fRunType == 0 || fRunType == 5))
    std::cout << " : " << PrtManager::Instance()->getEvent()->getHits().size() << std::endl;
  else if (eventNumber % 1000 == 0 && fRunType != 0)
    std::cout << " : " << PrtManager::Instance()->getEvent()->getHits().size() << std::endl;
  PrtManager::Instance()->fill();
}
