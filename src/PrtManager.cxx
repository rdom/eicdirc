#include "PrtManager.h"
#include "PrtHit.h"
#include "TInterpreter.h"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"

PrtManager * PrtManager::fInstance= NULL;

PrtManager::PrtManager(G4String outfile)
{
  TString filename = outfile;
  fRootFile = new TFile(filename,"RECREATE");
  fTree = new TTree("data","Prototype hits tree");
  fEvent = new PrtEvent();
  fTree->Branch("PrtEvent", "PrtEvent", &fEvent, 64000, 0);

  fHist = new TH1F("id", "name", 100, 0., 100);
  std::cout<<"PrtManager has been successfully initialized. " <<std::endl;
  fPhysList = 0;
  fParticle = 0;
  fMomentum = 0;
  fGeometry = 0;
  fAngle = 0;
  fLens = 0;
}

PrtManager* PrtManager::Instance(G4String outfile){
  if ( !fInstance){
    std::cout<<"Info in (PrtManager::Instance): Making a new instance. "<<std::endl;
    fInstance = new PrtManager(outfile);
  }
  return fInstance;
}
void PrtManager::AddEvent(PrtEvent event){
  fEvent = new PrtEvent(event);
  fEvent->SetPhysList(fPhysList);
  fEvent->SetAngle((180*deg-fAngle)/deg);
  fEvent->SetParticle(fParticle);
  fEvent->SetMomentum(fMomentum);
  fEvent->SetGeometry(fGeometry);
  fEvent->SetLens(fLens);
}


void PrtManager::AddHit(PrtHit hit){
  if ( fEvent ){
    fEvent->AddHit(hit);
  }else{
    std::cout<<"Event does not exist. Create it first. "<<std::endl;
  } 
}

void PrtManager::Fill(){
  std::cout<<"HHHHHHHHHHHHHHHHH   "<< fEvent->GetHitSize() <<std::endl;
  fTree->Fill();
  fEvent = NULL;
}
