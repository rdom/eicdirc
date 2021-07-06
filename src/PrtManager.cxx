#include "TInterpreter.h"
#include "G4SystemOfUnits.hh"

#include "PrtManager.h"

PrtManager *PrtManager::fInstance = NULL;

PrtManager::PrtManager(TString filename, PrtRun *run) {

  fOutName = filename;
  fRun = run;
  fRunType = fRun->getRunType();
  fEvent = new PrtEvent();
  
  if (fRunType < 2 || fRunType == 5){
    fRootFile = new TFile(filename, "RECREATE");
    fRun->setMc(true);
  }
  
  fRunTree = new TTree("header", "run info");
  fRunTree->Branch("PrtRun", "PrtRun", &fRun, 64000, 2);

  if (fRunType == 0 || fRunType == 5 || fRunType == 6) {
    fTree = new TTree("data", "Prototype hits tree");
    fTree->Branch("PrtEvent", "PrtEvent", &fEvent, 64000, 2);
  }

  if (fRunType == 1 || fRunType == 7 || fRunType == 11) {
    fLut = new TClonesArray("PrtLutNode");
    fTree = new TTree("prtlut", "Look-up table for the geometrical reconstruction.");
    fTree->Branch("LUT", &fLut, 256000, 2);
    TClonesArray &fLuta = *fLut;

    for (Long64_t i = 0; i < fRun->getNpmt()*fRun->getNpix(); i++) {
      new ((fLuta)[i]) PrtLutNode(i);
    }        
  }

  fnX1 = TVector3(1, 0, 0);
  fnY1 = TVector3(0, 1, 0);
  fCriticalAngle = asin(1.00028 / 1.47125);

  std::cout << "PrtManager has been successfully initialized. " << std::endl;
}

PrtManager *PrtManager::Instance(TString outfile, PrtRun *run) {
  if (!fInstance) {
    std::cout << "Info in (PrtManager::Instance): Making a new instance. " << std::endl;
    fInstance = new PrtManager(outfile, run);
  }
  return fInstance;
}

void PrtManager::addEvent(PrtEvent event) {
  if (fRunType == 0 || fRunType == 5 || fRunType == 6) {
    fEvent = new PrtEvent(event);
  }
}

void PrtManager::addHit(PrtHit hit) {
  fEvent->addHit(hit);
}

void PrtManager::addHit(PrtHit hit, TVector3 localpos, TVector3 vertex) {
  if (fRunType == 0 || fRunType == 5 || fRunType == 6) {
    fEvent->addHit(hit);
  } else if (fRunType == 1 || fRunType == 7 || fRunType == 11) {
    if (fMomentum.Angle(fnX1) > fCriticalAngle && fMomentum.Angle(fnY1) > fCriticalAngle) {
      int ch = hit.getChannel();
      double time = hit.getLeadTime();
      // if (fRunType == 11) time -= fTime;
      ((PrtLutNode *)(fLut->At(ch)))
        ->AddEntry(ch, fMomentum, hit.getPathInPrizm(), 0, time, localpos, vertex);
    }
  }
}

void PrtManager::fill() {
  if (fRunType == 0 || fRunType == 5 || fRunType == 6) {
    fTree->Fill();
    fEvent->Clear();
  }
}

void PrtManager::save(){
  if (fRootFile) {
    fRunTree->Fill();
    fRootFile->Write();
  }
}

void PrtManager::fillLut() {
  if (fRunType == 1 || fRunType == 7 || fRunType == 11) fTree->Fill();
}
