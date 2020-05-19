#include "TInterpreter.h"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"

#include "PrtManager.h"
#include "PrtHit.h"
#include "PrtLutNode.h"

PrtManager * PrtManager::fInstance= NULL;

PrtManager::PrtManager(G4String outfile, G4int runtype)
{
  TString filename  = outfile.c_str();
  fOutName = filename; 
  fOutName = fOutName.Remove(fOutName.Last('.'));
  fRunType = runtype;
  fRootFile = new TFile(filename,"RECREATE");
  
  if(fRunType==0 || fRunType==10){
    fTree = new TTree("data","Prototype hits tree");
    fEvent = new PrtEvent();
    fTree->Branch("PrtEvent", "PrtEvent", &fEvent, 64000, 2);
  }

  if(fRunType==1 || fRunType==5){
    fLut = new TClonesArray("PrtLutNode");
    fTree = new TTree("prtlut","Look-up table for the geometrical reconstruction.");
    fTree->Branch("LUT",&fLut,256000,2); 
    Int_t Nnodes = 100000;
    
    TClonesArray &fLuta = *fLut; 
    for (Long64_t n=0; n<Nnodes; n++) {
      new((fLuta)[n]) PrtLutNode(n);
    }    
  }

  // if(fRunType==2  || fRunType==3 || fRunType==4){
  //   fTree = new TTree("recodata","Reconstructed info for the prototype");
  //   fTrackInfoArray = new TClonesArray("PrtTrackInfo");
  //   fTree->Branch("PrtTrackInfo",&fTrackInfoArray,256000,3); 
  // }
  
  // fHist = new TH1F("id", "name", 100, 0., 100);

  fPhysList = 0;
  fParticle = 0;
  fMomentum = TVector3(2,0,0);
  fGeometry = 1;
  fEvType = 0;
  fAngle = 0;
  fShift = 150;
  fTest1 = 0;
  fTest2 = 0;
  fDispalyOpt = 0;
  fTimeRes = 0;
  fTimeCut = 0;
  fLens = 0;
  fMcpLayout = 0;
  fBeamDimension = 0;
  fZPos = 0.05;
  fRadiator = 1;
  fMix = 0; // no mix
  std::cout<<"PrtManager has been successfully initialized. " <<std::endl;
}

PrtManager* PrtManager::Instance(G4String outfile, G4int runtype){
  if ( !fInstance){
    std::cout<<"Info in (PrtManager::Instance): Making a new instance. "<<std::endl;
    fInstance = new PrtManager(outfile, runtype);
  }
  return fInstance;
}

void PrtManager::AddEvent(PrtEvent event){
  if(fRunType==0 || fRunType==10){
    fEvent = new PrtEvent(event);
    fEvent->SetPhysList(fPhysList);
    fEvent->SetAngle((180*deg-fAngle)/deg);
    // fEvent->SetRadiatorL(fRadiatorL);
    // fEvent->SetRadiatorW(fRadiatorW);
    // fEvent->SetRadiatorH(fRadiatorH);
    fEvent->SetTest1(fTest1);
    fEvent->SetTest2(fTest2);
    fEvent->SetParticle(fParticle);
    fEvent->SetMomentum(fMomentum);
    fEvent->SetGeometry(fGeometry);
    // fEvent->SetEvType(fEvType);
    fEvent->SetLens(fLens);
    fEvent->SetTimeRes(fTimeCut);
  }
}

void PrtManager::AddHit(PrtHit hit){
  if(fRunType==0 || fRunType==10){
    if ( fEvent ){
      fEvent->AddHit(hit);
    }else{
      std::cout<<"Event does not exist. Create it first. "<<std::endl;
    }
  }
  if(fRunType==1 || fRunType==5){

    if(fRunType==5) fMomentum = hit.GetMomentum();
      
    int id = 300*hit.GetMcpId() + hit.GetPixelId();
    ((PrtLutNode*)(fLut->At(id)))->
      AddEntry(id, fMomentum, hit.GetPathInPrizm(),
	       hit.GetNreflectionsInPrizm(),
	       hit.GetLeadTime(),hit.GetGlobalPos(),hit.GetDigiPos());
  } 
}

void PrtManager::AddTrackInfo(PrtTrackInfo trackinfo){
  new ((*fTrackInfoArray)[fTrackInfoArray->GetEntriesFast()]) PrtTrackInfo(trackinfo);
}

void PrtManager::Fill(){
  if(fRunType==0 || fRunType==10){
    fTree->Fill();
    fEvent->Clear();
  }
  if(fRunType==2 || fRunType==3 || fRunType==4){
    fTree->Fill();
    fTrackInfoArray->Clear();
  }
}

void PrtManager::FillLut(){
  if(fRunType==1 || fRunType==5) fTree->Fill();
}
