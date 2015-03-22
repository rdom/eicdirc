#include "PrtOpBoundaryProcess.h"

#include "PrtManager.h"

PrtOpBoundaryProcess::PrtOpBoundaryProcess()
  : G4OpBoundaryProcess()
{
  fLensId = PrtManager::Instance()->GetLens();
}

G4VParticleChange* PrtOpBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{  
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  G4VParticleChange* particleChange = G4OpBoundaryProcess::PostStepDoIt(aTrack, aStep); 
  
  // int parentId = aTrack.GetParentID();
  // std::cout<<"parentId   "<<parentId <<std::endl;
  // if(parentId==1) particleChange->ProposeTrackStatus(fStopAndKill);
  
  // ideal focusing
  if(PrtManager::Instance()->GetLens() == 10){
    double endofbar = 4200/2; //1250/2.;
    G4ThreeVector theGlobalPoint1 = pPostStepPoint->GetPosition();
    G4TouchableHistory* touchable = (G4TouchableHistory*)(pPostStepPoint->GetTouchable());
    G4ThreeVector lpoint =  touchable->GetHistory()->GetTransform( 1 ).TransformPoint(theGlobalPoint1);
    if(lpoint.getZ() < endofbar+0.0001 && lpoint.getZ() > endofbar-0.0001){
    
      G4ThreeVector ww  = pPreStepPoint->GetTouchableHandle()->GetHistory()->
	GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0,0,endofbar));
      if(aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()!="wBar") 
	particleChange->ProposeTrackStatus(fStopAndKill);
      else
	aParticleChange.ProposePosition(ww.getX(), ww.getY(),ww.getZ()-0.0005);
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
  }

  if(PrtManager::Instance()->GetRunType() == 1 && pPostStepPoint->GetPosition().z()<pPreStepPoint->GetPosition().z()){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // kill photons outside bar and prizm
  if(GetStatus() == FresnelRefraction 
     && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens1" 
      || aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens2"
      ) //  || aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3"
     &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  return particleChange;

}
