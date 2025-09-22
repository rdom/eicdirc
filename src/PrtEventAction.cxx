#include "G4TrajectoryContainer.hh"
#include "G4Event.hh"
#include "G4Trajectory.hh"

#include "PrtEventAction.h"

#include <fstream>

void PrtEventAction::EndOfEventAction(const G4Event* event){

  // // export trajectories
  // G4TrajectoryContainer* tc = event->GetTrajectoryContainer();
  // G4int ntr = 0;
  // if (tc) ntr = tc->entries();
  // std::cout << "ntr " << ntr << std::endl;
  
  // if(ntr>0){
  //   std::ofstream  file;
  //   file.open ("trackst.json");
  //   file<< "{\"event\": [";
  //   for (G4int i=0; i<ntr; i++){

  //     G4Trajectory* trj = (G4Trajectory*)((*tc)[i]);
  //     file<<"{\"energy\": "<<trj->GetInitialMomentum().mag() <<" }, \n";
  //     file<<"{\"path\": [";
  //     G4int np = trj->GetPointEntries();
  //     if(trj->GetPoint(np-1)->GetPosition().z()>2700)
  //     for(G4int p=0; p<np; p=p+2){
  // 	//     G4VTrajectoryPoint* tp = trj->GetPoint(p);
  // 	G4ThreeVector pos = trj->GetPoint(p)->GetPosition();
  // 	if(p<np-2) file<< "{\"v\": ["<< pos.x()<<","<<pos.y()<<","<<pos.z()<<"]},";
  // 	else file<< "{\"v\": ["<< pos.x()<<","<<pos.y()<<","<<pos.z()<<"]}";
  //     }
  //     if(i<ntr-1) file<< "]}, \n";
  //     else file<< "]} \n";
  //   }
  //   file<< "]} \n\n";

  //   file.close();
  // }
}
