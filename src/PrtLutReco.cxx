// -----------------------------------------
// PrtLutReco.cpp
//
// Created on: 13.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------

#include "PrtLutReco.h"

#include "PrtManager.h"

#include "PrtLutNode.h"
#include "PrtTrackInfo.h"
#include "PrtPhotonInfo.h"
#include "PrtAmbiguityInfo.h"

#include "TCanvas.h"
#include "TMath.h"

using std::cout;
using std::endl;

// -----   Default constructor   -------------------------------------------
PrtLutReco::PrtLutReco(TString infile, TString lutfile){
  fChain = new TChain("data");
  fChain->Add(infile);
  fChain->SetBranchAddress("PrtEvent", &fEvent);

  fFile = new TFile(lutfile);
  fTree=(TTree *) fFile->Get("prtlut") ;
  fLut = new TClonesArray("PrtLutNode");
  fTree->SetBranchAddress("LUT",&fLut); 
  fTree->GetEntry(0);

  fHist = new TH1F("chrenkov_angle_hist","chrenkov_angle_hist", 100,0.7,0.9);
  fFit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.35,0.85);
  fSpect = new TSpectrum(10);
 
  cout << "-I- PrtLutReco: Intialization successfull" << endl;
}

// -----   Destructor   ----------------------------------------------------
PrtLutReco::~PrtLutReco(){

}
 
//-------------- Loop over tracks ------------------------------------------
void PrtLutReco::Run(Int_t start, Int_t end){
 TVector3 dird, dir, momInBar(0,0,1),posInBar;
  Double_t cangle,tangle,boxPhi,evtime, bartime, directz,luttheta, barHitTime, hitTime;
  Int_t pdgcode, evpointcount=0;
  Bool_t reflected = kFALSE;

  TVector3 fnX1 = TVector3 (1,0,0);   
  TVector3 fnY1 = TVector3( 0,1,0);
  bool testTrRes = false;
  Double_t angdiv,dtheta,dtphi;
  
  std::cout<<"Run started " <<std::endl;

  for (Int_t ievent=0; ievent<fChain->GetEntries(); ievent++){
    fChain->GetEntry(ievent);
    std::cout<<"Event # "<< ievent << " has "<< fEvent->GetHitSize() <<" hits"<<std::endl;
    PrtTrackInfo trackinfo;
    trackinfo.AddInfo(fEvent->PrintInfo()+"\n Basic reco informaion: \n");
    
    Double_t minChangle = 0.35;
    Double_t maxChangle = 0.85;
    trackinfo.AddInfo(Form("Cerenkov angle selection: (%f,%f) \n",minChangle,maxChangle));

    for(Int_t h=0; h<fEvent->GetHitSize(); h++) {
      PrtPhotonInfo photoninfo;
      fHit = fEvent->GetHit(h);
      Double_t fAngle =  fEvent->GetAngle()-90;
      hitTime = fHit.GetLeadTime();
      if(hitTime>10000) reflected = kTRUE;
      else  reflected = kFALSE;

      Int_t sensorId = 100*fHit.GetMcpId()+fHit.GetPixelId();
   
      PrtLutNode *node = (PrtLutNode*) fLut->At(sensorId);
      Int_t size = node->Entries();

      TVector3 rotatedmom = momInBar;
      rotatedmom.RotateY(-fAngle/180.*TMath::Pi());
      //std::cout<<"fAngle   "<<fAngle <<std::endl;
      //rotatedmom.Print();
      for(int i=0; i<size; i++){
	dird = node->GetEntry(i);
	evtime = node->GetTime(i);
	for(int u=0; u<4; u++){
	  if(u == 0) dir = dird;
	  if(u == 1) dir.SetXYZ( dird.X(),-dird.Y(),  dird.Z());
	  if(u == 2) dir.SetXYZ( dird.X(), dird.Y(), -dird.Z());
	  if(u == 3) dir.SetXYZ( dird.X(),-dird.Y(), -dird.Z());

	  if(reflected) dir.SetXYZ( -dir.X(), dir.Y(), dir.Z());

	  //if(reflected) dir.RotateY(-2./180.*TMath::Pi());
	
	  // double criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
	  // if(dir.Angle(fnX1) < criticalAngle || dir.Angle(fnY1) < criticalAngle) continue;

	  // luttheta = dir.Theta();	
	  // if(luttheta > TMath::Pi()/2.) luttheta = TMath::Pi()-luttheta;
	  
	  // directz = posInBar.Z()+119;
	  // if(!reflected) bartime = directz/cos(luttheta)/19.8; 
	  // else bartime = ((240 - directz)*2 + directz)/cos(luttheta)/19.8; 
	
	  // if(fabs((bartime + evtime)-(hitTime-barHitTime))>2) continue;

	  tangle = rotatedmom.Angle(dir);
	  //if(  tangle>TMath::Pi()/2.) tangle = TMath::Pi()-tangle;

	  // std::cout<<"tangle  "<<tangle <<std::endl;
	  
	 
	  PrtAmbiguityInfo ambinfo;
	  ambinfo.SetBarTime(bartime);
	  ambinfo.SetEvTime(evtime);
	  ambinfo.SetCherenkov(tangle);
	  photoninfo.AddAmbiguity(ambinfo);
	  if(tangle > minChangle && tangle < maxChangle)
	    fHist->Fill(tangle);
	}
      }

      photoninfo.SetHitTime(hitTime);
      photoninfo.SetReflected(reflected);
      photoninfo.SetEvReflections(evpointcount);
      photoninfo.SetSensorId(sensorId);
      photoninfo.SetMcCherenkovInBar(fHit.GetCherenkovMC());
      

      trackinfo.AddPhoton(photoninfo);
    }

    //Double_t cherenkovreco = FindPeak();
    
    //Int_t pdgreco = FindPdg(fEvent->GetMomentum().Mag(), cherenkovreco);

    // if(fabs(cherenkovreco-cangle)>0.05 && fHist->GetEntries()>20){
    //	if(fVerbose>1){
    	  TCanvas* c = new TCanvas("c","c",0,0,800,1200);
    	  fHist->Draw();
    	  c->Modified();
    	  c->Update();
	  c->WaitPrimitive();
	  //	  	}
	//  }
    fHist->Reset();

    // if(testTrRes) trackinfo.SetMomentum(TVector3(dtheta,dtphi,0)); //track deviation
    trackinfo.SetMcMomentum(fEvent->GetMomentum());
    trackinfo.SetMcMomentumInBar(momInBar);
    trackinfo.SetMcPdg(0);
    trackinfo.SetPdg(0);
    trackinfo.SetAngle(fEvent->GetAngle());
    trackinfo.SetMcCherenkov(cangle);
    //trackinfo.SetCherenkov(cherenkovreco);
    trackinfo.SetMcTimeInBar(barHitTime);
    PrtManager::Instance()->AddTrackInfo(trackinfo);
    PrtManager::Instance()->Fill();
  }
  PrtManager::Instance()->Save();
}

Int_t g_num =0;
Double_t PrtLutReco::FindPeak(){
  Double_t cherenkovreco = -1;
  if(fHist->GetEntries()>20 ){
    Int_t nfound = fSpect->Search(fHist,1,"",0.6);
    Float_t *xpeaks = fSpect->GetPositionX();
    if(nfound>0) cherenkovreco = xpeaks[0];
    fFit->SetParameter(1,cherenkovreco);   // peak
    fFit->SetParameter(2,0.01); // width
    fHist->Fit("fgaus","Q","",cherenkovreco-0.02,cherenkovreco+0.02);
    cherenkovreco = fFit->GetParameter(1);
    if(cherenkovreco<0 || cherenkovreco>1 ) cherenkovreco = 0;
  
    if(fVerbose>1){
      TCanvas* c = new TCanvas("c","c",0,0,800,600);
      fHist->GetXaxis()->SetTitle("#theta_{C}, [rad]");
      fHist->GetYaxis()->SetTitle("Entries, [#]");
      fHist->Draw();
      c->Modified();
      c->Update();
      c->WaitPrimitive();
      c->Print(Form("pic/animpid/animpid_%d.png",g_num++));
    }
  }
  fHist->Reset();

  return cherenkovreco;
}

Int_t PrtLutReco::FindPdg(Double_t mom, Double_t cangle){
  Int_t pdg[]={11,13,211,321,2212};
  Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
  // Int_t pdg[]={211,321,2212};
  // Double_t mass[] = {0.139570,0.49368,0.9382723};
  Double_t tdiff, diff=100;
  Int_t minid=0;
  for(Int_t i=0; i<5; i++){
    tdiff = fabs(cangle - acos(sqrt(mom*mom + mass[i]*mass[i])/mom/1.46907)); //1.46907 - fused silica
    if(tdiff<diff){
      diff = tdiff;
      minid = i;
    }
  }
  return pdg[minid]; 
}


