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

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"

using std::cout;
using std::endl;

TH1F*  fHist1 = new TH1F("Time1","1", 1000,0,20);
TH1F*  fHist2 = new TH1F("Time2","2", 1000,0,20);
TH2F*  fHist3 = new TH2F("Time3","3", 500,5,80, 500,5,60);

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

  fHist = new TH1F("chrenkov_angle_hist","chrenkov_angle_hist", 200,0.7,0.9);
  fFit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +[3]",0.35,0.9);
  fSpect = new TSpectrum(10);
 
  cout << "-I- PrtLutReco: Intialization successfull" << endl;
}

// -----   Destructor   ----------------------------------------------------
PrtLutReco::~PrtLutReco(){

}
 
//-------------- Loop over tracks ------------------------------------------
void PrtLutReco::Run(Int_t start, Int_t end){
 TVector3 dird, dir, momInBar(0,0,1),posInBar;
 Double_t cangle,spr,tangle,boxPhi,evtime, bartime, lenz,dirz,luttheta, barHitTime, hitTime;
  Int_t pdgcode, evpointcount=0;
  Bool_t reflected = kFALSE;
  gStyle->SetOptFit(111);

  TVector3 fnX1 = TVector3 (1,0,0);   
  TVector3 fnY1 = TVector3( 0,1,0);
  bool testTrRes = false;
  Double_t angdiv,dtheta,dtphi;

  TString outFile = PrtManager::Instance()->GetOutName()+"_spr.root";
  Double_t theta(0),phi(0), trr(0),  nph(0),
    par1(0), par2(0), par3(0), par4(0), par5(0), par6(0), par7(0), par8(0);

  TFile f(outFile,"recreate");
  TTree tree("dirc","SPR");
  tree.Branch("spr", &spr,"spr/D");
  tree.Branch("trr", &trr,"trr/D");
  tree.Branch("nph",&nph,"nph/D");
  tree.Branch("cangle",&cangle,"cangle/D");
  tree.Branch("par3",&par3,"par3/D");
  tree.Branch("par4",&par4,"par4/D");
  tree.Branch("par5",&par5,"par5/D");
  tree.Branch("par6",&par6,"par6/D");
  tree.Branch("par7",&par7,"par7/D");
  tree.Branch("par8",&par8,"par8/D");
  tree.Branch("theta",&theta,"theta/D");
  tree.Branch("phi",&phi,"phi/D");
  
  std::cout<<"Run started " <<std::endl;
  Int_t ntotal=0;
  Int_t nEvents = fChain->GetEntries();
  for (Int_t ievent=0; ievent<nEvents; ievent++){
    fChain->GetEntry(ievent);
    Int_t nHits = fEvent->GetHitSize();
    ntotal+=nHits;
    std::cout<<"Event # "<< ievent << " has "<< nHits <<" hits"<<std::endl;
    PrtTrackInfo trackinfo;
    trackinfo.AddInfo(fEvent->PrintInfo()+"\n Basic reco informaion: \n");
    
    Double_t minChangle = 0.35;
    Double_t maxChangle = 0.9;
    trackinfo.AddInfo(Form("Cerenkov angle selection: (%f,%f) \n",minChangle,maxChangle));

    TVector3 rotatedmom = fEvent->GetMomentum().Unit();
    
    // lenz= 4200/2. - 1000/tan((180-fEvent->GetAngle())/180.*TMath::Pi());
    // std::cout<<"lenz "<<lenz <<std::endl;
 
    for(Int_t h=0; h<nHits; h++) {
      PrtPhotonInfo photoninfo;
      fHit = fEvent->GetHit(h);
      hitTime = fHit.GetLeadTime();
      lenz = 2100-fHit.GetPosition().Z();
      dirz = fHit.GetMomentum().Z();
      
      if(dirz<0) reflected = kTRUE;
      else reflected = kFALSE;

      Int_t sensorId = 100*fHit.GetMcpId()+fHit.GetPixelId();
   
      PrtLutNode *node = (PrtLutNode*) fLut->At(sensorId);
      Int_t size = node->Entries();

      // Double_t fAngle =  fEvent->GetAngle()-90;
      // TVector3 rotatedmom = momInBar;
      // rotatedmom.RotateY(-fAngle/180.*TMath::Pi());
      // std::cout<<"fAngle   "<<fAngle <<std::endl;
      // rotatedmom.Print();
    
      for(int i=0; i<size; i++){
	dird = node->GetEntry(i);
	evtime = node->GetTime(i);
	for(int u=0; u<4; u++){
	  if(u == 0) dir = dird;
	  if(u == 1) dir.SetXYZ( dird.X(),-dird.Y(), dird.Z());
	  if(u == 2) dir.SetXYZ(-dird.X(), dird.Y(), dird.Z());
	  if(u == 3) dir.SetXYZ(-dird.X(),-dird.Y(), dird.Z());
	  if(reflected) dir.SetXYZ( dir.X(), dir.Y(),-dir.Z());

	  //if(reflected) dir.RotateY(-2./180.*TMath::Pi());
	
	  double criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
	  //if(dir.Angle(fnX1) < criticalAngle || dir.Angle(fnY1) < criticalAngle) continue;

	  luttheta = dir.Theta();	
	  if(luttheta > TMath::Pi()/2.) luttheta = TMath::Pi()-luttheta;
	  
	  if(!reflected) bartime = lenz/cos(luttheta)/198.; 
	  else bartime = (2*4200 - lenz)/cos(luttheta)/198.; 
	
	  fHist1->Fill(hitTime);
	  fHist2->Fill(bartime + evtime);
 
	  // if(fabs((bartime + evtime)-hitTime)>0.5) continue;
	  fHist3->Fill(fabs((bartime + evtime)),hitTime);
	  tangle = rotatedmom.Angle(dir);
	  if(  tangle>TMath::Pi()/2.) tangle = TMath::Pi()-tangle;

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

    //FindPeak(cangle,spr);
    std::cout<<"RES   "<<spr*1000 << "   N "<<nHits << "  "<<spr/sqrt(nHits)*1000<<std::endl;
    
    
    //Int_t pdgreco = FindPdg(fEvent->GetMomentum().Mag(), cherenkovreco);

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
 
  FindPeak(cangle,spr,fEvent->GetAngle()+0.01);
  Double_t aEvents = ntotal/(Double_t)nEvents;

  nph = ntotal/(Double_t)nEvents;
  spr = spr*1000;
  trr = spr/sqrt(aEvents);
  theta = fEvent->GetAngle();
  par3 = fEvent->GetTest();
  
  std::cout<<"RES   "<<spr << "   N "<< nph << " trr  "<<trr<<std::endl;
    
  tree.Fill();
  tree.Write();

  PrtManager::Instance()->Save();
}

Int_t g_num =0;
Bool_t PrtLutReco::FindPeak(Double_t& cherenkovreco, Double_t& spr, Int_t a){
  cherenkovreco=0;
  spr=0;

  if(fHist->GetEntries()>20 ){
    gROOT->SetBatch(1);
    Int_t nfound = fSpect->Search(fHist,1,"",0.6);
    Float_t *xpeaks = fSpect->GetPositionX();
    if(nfound>0) cherenkovreco = xpeaks[0];
    fFit->SetParameters(100,cherenkovreco,0.005,10);   // peak
    // fFit->SetParameter(1,cherenkovreco);   // peak
    // fFit->SetParameter(2,0.005); // width

    fFit->FixParameter(2,0.009); // width
    fHist->Fit("fgaus","","",cherenkovreco-0.07,cherenkovreco+0.07);
    fFit->ReleaseParameter(2); // width
    fHist->Fit("fgaus","M","",cherenkovreco-0.07,cherenkovreco+0.07);
    cherenkovreco = fFit->GetParameter(1);
    spr = fFit->GetParameter(2);
    //gROOT->SetBatch(0);
    
    Int_t fVerbose=1;
    if(fVerbose>0){
      TCanvas* c = new TCanvas("c","c",0,0,800,500);
      fHist->GetXaxis()->SetTitle("#theta_{C} [rad]");
      fHist->GetYaxis()->SetTitle("Entries [#]");
      fHist->SetTitle(Form("theta %d", a));
      fHist->Draw();
      // fHist1->SetLineColor(2);
      // fHist1->Draw();
      // fHist2->Draw("same");

      c->Modified();
      c->Update();
      c->Print(Form("spr/tangle_%d.png", a));
      c->WaitPrimitive();

      TCanvas* c2 = new TCanvas("c2","c2",0,0,800,600);
      fHist3->GetXaxis()->SetTitle("calculated time [ns]");
      fHist3->GetYaxis()->SetTitle("measured time [ns]");
      fHist3->SetTitle(Form("theta %d", a));
      fHist3->Draw("colz");
      c2->Print(Form("spr/tcorr_%d.png", a));
      c2->Modified();
      c2->Update();
      c2->WaitPrimitive();
    
    }
  }

  gROOT->SetBatch(0);
  fHist->Reset();
  fHist1->Reset();
  fHist2->Reset();
  fHist3->Reset();

  return (cherenkovreco>0 && cherenkovreco<1);
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


