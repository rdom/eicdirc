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
#include "TRotation.h"
#include "TGraph.h"
#include <TVirtualFitter.h>
#include <TArc.h>
#include <TLegend.h>

#define eic__sim
#include "../../prttools/prttools.C"


using std::cout;
using std::endl;

TH1F*  fHist1 = new TH1F("Time1","1", 1000,0,20);
TH1F*  fHist2 = new TH1F("Time2","2", 1000,-10,10);
TH1F*  fHistDiff = new TH1F("TimeDiff","calculated time - measured time [ns];entries [#]", 500,-10,10);
TH2F*  fHist3 = new TH2F("Time3","3", 500,5,80, 500,5,60);
TH2F*  fHist4 = new TH2F("Time4","4", 200,-1,1, 200,-1,1);
TH2F*  fHist5 = new TH2F("Time5","5", 200,-1,1, 200,-1,1);


Int_t gg_i(0);
TGraph gg_gr;

// -----   Default constructor   -------------------------------------------
PrtLutReco::PrtLutReco(TString infile, TString lutfile, Int_t verbose){
  fVerbose = verbose;
  fChain = new TChain("data");
  fChain->Add(infile);
  fChain->SetBranchAddress("PrtEvent", &fEvent);

  fFile = new TFile(lutfile);
  fTree=(TTree *) fFile->Get("prtlut") ;
  fLut = new TClonesArray("PrtLutNode");
  fTree->SetBranchAddress("LUT",&fLut); 
  fTree->GetEntry(0);

  fHist = new TH1F("chrenkov_angle_hist","chrenkov_angle_hist", 100,0.7,0.9); //200
  fFit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +[3]",0.35,0.9);
  fSpect = new TSpectrum(10);
  fMethod = PrtManager::Instance()->GetRunType();
  fSavePath="data/reco";

  for(Int_t i=0; i<5; i++){
    fLnDiff[i] = new TH1F(Form("LnDiff_%d",i),  ";ln L(K) - ln L(#pi);entries [#]",100,-500,500);
    fFunc[i] = new TF1(Form("gaus_%d",i),"[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
  }

  cout << "-I- PrtLutReco: Intialization successfull" << endl;
}

// -----   Destructor   ----------------------------------------------------
PrtLutReco::~PrtLutReco(){

}
 
//-------------- Loop over tracks ------------------------------------------
void PrtLutReco::Run(Int_t start, Int_t end){
  TVector3 dird, dir, momInBar(0,0,1),posInBar;
  Double_t cangle,spr,tangle,boxPhi,evtime, bartime, lenz,dirz,luttheta, barHitTime, hitTime;
  Int_t pdgcode, evpointcount(0), tofPid(0),distPid(0),likePid(0);
  Bool_t reflected = kFALSE;
  gStyle->SetOptFit(111);
  
  TVector3 fnX1 = TVector3 (1,0,0);   
  TVector3 fnY1 = TVector3( 0,1,0);
  bool testTrRes = false;
  Double_t angdiv,dtheta,dtphi,mom;
  Int_t nsHits(0),nsEvents(0),studyId(0), nHits(0), ninfit(1);
  
  TString outFile = PrtManager::Instance()->GetOutName()+"_spr.root";
  Double_t theta(0),phi(0), trr(0),  nph(0),
    par1(0), par2(0), par3(0), par4(0), par5(0), par6(0), timeRes(0), test2(0),test3(0),separation(0),likelihood(0);
 
  TFile f(outFile,"recreate");
  TTree tree("dirc","SPR");
  tree.Branch("mom", &mom,"mom/D");
  tree.Branch("tofPid", &tofPid,"tofPid/I");
  tree.Branch("distPid", &distPid,"distPid/I");
  tree.Branch("likePid", &likePid,"likePid/I");
  tree.Branch("spr", &spr,"spr/D");
  tree.Branch("trr", &trr,"trr/D");
  tree.Branch("nph",&nph,"nph/D");
  tree.Branch("cangle",&cangle,"cangle/D");
  tree.Branch("likelihood",&likelihood,"par3/D");
  tree.Branch("separation",&separation,"separation/D");
  tree.Branch("par5",&par5,"par5/D");
  tree.Branch("par6",&par6,"par6/D");
  tree.Branch("timeRes",&timeRes,"timeRes/D");
  tree.Branch("test2",&test2,"test2/D");
  tree.Branch("test3",&test3,"test3/D");
  tree.Branch("theta",&theta,"theta/D");
  tree.Branch("phi",&phi,"phi/D");

  timeRes = PrtManager::Instance()->GetTimeRes();
  Int_t nEvents = fChain->GetEntries();
  if(end==0) end = nEvents;
  
  std::cout<<"Run started for ["<<start<<","<<end <<"]"<<std::endl;

  Int_t ntotal=0;
  for (Int_t ievent=0; ievent<nEvents; ievent++){
    fChain->GetEntry(ievent);
    Int_t nHits = fEvent->GetHitSize();
    ntotal+=nHits;
    mom=fEvent->GetMomentum().Mag()/1000.;    
    tofPid = fEvent->GetParticle();

    if(fMethod==2 && tofPid!=211) continue;
    
    std::cout<<"Event # "<< ievent << " has "<< nHits <<" hits"<<std::endl;
    PrtTrackInfo trackinfo;
    trackinfo.AddInfo(fEvent->PrintInfo()+"\n Basic reco informaion: \n");
    
    Double_t minChangle = 0.35;
    Double_t maxChangle = 0.9;
    trackinfo.AddInfo(Form("Cerenkov angle selection: (%f,%f) \n",minChangle,maxChangle));
    
    TVector3 rotatedmom = fEvent->GetMomentum().Unit();
    Int_t pdg[]={11,13,211,321,2212};
    Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
    Double_t angle1(0), angle2(0),sum1(0),sum2(0), sigma(0.007),range(5*sigma),noise(0.3);

    for(Int_t i=0; i<5; i++){
      fAngle[i] = acos(sqrt(mom*mom + mass[i]*mass[i])/mom/1.4738)+0.002; //1.4738 = 370 = 3.35
      fFunc[i]->SetParameter(0,1);
      fFunc[i]->SetParameter(1,fAngle[i]);
      fFunc[i]->SetParameter(2,sigma);
    }
 
    for(Int_t h=0; h<nHits; h++) {
      PrtPhotonInfo photoninfo;
      fHit = fEvent->GetHit(h);
      hitTime = fHit.GetLeadTime();
      lenz = 2100-fHit.GetPosition().Z();
      dirz = fHit.GetMomentum().Z();

      TVector3 dir0 = fHit.GetMomentum().Unit();
      
      TVector3 cz = TVector3(-rotatedmom.X(),rotatedmom.Y(),rotatedmom.Z());
      TVector3 cd = TVector3(-dir0.X(),dir0.Y(),dir0.Z());
    
      TVector3 unitdir1 = rotatedmom.Unit();
      TVector3 unitdir2 = rotatedmom.Unit();
      cz.RotateUz(unitdir1);
      cd.RotateUz(unitdir2);
      
      Double_t phi0 =  cd.Phi();
      if(dirz<0) reflected = kTRUE;
      else reflected = kFALSE;

      Double_t theta0 = rotatedmom.Angle(dir0);
      fHist5->Fill(theta0*TMath::Sin(phi0),theta0*TMath::Cos(phi0));      

      Int_t sensorId = 100*fHit.GetMcpId()+fHit.GetPixelId();

      PrtLutNode *node = (PrtLutNode*) fLut->At(sensorId);
      Int_t size = node->Entries();
      Bool_t isGoodHit(false);
      
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

	  double criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
	  if(dir.Angle(fnX1) < criticalAngle || dir.Angle(fnY1) < criticalAngle) continue;

	  luttheta = dir.Theta();	
	  if(luttheta > TMath::Pi()/2.) luttheta = TMath::Pi()-luttheta;
	  
	  if(!reflected) bartime = lenz/cos(luttheta)/198.; 
	  else bartime = (2*4200 - lenz)/cos(luttheta)/198.; 
	
	  fHist1->Fill(hitTime);
	  par4=(bartime + evtime)-hitTime;
	  fHistDiff->Fill(par4);
	  fHist3->Fill(fabs((bartime + evtime)),hitTime);
	  if(fabs((bartime + evtime)-hitTime)>timeRes) continue;
	  tangle = rotatedmom.Angle(dir);
	  //if(  tangle>TMath::Pi()/2.) tangle = TMath::Pi()-tangle;

	  if(tangle < 0.6 ||  tangle > 0.9) continue; 
	  PrtAmbiguityInfo ambinfo;
	  ambinfo.SetBarTime(bartime);
	  ambinfo.SetEvTime(evtime);
	  ambinfo.SetCherenkov(tangle);
	  photoninfo.AddAmbiguity(ambinfo);

	  if(tangle > minChangle && tangle < maxChangle){
	    fHist->Fill(tangle);

	    // Double_t phi = rotatedmom.Phi()-dir.Phi();
	    // if(TMath::Abs(tangle-0.82)<0.015)  fHistDiff->Fill(par4);

	    TVector3 rdir = TVector3(-dir.X(),dir.Y(),dir.Z());
	    TVector3 unitdir3 = rotatedmom.Unit();
	    rdir.RotateUz(unitdir3);
	    Double_t phi =  rdir.Phi();
	    fHist4->Fill(tangle*TMath::Sin(phi),tangle*TMath::Cos(phi));
	    gg_gr.SetPoint(gg_i,tangle*TMath::Sin(phi),tangle*TMath::Cos(phi));
	    gg_i++;
	  }
	  isGoodHit=true;
	  
	  sum1 += -TMath::Log(fFunc[2]->Eval(tangle)+noise);
	  sum2 += -TMath::Log(fFunc[3]->Eval(tangle)+noise);
	}
      }

      //if(isGoodHit)
	nsHits++;

      photoninfo.SetHitTime(hitTime);
      photoninfo.SetReflected(reflected);
      photoninfo.SetEvReflections(evpointcount);
      photoninfo.SetSensorId(sensorId);
      photoninfo.SetMcCherenkovInBar(fHit.GetCherenkovMC());      

      trackinfo.AddPhoton(photoninfo);
    }

    Double_t sum = sum1-sum2;
    if(sum!=0){
      if(tofPid==211) fLnDiff[2]->Fill(sum);
      if(tofPid==321) fLnDiff[3]->Fill(sum);
    }

    if(fVerbose==1){
      FindPeak(cangle,spr,fEvent->GetAngle()+0.01);
      nph = nsHits/(Double_t)nsEvents;
      spr = spr*1000;
      trr = spr/sqrt(nph);
      theta = fEvent->GetAngle();
      par3 = fEvent->GetTest2();
      
      // canvasAdd("ff",800,400);
      // //fHist->Scale(1/fHist->GetMaximum());
      // fHist->Draw();
      // fFunc[2]->Draw("same");
      // fFunc[3]->SetLineColor(4);
      // fFunc[3]->Draw("same");
      
      // waitPrimitive("ff");
      // canvasDel("ff");
      // fHist->Reset();
      tree.Fill();
    }
    
    //FindPeak(cangle,spr);
    
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
    if(++nsEvents>=end) break;	
  }

  if(fMethod==2){
    FindPeak(cangle,spr,fEvent->GetAngle()+0.01);
    std::cout<<"nsEvents "<<nsEvents <<"  nsHits "<< nsHits <<std::endl;
    
    nph = nsHits/(Double_t)nsEvents;
    spr = spr*1000;
    trr = spr/sqrt(nph);
    theta = fEvent->GetAngle();
    par3 = fEvent->GetTest2();
    std::cout<<Form("SPR=%2.2F N=%2.2f",spr,nph)<<std::endl; 
  }else{
    //if(!fVerbose) gROOT->SetBatch(1);
    canvasAdd("r_lhood",800,400);
    prt_normalize(fLnDiff[2],fLnDiff[3]);
    fLnDiff[3]->SetLineColor(2);

    TF1 *ff;
    Double_t m1,m2,s1,s2; 
    if(fLnDiff[3]->GetEntries()>10){
      fLnDiff[3]->Fit("gaus","S");
      ff = fLnDiff[3]->GetFunction("gaus");
      m1=ff->GetParameter(1);
      s1=ff->GetParameter(2);
    }
    if(fLnDiff[2]->GetEntries()>10){
      fLnDiff[2]->Fit("gaus","S");
      ff = fLnDiff[2]->GetFunction("gaus");
      m2=ff->GetParameter(1);
      s2=ff->GetParameter(2);
    }
    separation = (fabs(m2-m1))/(0.5*(s1+s2));
    std::cout<<"separation "<< separation <<std::endl;
    
    fLnDiff[3]->SetName(Form("s_%2.2f",separation));
    fLnDiff[3]->Draw();
    fLnDiff[2]->SetLineColor(4);
    fLnDiff[2]->Draw("same");
    canvasSave(1,0);
    //waitPrimitive("r_lhood","w");
    if(fVerbose) gROOT->SetBatch(0);
  }
   
  tree.Fill();
  tree.Write();
  PrtManager::Instance()->Save();
}

Int_t g_num =0;
Bool_t PrtLutReco::FindPeak(Double_t& cherenkovreco, Double_t& spr, Int_t a){
  cherenkovreco=0;
  spr=0;
  //  gStyle->SetCanvasPreferGL(kTRUE);
  
  if(fHist->GetEntries()>20 ){
     gROOT->SetBatch(1);
    Int_t nfound = fSpect->Search(fHist,1,"",0.9); //0.6
    if(nfound>0) cherenkovreco = fSpect->GetPositionX()[0];
    else cherenkovreco =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());

    fFit->SetParameters(100,cherenkovreco,0.005,10);   // peak
    // fFit->SetParameter(1,cherenkovreco);   // peak
    // fFit->SetParameter(2,0.005); // width

    fFit->FixParameter(2,0.009); // width
    fHist->Fit("fgaus","","",cherenkovreco-0.07,cherenkovreco+0.07);
    fFit->ReleaseParameter(2); // width
    fHist->Fit("fgaus","M","",cherenkovreco-0.07,cherenkovreco+0.07);
    cherenkovreco = fFit->GetParameter(1);
    spr = fFit->GetParameter(2); 
    if(fVerbose>2) gROOT->SetBatch(0);
    
    Bool_t storePics(false);
    if(fVerbose>1){
      canvasAdd(Form("tdiff_%d",a),800,400);
      fHistDiff->GetXaxis()->SetTitle("calculated time - measured time [ns]");
      fHistDiff->GetYaxis()->SetTitle("entries [ns]");
      fHistDiff->SetTitle(Form("theta %d", a));
      fHistDiff->Draw();
      
      canvasAdd(Form("tangle_%d",a),800,400);
      fHist->GetXaxis()->SetTitle("#theta_{C} [rad]");
      fHist->GetYaxis()->SetTitle("Entries [#]");
      fHist->SetTitle(Form("theta %d", a));
      fHist->Draw();

      canvasSave(1,0);
      
      if(fVerbose==3){
	TCanvas* c2 = new TCanvas("c2","c2",0,0,800,500);
	// c2->Divide(2,1);
	//c2->cd(1);
	fHist3->GetXaxis()->SetTitle("calculated time [ns]");
	fHist3->GetYaxis()->SetTitle("measured time [ns]");
	fHist3->SetTitle(Form("theta %d", a));

	fHist4->SetStats(0);
	fHist4->GetXaxis()->SetTitle("#theta_{c}sin(#varphi_{c})");
	fHist4->GetYaxis()->SetTitle("#theta_{c}cos(#varphi_{c})");
	fHist4->SetTitle(Form("Calculated from LUT, #theta = %d#circ", a));
	fHist4->Draw("colz");
	Double_t x0(0), y0(0), theta(cherenkovreco);
	FitRing(x0,y0,theta);
	TVector3 corr(x0,y0,1-TMath::Sqrt(x0*x0+y0*y0));
	std::cout<<"Tcorr "<< corr.Theta()*1000<< "  Pcorr "<< corr.Phi() <<std::endl;

	TLegend *leg = new TLegend(0.5,0.7,0.85,0.87);
	//      leg->SetFillColor(0);
	//leg->SetFillColorAlpha(0,0.8);
	//leg->SetFillStyle(0);
	leg->SetFillStyle(4000); 
	leg->SetBorderSize(0);
	leg->AddEntry((TObject*)0,Form("Entries %0.0f",fHist4->GetEntries()),"");
	leg->AddEntry((TObject*)0,Form("#Delta#theta_{c} %f [mrad]",corr.Theta()*1000),"");
	leg->AddEntry((TObject*)0,Form("#Delta#varphi_{c} %f [rad]",corr.Phi()),"");
	leg->Draw();

	TArc *arc = new TArc(x0,y0,theta);
	arc->SetLineColor(kRed);
	arc->SetLineWidth(1);
	arc->SetFillStyle(0);
	arc->Draw();
	gg_i=0;
	gg_gr.Set(0);

	// c2->cd(2);
	// gStyle->SetOptStat(1110); 
	// fHist5->GetXaxis()->SetTitle("#theta_{c}sin(#varphi_{c})");
	// fHist5->GetYaxis()->SetTitle("#theta_{c}cos(#varphi_{c})");
	// fHist5->SetTitle(Form("True from MC, #theta = %d#circ", a));
	// fHist5->Draw("colz");

	c2->Print(Form("spr/tcorr_%d.png", a));
	c2->Modified();
	c2->Update();
	c2->WaitPrimitive();
      }
    }
  }

  if(fVerbose<2) gROOT->SetBatch(0);
  fHist->Reset();
  fHist1->Reset();
  fHist2->Reset();
  fHist3->Reset();
  fHist4->Reset();
  fHistDiff->Reset();
  
  return (cherenkovreco>0 && cherenkovreco<1);
}

void circleFcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
  Int_t np = gg_gr.GetN();
  f = 0;
  Double_t *x = gg_gr.GetX();
  Double_t *y = gg_gr.GetY();
  for (Int_t i=0;i<np;i++) {
    Double_t u = x[i] + par[0];
    Double_t v = y[i] + par[1];
    Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
    f += dr*dr;  
  }
  std::cout<<"fcn  "<< f<<std::endl;
  
}

void circleFcn2(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
  Int_t np = gg_gr.GetN();
  f = 0;
  Double_t *x = gg_gr.GetX();
  Double_t *y = gg_gr.GetY();
  for (Int_t i=0;i<np;i++) {
    Double_t u = x[i] + par[0];
    Double_t v = y[i] + par[1];
    Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
    if(dr>0.07) f += dr*dr; 
    else f += fabs(dr);
  }
}

void PrtLutReco::FitRing(Double_t& x0, Double_t& y0, Double_t& theta){


  TGraph ff_gr;
  Int_t ff_i(0);
  Int_t np = gg_gr.GetN();
  Double_t *x = gg_gr.GetX();
  Double_t *y = gg_gr.GetY();
  for (Int_t i=0;i<np;i++) {
    if( fabs(theta - TMath::Sqrt(x[i]*x[i]+y[i]*y[i]))<0.05) {
      ff_gr.SetPoint(ff_i,x[i],y[i]);
      ff_i++;
    }
  }
  gg_gr = ff_gr;

  //Fit a circle to the graph points
  TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
  TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
  fitter->SetPrecision(0.00000001);
  fitter->SetMaxIterations(1000);

  fitter->SetFCN(circleFcn);
  fitter->SetParameter(0, "x0",   0.03, 0.01, -0.05,0.05);
  fitter->SetParameter(1, "y0",   0, 0.01, -0.05,0.05);
  fitter->SetParameter(2, "R",    theta, 0.01, theta-0.05,theta+0.05);

  //fitter->FixParameter(0);
  //fitter->FixParameter(1);
  fitter->FixParameter(2);
  Double_t arglist[1] = {0};
  fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  // fitter->SetFCN(circleFcn2);
  // fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  x0 = fitter->GetParameter(0);
  y0 = fitter->GetParameter(1);
  theta = fitter->GetParameter(2);
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


