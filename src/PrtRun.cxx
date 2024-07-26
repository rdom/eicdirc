#include "PrtRun.h"

PrtRun::PrtRun()
  : fShortInfo(""), fName(""), fId(0), fRunType(0), fStudy(0), fMc(0), fPhysList(0), fPid(0),
    fField(0), fGeometry(0), fEv(0), fLens(0), fRadiator(0), fPmtLayout(0), fTrigger(0), fNpmt(0),
    fNpix(0), fCorrection(0), fTheta(0), fPhi(0), fMomentum(0), fPrismStepX(0), fPrismStepY(0),
    fBeamX(0), fBeamZ(0), fBeamSize(0), fTrackingResTheta(0), fTrackingResPhi(0), fTimeSigma(0),
    fTimeCut(0), fSimOffset(0), fRadiatorL(0), fRadiatorW(0), fRadiatorH(0), fDarkNoise(0),
    fTest1(0), fTest2(0), fTest3(0) {}

void PrtRun::setPmtLayout(Int_t v) {
  if (v == 2030) {
    fNpmt = 12;
    fNpix = 256;
  }
  if (v == 2038) {
    fNpmt = 24;
    fNpix = 64;
    v = 2031;
  }
  fPmtLayout = v;
}

TString PrtRun::getInfo() {
  int b = 10;
  TString info =  "------------------------------------------------ run info ---\n";
  info += fShortInfo;
  info += Form("%-50s %*d \n", "Run type",b, fRunType);
  info += Form("%-50s %*d \n", "Study", b, fStudy);
  info += Form("%-50s %*d \n", "Filed id", b, fField);
  info += Form("%-50s %*d \n", "Geometry id", b, fGeometry);
  info += Form("%-50s %*d \n", "EV id", b, fEv);
  info += Form("%-50s %*d \n", "Pmt layout id", b, fPmtLayout);
  info += Form("%-50s %*d \n", "Lens  id ",  b,fLens);
  info += Form("%-50s %*d \n", "Radiator  id", b, fRadiator);
  info += Form("%-50s %*d \n", "Physics list", b, fPhysList);
  info += Form("%-50s %*d \n", "Particle  id", b, fPid);
  info += Form("%-50s %*d \n", "Correction type", b, fCorrection);
  info += Form("%-50s %*.4f \n", "Particle momentum", b, fMomentum);
  info += Form("%-50s %*.4f \n", "Theta", b, fTheta);
  info += Form("%-50s %*.4f \n", "Phi", b, fPhi);
  info += Form("%-50s %*.4f \n", "Pmt's time resolution ", b, fTimeSigma);
  info += Form("%-50s %*.4f \n", "Time cut constant",  b,fTimeCut);
  info += Form("%-50s %*.4f \n", "Dark Noise ", b, fDarkNoise);

  info += Form("%-40s %*.4f%*.4f \n", "beam (x,z)", b, fBeamX, b, fBeamZ);
  info += Form("%-40s %*.4f%*.4f \n", "Prism step (x,y)", b, fPrismStepX, b, fPrismStepY);
  info += Form("%-40s %*.4f%*.4f \n", "Tracking resolution (theta, phi)", b, fTrackingResTheta, b,
               fTrackingResPhi);
  info += "-------------------------------------------------------------\n";
  return info;
}

ClassImp(PrtRun)
