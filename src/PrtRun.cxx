#include "PrtRun.h"

PrtRun::PrtRun()
  : fShortInfo(""), fName(""), fId(0), fRunType(0), fStudy(0), fMc(0), fPhysList(0), fPid(0),
    fField(0), fGeometry(0), fEv(0), fLens(0), fRadiator(0), fPmtLayout(0), fTrigger(0), fNpmt(0),
    fNpix(0), fTheta(0), fPhi(0), fMomentum(0), fPrismStepX(0), fPrismStepY(0), fBeamX(0),
    fBeamZ(0), fBeamSize(0), fTimeSigma(0), fSimOffset(0), fRadiatorL(0), fRadiatorW(0),
    fRadiatorH(0), fTest1(0), fTest2(0), fTest3(0) {}

void PrtRun::setPmtLayout(Int_t v) {
  if (v == 2030) {
    fNpmt = 12;
    fNpix = 256;
  }
  fPmtLayout = v;
}

TString PrtRun::getInfo() {
  TString info = fShortInfo;
  info += Form("Run type %d \n", fRunType);
  info += Form("Study %d \n", fStudy);
  info += Form("Filed id %d \n", fField);
  info += Form("Geometry id %d \n", fGeometry);
  info += Form("EV id %d \n", fEv);
  info += Form("Pmt layout id %d \n", fPmtLayout);
  info += Form("Lens  id %d \n", fLens);
  info += Form("Radiator  id %d \n", fRadiator);
  info += Form("Physics list %d \n", fPhysList);
  info += Form("Particle  id %d \n", fPid);
  info += Form("Particle momentum %f \n", fMomentum);
  info += Form("Theta %f \n", fTheta);
  info += Form("Phi %f \n", fPhi);
  info += Form("beam X:Z   %f : %f \n", fBeamX, fBeamZ);
  info += Form("Prism step X:Y  %f :%f \n", fPrismStepX, fPrismStepY);
  info += Form("Pmt's time resolution %f \n", fTimeSigma);
  return info;
}

ClassImp(PrtRun)
