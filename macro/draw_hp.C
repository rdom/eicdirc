#if defined(__ACLIC__)
#include "../src/PrtTools.h"
#else
R__LOAD_LIBRARY(../build/libPrt.so)
#endif

void draw_hp(TString infile = "../build/hits.root") {

  PrtTools t(infile);
  double t1 = t.run()->getTest1();

  while (t.next() && t.i() < 100000) {
    double th = t.event()->getMomentum().Theta();

    for (auto& hit : t.event()->getHits()) {

      // Long_t hpath = hit.getPathInPrizm();
      // TString spath = Form("%ld", hpath);
      // if(!spath.Contains("5") && !spath.Contains("6") && !spath.Contains("7") && !spath.Contains("8")) continue;
       // if(spath.Contains("6")) continue;
      // if(spath.Contains("7")) continue;
      // if(spath.Contains("8")) continue;
    
      int ch = hit.getChannel();
      int pmt = hit.getPmt();
      int pix = hit.getPixel();
      double time = hit.getLeadTime();
      if (t.pid() == 2) t.fill_digi(pmt, pix);
    }
  }

  auto cdigi = t.draw_digi();
  // auto cdigi = t.draw_digi(1, 0, new TCanvas("hp", "hp", 1200, 600));
  cdigi->SetName(Form("hp_%d",(int) t1));
  t.add_canvas(cdigi);
  t.save_canvas("data/draw_hp_l3_w10", 0);
}
