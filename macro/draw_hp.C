#if defined(__ACLIC__)
#include "PrtTools.h"
#else
R__LOAD_LIBRARY(../build/libPrt.so)
#endif

void draw_hp(TString infile = "../build/hits.root") {

  PrtTools t(infile);

  while (t.next() && t.i() < 10000) {
    for (auto hit : t.event()->getHits()) {

      int ch = hit.getChannel();
      int pmt = hit.getPmt();
      int pix = hit.getPixel();
      double time = hit.getLeadTime();
      if (t.pid() == 2) t.fill_digi(pmt, pix);
    }
  }
  auto cdigi = t.draw_digi();
  //auto cdigi = t.draw_digi(1, 0, new TCanvas("hp", "hp", 1200, 600));
  t.add_canvas(cdigi);
  t.save_canvas("data/draw_hp", 0);
}
