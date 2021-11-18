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
      int pix = hit.getpixel();
      double time = hit.getleadtime();

      if (t.pid() == 2) t.fill_digi(pmt, pix);
    }
  }

  auto cdigi = t.draw_digi(0, 0);
  t.add_canvas(cdigi);
  t.save_canvas("data/draw_hp", 0);
}
