#if defined(__ACLIC__)
#include "PrtTools.h"
#else
R__LOAD_LIBRARY(../build/libPrt.so)
#endif

void draw_digi(TString infile = "../build/hits.root") {

  PrtTools t(infile);

  TH1D *hist_mcp_num_default = new TH1D("mcp_num_default", ";MCP id; Number of hits", 24, 0, 24);
  TH1D *hist_mcp_num_new = new TH1D("mcp_num_new", ";MCP id; Number of hits", 24, 0, 24);

  TH1D *hist_pixel_num_default =
    new TH1D("pixel_num_default", ";Pixel id; Number of hits", 256, 0, 256);
  TH1D *hist_pixel_num_new = new TH1D("pixel_num_new", ";Pixel id; Number of hits", 256, 0, 256);

  std::cout << std::setprecision(12);

  double prism_width = 350;//1.35;
  double lens_height = 50.0;
  double prism_x_edge = -lens_height / 2;
  double prism_y_edge = -prism_width / 2;
  double NCol = 4, NRow = 6;
  double prism_x_dim = lens_height + (300 * tan(32 * TMath::DegToRad()));
  double prism_y_dim = prism_width;
  double MCP_total_dim = 57.0;
  double MCP_active_dim = 53.0;
  double pixel_dim = MCP_active_dim / 16.0;

  double gap_x = (prism_x_dim - NCol * MCP_total_dim) / (NCol + 1);
  double gap_y = (prism_y_dim - NRow * MCP_total_dim) / (NRow + 1);

  double mcp_active_x_edge = prism_x_edge + gap_x + 2;
  double mcp_active_y_edge = prism_y_edge + gap_y + 2;

  while (t.next() && t.i() < 10000) {
    for (auto hit : t.event()->getHits()) {

      int ch = hit.getChannel();
      int pmt = hit.getPmt();
      int pix = hit.getPixel();
      double time = hit.getLeadTime();

      if (t.pid() == 2) {
        double hit_pos_x = hit.getPosition().X();
        double hit_pos_y = hit.getPosition().Y();
     
        int mx = (hit_pos_x - (prism_x_edge + gap_x)) / (MCP_total_dim + gap_x);
        int my = (hit_pos_y - (prism_y_edge + gap_y)) / (MCP_total_dim + gap_y);
        int mcp_num = (6 * mx) + my;

        double hit_mcp_x_edge = mcp_active_x_edge + (MCP_active_dim + gap_x + 4) * mx;
        double hit_mcp_y_edge = mcp_active_y_edge + (MCP_active_dim + gap_y + 4) * my;

        int pixel_x = (hit_pos_x - hit_mcp_x_edge) / pixel_dim;
        int pixel_y = (hit_pos_y - hit_mcp_y_edge) / pixel_dim;

        if (((hit_pos_x - hit_mcp_x_edge) > 0) &&
            (fmod(hit_pos_x - hit_mcp_x_edge, pixel_dim) == 0))
          pixel_x -= 1;
        if (((hit_pos_y - hit_mcp_y_edge) > 0) &&
            (fmod(hit_pos_y - hit_mcp_y_edge, pixel_dim) == 0))
          pixel_y -= 1;

        int pixel_num = (16 * pixel_x) + pixel_y;

        hist_mcp_num_default->Fill(pmt);
        hist_mcp_num_new->Fill(mcp_num);

        hist_pixel_num_default->Fill(pix);
        hist_pixel_num_new->Fill(pixel_num);
      }

      if (t.pid() == 2) t.fill_digi(pmt, pix);
    }
  }

  t.add_canvas("pix");
  hist_pixel_num_default->SetLineColor(1);
  hist_pixel_num_default->Draw();
  hist_pixel_num_new->SetLineColor(2);
  hist_pixel_num_new->Draw("same");
  
  auto cdigi = t.draw_digi(0, 0);
  t.add_canvas(cdigi);
  t.save_canvas("data/draw_digi", 0);
}
