#if defined(__ACLIC__)
#include "PrtTools.h"
#else
R__LOAD_LIBRARY(../build/libPrt.so)
#endif
using namespace std;

void set_style(TGraph *g, int id){
  int cl[] = {kBlack, kRed, kBlue, kGreen + 1, kCyan - 3, kOrange + 6, kMagenta - 3, kAzure + 2};
  int cm[] = {kBlack,    kRed + 1,    kBlue + 1,    kGreen + 2,
              kCyan + 1, kOrange + 7, kMagenta + 1, kAzure + 2};

  g->SetLineColor(cl[id]);
  g->SetMarkerColor(cm[id]);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.4);  
  gPad->SetGrid();
}

TGraph *get_graph(TString in = "reco.root", TString sx = "theta", TString sy = "",
                  TString scut = "", TString name = "") {

  struct ainfo {
    TString name;
    TString axis;
    double min;
    double max;
    ainfo(TString n = "", TString a = "", double mi = 0, double ma = 0) {
      name = n;
      axis = a;
      min = mi;
      max = ma;
    };
  };

  double x, y, y_err(0),ya[5],ya_err[5];
  ainfo ax, ay;
  TString sall = "spr,nph_gr,nph_ti,cangle,trr";
  TChain ch("reco");
  ch.Add(in);

  vector<ainfo> vainfo = {
    ainfo("theta", "polar angle [deg]", 20, 160),
    ainfo("phi", "azimuthal angle [deg]", 20, 160),
    ainfo("nph_gr", "detected photons [#]", 0, 200),
    ainfo("nph_ti", "detected photons [#]", 0, 200),
    ainfo("spr", "SPR [mrad]", 0, 10),
    ainfo("sep_gr", "separation [s.d.]", 0, 10),
    ainfo("sep_ti", "separation [s.d.]", 0, 10),
  };
  
  for (auto v : vainfo){
    if (v.name.Contains(sx)) ax = v;
    if (v.name.Contains(sy)) ay = v;
  }
  
  ch.SetBranchAddress(sx, &x);

  if(sx == sy) y = x;
  else ch.SetBranchAddress(sy, &ya);
  ch.SetBranchAddress(sy + "_err", &ya_err);
  
  auto gg = new TGraphAsymmErrors();

  int nent = ch.GetEntries();
  int is(-1);
  ch.Draw(">>cutlist", TCut(scut));
  auto elist = (TEventList *)gDirectory->Get("cutlist");

  for (int i = 0; i < elist->GetN(); i++) {
    ch.GetEvent(elist->GetEntry(i));
    if(sall.Contains(sy)){
      y = ya[2];
      y_err = ya_err[2];
    } else {
      y = ya[0];
      y_err = ya_err[0];
    }

    is++;
    if(sx == sy) y = x;

    gg->SetPoint(is, x, y);
    gg->SetPointEYhigh(is, y_err);
    gg->SetPointEYlow(is, y_err);
  }

  gg->Sort();
  gg->GetXaxis()->SetTitle(ax.axis);
  gg->GetYaxis()->SetTitle(ay.axis);
  gg->GetXaxis()->SetRangeUser(ax.min,ax.max);
  gg->GetYaxis()->SetRangeUser(ay.min,ay.max);
  gg->SetTitle(name);

  return gg;
}

void dv_phi(TString in = "res.f0.root"){
  
  gStyle->SetOptTitle(0);  
  PrtTools t;

  vector<array<TString, 5>> wid = {
    {in, "theta", "sep_gr", "barid == 5", "bar 5"},
    {in, "theta", "sep_gr", "barid == 0", "bar 0"},
    {in, "theta", "sep_gr", "barid == 2", "bar 2"},
    {in, "theta", "sep_gr", "barid == 7", "bar 7"},
    {in, "theta", "sep_gr", "barid == 9", "bar 9"},
  };

  TString nid[] = {"sep_gr", "sep_ti", "nph_gr", "nph_ti"};

  for (auto name : nid) {
    t.add_canvas(name + "_f0", 1200, 600);

    TLegend *leg = new TLegend(0.15, 0.7, 0.85, 0.85);
    leg->SetNColumns(3);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->SetBorderSize(0);

    int ind = 0;
    for (auto w : wid) {
      auto g = get_graph(w[0], w[1], name, w[3], w[4]);
      set_style(g, ind);
      leg->AddEntry(g);
      g->Draw((ind++ == 0) ? "APL" : "PL same");
    }
    leg->Draw();

    gPad->Update();
    // 20-160 deg
    auto f2 = new TF1("f2", " 2.0*atan(exp(-x))*TMath::RadToDeg()", 1.7354152, -1.7354152);
    TGaxis *A2 = new TGaxis(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax(),
                            "f2", 510, "-");
    A2->SetTitle("#eta    ");
    A2->SetTitleFont(42);
    A2->SetTitleSize(0.06);
    A2->SetTitleOffset(0.8);
    A2->SetLabelOffset(0.01);
    A2->SetLabelFont(42);
    A2->SetLabelSize(0.05);
    A2->Draw();
  }

  t.save_canvas("data/dv_phi", 1);  
}
