#define eic__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

void draw_hp(TString infile="../build/hits.root",TString id = ""){

  if(!prt_init(infile,2)) return;
  for (int ievent=0; ievent<prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    int pdg = prt_event->GetParticle();
    for(auto hit : prt_event->GetHits()){
      int mcp = hit.GetMcpId();
      int pix = hit.GetPixelId()-1;
      double time = hit.GetLeadTime();
      // if(pdg==211) prt_hdigi[mcp]->Fill(pix/16, pix%16);
      if(pdg==211) prt_hdigi[mcp]->Fill(pix%16, pix/16);
    }
  }

  TGaxis::SetMaxDigits(3);
    
  auto cdigi = prt_drawDigi(2030); //2031 //2032
  cdigi->SetName("hp_"+id);
  prt_canvasAdd(cdigi);
  prt_canvasSave("data/drawHP",0);
  auto s = prt_getPixData("m,p,v\n",2030);
  prt_writeString("digi.csv", s);
}

