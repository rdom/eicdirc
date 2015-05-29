#include "../../prttools/prttools.C"
void addcanvases(){
  fSavePath = "data/perfLL";
  const Int_t narr = 20;
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 

  TFile *f1 = TFile::Open("c_l3.root");
  TIter next1(f1->GetListOfKeys());
  TKey *key1;
  Int_t it1 = 0;
  TCanvas *carr1[narr];

  while((key1 = (TKey*)next1())) {
    TClass *cl = gROOT->GetClass(key1->GetClassName());
    if (!cl->InheritsFrom("TCanvas")) continue;
    carr1[it1] = (TCanvas*)key1->ReadObj();
    it1++;
  }

  TFile *f2 = TFile::Open("c_l0.root");
  TIter next2(f2->GetListOfKeys());
  TKey *key2;
  Int_t it2 = 0;
  TCanvas *carr2[narr];

  while ((key2 = (TKey*)next2())) {
    TClass *cl = gROOT->GetClass(key2->GetClassName());
    if (!cl->InheritsFrom("TCanvas")) continue;
    carr2[it2] = (TCanvas*)key2->ReadObj();
    it2++;
  }

  TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);


  for(Int_t i=0; i<it2; i++){
    carr1[i]->Draw();
    //leg->AddEntry( carr1[i],"3x3 full coverage","l");
    //leg->Draw();
    canvasAdd(carr1[i]);
    TIter next(carr2[i]->GetListOfPrimitives());
    TObject *obj;

    while((obj = next())){
      // if(obj->InheritsFrom("TH1F")){
      // 	TH1F *h = (TH1F*)obj;
      // 	std::cout<<"name "<< h->GetName() <<std::endl;      
      // 	h->SetLineStyle(7);
      // 	h->SetLineWidth(2);
      //   h->Draw("same");
      // }
      if(obj->InheritsFrom("TGraph")){
	TGraph *h = (TGraph*)obj;
	std::cout<<"name "<< h->GetName() <<std::endl;      
	h->SetLineColor(32);
	h->SetMarkerColor(2);
	//	h->SetLineWidth(2);
        h->Draw("same PL");
	// leg->AddEntry(h,"6.5x6.5 MCP PMTs coverage","lp");
	// leg->Draw();
      }

    }
  }
  std::cout<<"save all  " <<std::endl;
  
  canvasSave(0,1);
}
